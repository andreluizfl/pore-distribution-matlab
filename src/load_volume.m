function C = load_volume(imgDir, extType, volumetricSize, binThreshold, useParallel)
% LOAD_VOLUME  Read an image stack and produce a binary 3D volume.
%
%   C = LOAD_VOLUME(imgDir)
%   C = LOAD_VOLUME(imgDir, extType, volumetricSize, binThreshold, useParallel)
%
% INPUTS:
%   imgDir         - (char) path to the folder containing the image stack.
%                     A path separator at the end is optional.
%   extType        - (char, optional) image extension to search for,
%                     e.g. '.bmp' (default), '.tif', '.png'.
%   volumetricSize - (optional) defines the sub-volume to read. Accepts:
%                     [] (default)         -> full image extents and all slices
%                     string or char 'fit' -> automatically fits a centered cubic volume 
%                     whose side length equals the smallest image dimension 
%                     scalar N             -> NxNxN volume starting at (1,1,1)
%                     3x2 matrix [x1 x2; y1 y2; z1 z2] -> explicit ranges
%   binThreshold   - (numeric, optional) binarization threshold.
%                     -1 (default) => compute global Otsu threshold from the
%                     entire loaded volume (intensity normalized to [0,1]).
%                     Otherwise use the provided numeric threshold (in [0,1]).
%   useParallel     - (logical, optional) request parallel processing.
%                     If omitted or empty, true is assumed and the function
%                     will attempt to use a parpool if available. If no
%                     parallel pool can be created, the function falls back
%                     to serial processing.
%
% OUTPUT:
%   C - logical 3D array (rows x cols x slices) containing the binary
%       pore mask (true => pore / foreground).
%
% NOTES & IMPLEMENTATION DETAILS:
%   - All input images must share identical dimensions. The function checks
%     and throws an error if a mismatch is detected.
%   - Images read as RGB will be converted to grayscale via rgb2gray.
%   - Integer images (uint8, uint16) are normalized to [0,1] prior to
%     Otsu thresholding. Floating point images are left as-is.
%   - When binThreshold == -1 the function computes a global Otsu level
%     using graythresh on the flattened normalized intensities of the
%     loaded volume.
%   - The function uses parfor to parallelize both image reading and
%     per-slice binarization when a parallel pool is available and
%     useParallel==true.
%
% EXAMPLE:
%   C = load_volume('data/images', '.tif', [], -1, true);
%
% Author: adapted and modularized
% Date:   (no date hard-coded here)
%

% -------------------- Defaults and input normalization --------------------
if nargin < 2 || isempty(extType)
    extType = '.*';
end
if nargin < 3
    volumetricSize = [];
end
if nargin < 4 || isempty(binThreshold)
    binThreshold = -1;
end
if nargin < 5 || isempty(useParallel)
    useParallel = true;
end

% Ensure imgDir ends with filesep for safe concatenation
if isempty(imgDir)
    error('imgDir must be a non-empty directory path.');
end
if imgDir(end) ~= filesep
    imgDir = [imgDir filesep];
end

% -------------------- Find and sort files --------------------------------
files = dir([imgDir '*' extType]);
files = files(~[files.isdir]);  % remove directories

if isempty(files)
    error('No image files found in the specified directory: %s', imgDir);
end

% Numeric-aware sorting: extract trailing number groups and sort by last group
file_numbers = zeros(1, numel(files));
for ff = 1:numel(files)
    nums = regexp(files(ff).name, '\d+', 'match');
    if isempty(nums)
        file_numbers(ff) = 0;
    else
        file_numbers(ff) = str2double(nums{end});
    end
end
[~, sort_idx] = sort(file_numbers);
files = files(sort_idx);
num_images = numel(files);

% -------------------- Validate consistent image geometry -----------------
first_info = imfinfo([imgDir files(1).name]);
orig_rows = first_info.Height;
orig_cols = first_info.Width;

for ff = 2:numel(files)
    info_ff = imfinfo([imgDir files(ff).name]);
    if info_ff.Height ~= orig_rows || info_ff.Width ~= orig_cols
        error('All images must have identical dimensions. File "%s" differs.', files(ff).name);
    end
end

% -------------------- Determine volumetric ranges ------------------------
if isempty(volumetricSize)
    rangeX = [1 orig_cols];
    rangeY = [1 orig_rows];
    rangeZ = [1 num_images];
elseif isstring(volumetricSize) || ischar(volumetricSize)
    if  string(volumetricSize) == "fit"
        N = min([orig_rows,orig_cols,num_images]);
        diff_x = round((orig_cols-N)/2);
        diff_y = round((orig_rows-N)/2);
        diff_z = round((num_images-N)/2);
        rangeX = [diff_x+1 diff_x+N];
        rangeY = [diff_y+1 diff_y+N];
        rangeZ = [diff_z+1 diff_z+N];
    else
        error("Invalid volumetricSize. Provide a invalid name. Must be 'fit'");
    end
elseif isscalar(volumetricSize)
    N = round(volumetricSize);
    if N <= 0, error('Scalar volumetricSize must be positive.'); end
    rangeX = [1 N]; rangeY = [1 N]; rangeZ = [1 N];
elseif ismatrix(volumetricSize) && all(size(volumetricSize) == [3 2])
    rangeX = volumetricSize(1, :);
    rangeY = volumetricSize(2, :);
    rangeZ = volumetricSize(3, :);
else
    error('Invalid volumetricSize. Provide [] | scalar | 3x2 matrix.');
end

% Clamp user ranges to available image dimensions / slice count
rangeX(2) = min(rangeX(2), orig_cols);
rangeY(2) = min(rangeY(2), orig_rows);
rangeZ(2) = min(rangeZ(2), num_images);

rows_range = rangeY(1):rangeY(2);
cols_range = rangeX(1):rangeX(2);
imgs_range = rangeZ(1):rangeZ(2);

rows = numel(rows_range);
cols = numel(cols_range);
num_imgs = numel(imgs_range);

fprintf('Using volume range: X=[%d %d], Y=[%d %d], Z=[%d %d]\n', ...
    rangeX(1), rangeX(2), rangeY(1), rangeY(2), rangeZ(1), rangeZ(2));

% -------------------- Parallel setup (best-effort) -----------------------
if useParallel
    try
        % try to ensure a pool exists; if it fails, fall back to serial
        hasDCT = license('test', 'Distrib_Computing_Toolbox');
        if hasDCT
            delete(gcp("nocreate"));
            pool = parpool("Threads");
            if isempty(pool)
                try
                    parpool('local'); %#ok<*PFBNS>
                catch
                    % fallback silently to serial if pool cannot be created
                end
            end
            % pool = gcp('nocreate');
            if isempty(pool)
                useParallel = false;
            else
                useParallel = true;
            end
        else
            useParallel = false;
        end
    catch
        useParallel = false;
    end
end

% -------------------- Load images into a single-precision volume ---------
% We'll normalize integer types to [0,1] for consistent Otsu behavior.
grayVolume = zeros(rows, cols, num_imgs, 'single');

if useParallel
    % Use parfor to parallelize reading and normalization per slice
    parfor ii = 1:num_imgs
        idxFile = imgs_range(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rows_range, cols_range);
        origClass = class(I);
        I_single = single(I);
        switch origClass
            case 'uint8'
                I_single = I_single / 255;
            case 'uint16'
                I_single = I_single / 65535;
            otherwise
                % assume already in [0,1] or floating-point; leave as-is
        end
        grayVolume(:, :, ii) = I_single;
    end
else
    for ii = 1:num_imgs
        idxFile = imgs_range(ii);
        I = imread([imgDir files(idxFile).name]);
        if ndims(I) == 3
            I = rgb2gray(I);
        end
        I = I(rows_range, cols_range);
        I_single = single(I);
        % This normalization was writen for safety, but its kind redundant
        % with graythresh and imbinarize internal operations
        % origClass = class(I);
        % switch origClass
        %     case 'uint8'
        %         I_single = I_single / 255;
        %     case 'uint16'
        %         I_single = I_single / 65535;
        %     otherwise
        % end
        grayVolume(:, :, ii) = I_single;
    end
end

% -------------------- Compute global Otsu threshold if requested ---------
if binThreshold == -1
    minVal = min(grayVolume(:));
    maxVal = max(grayVolume(:));
    if maxVal > minVal
        normalizedVol = (grayVolume - minVal) / (maxVal - minVal);
    else
        normalizedVol = grayVolume;
    end
    globalLevel = graythresh(normalizedVol(:));
else
    globalLevel = binThreshold;
end

% -------------------- Binarize the volume -------------------------------
C = false(rows, cols, num_imgs);
if useParallel
    parfor ii = 1:num_imgs
        C(:, :, ii) = imbinarize(grayVolume(:, :, ii), globalLevel);
    end
else
    for ii = 1:num_imgs
        C(:, :, ii) = imbinarize(grayVolume(:, :, ii), globalLevel);
    end
end

end
