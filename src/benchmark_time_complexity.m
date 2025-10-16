clear all; clc;

dir_src = '..\data\CT_01';
extension = '.bmp';

C = load_volume(dir_src,extension);

sizes = 10:10:100;
repeats = 5;

ts = zeros([length(sizes),3]);
check_res = false([length(sizes),6]);

i = 1;
for N = sizes
    C_sub = remap_volume(C, N, 'top-left-front');
    t_ori=0;t_opt=0;t_optp=0;
    for r = 1:repeats
        tic;
        [C0_ori, C1_ori, Re_ori] = poredistribution_yang_original(C_sub);
        t_ori = t_ori+toc;
        
        tic;
        [C0_opt, C1_opt, Re_opt] = poredistribution_yang_optimized(C_sub,false);
        t_opt = t_opt+toc;
    
        tic;
        [C0_optp, C1_optp, Re_optp] = poredistribution_yang_optimized(C_sub,true);
        t_optp = t_optp+toc;
                
        check_res(i,:) = [isequal(C0_ori,C0_opt), isequal(C1_ori,C1_opt), isequal(Re_ori,Re_opt), isequal(C0_ori,C0_optp), isequal(C1_ori,C1_optp), isequal(Re_ori,Re_optp)];
        
    end
    ts(i,:) = [t_ori, t_opt, t_optp]/repeats;
    i=i+1;
end

A = cat(2,sizes', ts);
T = array2table(A);
T.Properties.VariableNames(1:4) = {'n','original','opt_seq','opt_par'};
writetable(T,'..\results\ts.csv')

fprintf("Checking if all results are identical to the original: %s\n",string(all(check_res(:))));

figure; grid on; hold on;
plot(ts(:,1),'-r'); 
plot(ts(:,2),'-b'); 
plot(ts(:,3),'-g');
title('Time Complexity between 3 algorithms');
ylabel('Time (s)');
xlabel('N');
legend('Original','Optimized Sequential','Optimized Parallel');
saveas(gcf,'..\results\time_complexity.png');

figure; grid on; hold on;
plot(ts(:,2),'-b'); 
plot(ts(:,3),'-g');
title('Time Complexity between optimized  algorithms versions');
ylabel('Time (s)');
xlabel('N');
legend('Optimized Sequential','Optimized Parallel');
saveas(gcf,'..\results\time_complexity_opt.png');

slice = sizes(end)/2;
figure;
subplot(2,2,1); imshow(C(:,:,slice)); title('C Image');
subplot(2,2,2); plot(Re_ori); title('Pore Distribution');
subplot(2,2,3); imshow(C0_ori(:,:,slice), []); title('C0');
subplot(2,2,4); imshow(C1_ori(:,:,slice), []); title('C1');
saveas(gcf,'..\results\results_original_alg.png');

figure;
subplot(2,2,1); imshow(C(:,:,slice)); title('C Image');
subplot(2,2,2); plot(Re_opt); title('Pore Distribution');
subplot(2,2,3); imshow(C0_opt(:,:,slice), []); title('C0');
subplot(2,2,4); imshow(C1_opt(:,:,slice), []); title('C1');
saveas(gcf,'..\results\results_optimized_alg.png');
