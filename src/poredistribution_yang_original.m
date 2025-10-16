function [C0, C1, Re] = poredistribution_yang_original(C)

a1=1;
a2=size(C,1); % pixel number in the x-direction
b1=1;
b2=size(C,2); % pixel number in the y-directioin
c1=1;
c2=size(C,3); % total number of the images in a series

Re(100)=0; % allocated stream for data storage

C0=0*C;
C1=0*C;
for i=a1:a2 % begin to find the critic radius for each unit-valued pixel, and store the value of critic radius in the matrix “C1”
    for j=b1:b2
        for k=c1:c2
            if C(i,j,k)~=0
                mark=1;
                l=0;
                while mark>0
                    C0(i,j,k)=l+0.5;
                    l=l+1;
                    if ((i-l)<=(a1-1)||(j-l)<=(b1-1)||(k-l)<=(c1-1)||(i+l)>=(a2+1)||(j+l)>=(b2+1)||(k+l)>=(c2+1))
                        mark=0;
                    end
                    if (mark~=0)
                        for aa=(i-l):(i+l)
                            for bb=(j-l):(j+l)
                                for cc=(k-l):(k+l)
                                    if sqrt((aa-i)^2+(bb-j)^2+(cc-k)^2)<=l
                                        if C(aa,bb,cc)==0
                                            mark=0;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end % finish the critic radius-finding procedure

dpm=max(max(max(C0)))-0.5;
for dp=dpm:-1:0 % begin to locate the surrounding pixels
    for i=a1:a2
        for j=b1:b2
            for k=c1:c2
                if C0(i,j,k)==dp+0.5
                    for aa=(i-dp):(i+dp)
                        for bb=(j-dp):(j+dp)
                            for cc=(k-dp):(k+dp)
                                if sqrt((aa-i)^2+(bb-j)^2+(cc-k)^2)<=dp
                                    if C1(aa,bb,cc)==0
                                        C1(aa,bb,cc)=dp+1;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end % finish the surrounding pixel-finding procedure

for dp=0:dpm % begin to store the pixel number at different critic radius into a stream named “Re”
    for i=a1:a2
        for j=b1:b2
            for k=c1:c2
                if C1(i,j,k)==dp+1
                    Re(dp+1)= Re(dp+1)+1;
                    Re=Re';
                end
            end
        end
    end
end % complete the procedure
Re = Re(:);
% RGB_label = label2rgb(C1(:,:,round((c1+c2)/2)),'jet'); % paint pore area with different colors according to pore size
% figure, imshow(RGB_label) % show the colored pore configuration

end