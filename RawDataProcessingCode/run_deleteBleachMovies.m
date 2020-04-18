% delete all bleach movies except the first FOV of each hyb
% updated on 190314.
clear all
close all
NFOV = 14; % number of fields of views
NumHybs = 40; % number of secondary hybs

for jj = 1:NFOV-1
    if NFOV<=10
        FOVid = num2str(jj);
    elseif NFOV>10 && NFOV<=100
        if jj<10
            FOVid = ['0' num2str(jj)];
        else
            FOVid = [num2str(jj)];
        end
    elseif NFOV>100
        if jj<10
            FOVid = ['00' num2str(jj)];
        elseif jj<100
            FOVid = ['0' num2str(jj)];
        else
            FOVid = [num2str(jj)];
        end
    end
    for i = 0:NumHybs
        if i<10
            FileName = ['sequential/Bleach_0' num2str(i) '_' FOVid '.*'];
        else
            FileName = ['sequential/Bleach_' num2str(i) '_' FOVid '.*'];
        end
        delete(FileName)
    end
end
    