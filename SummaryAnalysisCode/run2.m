clear all
close all
% plot the combined FPKM versus perfect match copy number

% List datasets to include:
Folder{1} = 'G:\Miao\190622';
NFOV(1) = 12; % number of fields of views
Folder{2} = 'G:\Miao\190706';
NFOV(2) = 12; % number of fields of views
Folder{3} = 'G:\Miao\190715';
NFOV(3) = 12; % number of fields of views
Folder{4} = 'G:\Miao\190725';
NFOV(4) = 14; % number of fields of views

ParameterRange = 2:7; % range of MERFISH parameter sets
NumSteps = 11; % number of z steps (z sections) analyzed in MERFISH
load('CodeBookSubPool3_190602.mat');
load('BestParameterSets.mat');

%% Load the bulk RNA seq results.
N = length(Codebook); % number of genes
load('AllCodes.mat')
NumCodesAll = length(GoodCodes); % number of all posible codes
if exist('FPKMs.mat')
    load('FPKMs.mat');
    display('Load existing FPKMs mat file.')
else
    display('Load from genes.fpkm_tracking.')
    FPKMs = zeros(N,1);
    fid = fopen(['genes.fpkm_tracking'], 'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    i = 0;
    while ischar(tline)
        i = i+1;
        C = strsplit(tline);
        Gene(i).Name = str2num(C{4}(8:end));
        Gene(i).ShortName = C{5};
        Gene(i).FPKM = str2num(C{10});
        tline = fgetl(fid);
    end
    fclose(fid);
    display('Gene names and FPKM loaded.');
    for i = 1:N
        flag = 0;
        for j = 1:length(Gene)
            if strcmpi(Codebook(i).GeneShortName,Gene(j).ShortName)
                FPKMs(i) = FPKMs(i) + Gene(j).FPKM;
                flag = flag+1;
            end
        end
        if flag == 0
            display(['No FPKM found for ' Codebook(i).GeneShortName]);
            FPKMs(i) = nan;
        elseif flag > 1
            display(['More than 1 FPKM found for ' Codebook(i).GeneShortName]);
        end
    end
    save('FPKMs.mat','FPKMs');
end

%% Combine the best parameter sets from each dataset
CopyNumbers = zeros(NumCodesAll,1);
CopyNumbers_perfect = zeros(NumCodesAll,1);
for FolderID = 1:length(Folder)
    iii = BestParameterSet(FolderID);
    for jj = 0:NFOV(FolderID)-1
        if NFOV(FolderID)<=10
            FOVid = num2str(jj);
        elseif NFOV(FolderID)>10 && NFOV(FolderID)<=100
            if jj<10
                FOVid = ['0' num2str(jj)];
            else
                FOVid = [num2str(jj)];
            end
        elseif NFOV(FolderID)>100
            if jj<10
                FOVid = ['00' num2str(jj)];
            elseif jj<100
                FOVid = ['0' num2str(jj)];
            else
                FOVid = [num2str(jj)];
            end
        end
        for kk = 1:NumSteps
            if exist([Folder{FolderID} '/results' num2str(iii) '/CC_' FOVid '_' num2str(kk) '.mat'])
                load([Folder{FolderID} '/results' num2str(iii) '/CC_' FOVid '_' num2str(kk) '.mat']);
                load([Folder{FolderID} '/results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat']);
                CopyNumbers = CopyNumbers + [CC.NumObjects]';
                CopyNumbers_perfect = CopyNumbers_perfect + [CC_perfect.NumObjects]';
            end
        end
    end
end
figure(1)
CopyNumbers_Gene = CopyNumbers(1:N);
Ind = find(FPKMs>0 & CopyNumbers_Gene>0);
loglog(FPKMs(Ind), CopyNumbers_Gene(Ind), 'o');
title(['MHD4, Count = ' num2str(sum(CopyNumbers_Gene)) ', R = ' num2str(corr(log10(FPKMs(Ind)), log10(CopyNumbers_Gene(Ind))))])
xlabel('FPKM');
ylabel('Copy number')
savefig('Correlation coefficient for best parameter set (MHD4).fig')
figure(2)
CopyNumbers_perfect_Gene = CopyNumbers_perfect(1:N);
Ind = find(FPKMs>0 & CopyNumbers_perfect_Gene>0);
loglog(FPKMs(Ind), CopyNumbers_perfect_Gene(Ind), 'ko','MarkerFaceColor', 'k');
title(['Count = ' num2str(sum(CopyNumbers_perfect_Gene)) ', R = ' num2str(corr(log10(FPKMs(Ind)), log10(CopyNumbers_perfect_Gene(Ind))))])
xlabel('FPKM');
ylabel('Copy number (perfect match)')
savefig('Correlation coefficient for best parameter set (perfect match).fig')




