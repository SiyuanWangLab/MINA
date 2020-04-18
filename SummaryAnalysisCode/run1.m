clear all
close all
% Combine data from mutliple datasets and analyze together.

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

%% First analyze the MERFISH results and determine the best MERFISH parameter sets to use.
% Load the bulk RNA seq results.
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

% Analyze MERFISH data
for FolderID = 1:length(Folder)
    Fig1Corr{FolderID} = []; % this is the MHD4 copy number correlation with FPKMs
    Fig2Corr{FolderID} = []; % this is the Perfect Match copy number correlation with FPKMs
    ToalCopyNumbers{FolderID} = []; % total Perfect Match RNA copy numbers for each parameter set
    FoldDifferenceSlope{FolderID} = []; % slope of Perfect Match copy number fold difference over FPKM fold difference
    for iii = ParameterRange
        CopyNumbers = zeros(NumCodesAll,1);
        CopyNumbers_perfect = zeros(NumCodesAll,1);
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
        CopyNumbers_Gene = CopyNumbers(1:N);
        Ind = find(FPKMs>0 & CopyNumbers_Gene>0);
        Fig1Corr{FolderID} =[Fig1Corr{FolderID} corr(log10(FPKMs(Ind)), log10(CopyNumbers_Gene(Ind)))];
        CopyNumbers_perfect_Gene = CopyNumbers_perfect(1:N);
        Ind = find(FPKMs>0 & CopyNumbers_perfect_Gene>0);
        Fig2Corr{FolderID} =[Fig2Corr{FolderID} corr(log10(FPKMs(Ind)), log10(CopyNumbers_perfect_Gene(Ind)))];
        ToalCopyNumbers{FolderID} = [ToalCopyNumbers{FolderID} sum(CopyNumbers_perfect_Gene)];
        b = regress(log10(CopyNumbers_perfect_Gene(Ind)),[ones(size(FPKMs(Ind))) log10(FPKMs(Ind))]);
        FoldDifferenceSlope{FolderID} = [FoldDifferenceSlope{FolderID} b(2)];    
    end
    figure(1)
    plot(ParameterRange, Fig1Corr{FolderID});
    title('MHD4 correlations')
    ylabel('Correlation coef')
    xlabel('Pameter set')
    legend;
    hold on
    figure(2)
    plot(ParameterRange, Fig2Corr{FolderID});
    title('Perfect Match correlations')
    ylabel('Correlation coef')
    xlabel('Pameter set')
    legend;
    hold on
    figure(3)
    plot(ParameterRange, ToalCopyNumbers{FolderID});
    title('Total copy numbers')
    ylabel('Total perfect match copy number')
    xlabel('Pameter set')
    legend;
    hold on
    figure(4)
    plot(ParameterRange, FoldDifferenceSlope{FolderID});
    title('Fold Difference Slopes')
    ylabel('Fold Difference Slope')
    xlabel('Pameter set')
    legend;
    hold on
end
figure(1)
hold off
savefig('MHD4 correlations.fig')
figure(2)
hold off
savefig('Perfect Match correlations.fig')
figure(3)
hold off
savefig('Total copy numbers.fig')
figure(4)
hold off
savefig('Fold Difference Slopes.fig')
% define the best parameter set for each dataset
for FolderID = 1:length(Folder)
    [M, I] = max(Fig1Corr{FolderID});
    BestParameterSet(FolderID) = ParameterRange(I(1));
end
save('BestParameterSets.mat','BestParameterSet');


