% updated on 190220 for the combined analysis the RNA MERFISH and DNA
% tracing. The RNA MERFISH uses 16 rounds of readout hybs and all readout
% hybs are in one (Cy7) fluorescent channel.

clear all
close all
NFOV = 14; % number of fields of views
ParameterRange = 2:7;
NumSteps = 11;
load('CodeBookSubPool3_190602.mat');

%%
N = length(Codebook); % number of genes
load('AllCodes.mat')
NumCodesAll = length(GoodCodes); % number of all posible codes

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
%         FPKMs(i) = nan;
    end
end
%%
Fig1Corr = [];
Fig2Corr = [];
Fig3Corr = [];
for iii = ParameterRange
    CopyNumbers = zeros(NumCodesAll,1);
    CopyNumbers_perfect = zeros(NumCodesAll,1);
    for jj = 0:NFOV-1
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
        
        for kk = 1:NumSteps
            if exist(['results' num2str(iii) '/CC_' FOVid '_' num2str(kk) '.mat'])
                load(['results' num2str(iii) '/CC_' FOVid '_' num2str(kk) '.mat']);
                load(['results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat']);
                CopyNumbers = CopyNumbers + [CC.NumObjects]';
                CopyNumbers_perfect = CopyNumbers_perfect + [CC_perfect.NumObjects]';
            end
        end
    end

    
    CopyNumbers_Gene = CopyNumbers(1:N);
    Ind = find(FPKMs>0 & CopyNumbers_Gene>0);
%     figure(1)
%     loglog(FPKMs(Ind), CopyNumbers_Gene(Ind), 'o');
%     title(['R = ' num2str(corr(log10(FPKMs(Ind)), log10(CopyNumbers_Gene(Ind))))])
%     xlabel('FPKM');
%     ylabel('Copy number')
    
    Fig1Corr =[Fig1Corr corr(log10(FPKMs(Ind)), log10(CopyNumbers_Gene(Ind)))];

    
    CopyNumbers_perfect_Gene = CopyNumbers_perfect(1:N);
    Ind = find(FPKMs>0 & CopyNumbers_perfect_Gene>0);
%     figure(2)
%     loglog(FPKMs(Ind), CopyNumbers_perfect_Gene(Ind), 'o');
%     title(['R = ' num2str(corr(log10(FPKMs(Ind)), log10(CopyNumbers_perfect_Gene(Ind))))])
%     xlabel('FPKM');
%     ylabel('Copy number (perfect match)')
    
    Fig2Corr =[Fig2Corr corr(log10(FPKMs(Ind)), log10(CopyNumbers_perfect_Gene(Ind)))];


end
figure(4)
plot(ParameterRange,Fig1Corr)
xlabel('Parameter set');
ylabel('Correlation coefficient');
title('MHD4');
saveas(gcf, 'Correlation coefficients for different parameter sets (MHD4).jpg')

figure(5)
plot(ParameterRange,Fig2Corr)
xlabel('Parameter set');
ylabel('Correlation coefficient');
title('perfect match');
saveas(gcf, 'Correlation coefficients for different parameter sets (perfect match).jpg')


[M, I] = max(Fig1Corr);

BestParameterSet = ParameterRange(I(1));
save('BestParameterSet.mat','BestParameterSet');



