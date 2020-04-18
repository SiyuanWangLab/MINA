clear all
close all
% Analyze the fine scale chromatin traces. 
% List datasets to include:
Folder{1} = 'G:\Miao\190622';
NFOV(1) = 12; % number of fields of views
Folder{2} = 'G:\Miao\190706';
NFOV(2) = 12; % number of fields of views
Folder{3} = 'G:\Miao\190715';
NFOV(3) = 12; % number of fields of views
Folder{4} = 'G:\Miao\190725';
NFOV(4) = 14; % number of fields of views

FolderIDwithFineTracing = [3, 4]; % this indicate which of all folders has the fine tracing data
ImageSize = 1536; % number of pxls
TotalNumTADs = 19;
ThresholdOfContact = 0.15; %um, distance shorter than this is considered contact

mkdir figures_TraceAnalysisByClusters_SmallScale
%% load in hi-c result
fid = fopen('DomainStartsAndEnds19_5kb_19loci_mm9.txt', 'r');
tline = fgetl(fid);
i = 0;
while ischar(tline)
    i = i+1;
    C = strsplit(tline);
    DomainStart(i) = str2num(C{2});
    DomainEnd(i) = str2num(C{3});
    DomainCenter(i) = (DomainStart(i)+DomainEnd(i))/2;
    tline = fgetl(fid);
end
fclose(fid);
if exist(['HiC_5kb.mat'])==2
    load('HiC_5kb.mat');
else
    HiC = zeros(20,20);
    ChrName = '19';
    fid = fopen(['GSM1718024_mouse1_INL_GCCAAT_L001_R.raw.txt'], 'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    while ischar(tline)
        C = strsplit(tline);
        if strcmp(C{1},ChrName) && strcmp(C{4},ChrName)
            if str2num(C{2})>DomainStart(1) && str2num(C{2})<DomainEnd(end) && ...
                    str2num(C{5})>DomainStart(1) && str2num(C{5})<DomainEnd(end)
                Bin1 = ceil((str2num(C{2})-DomainStart(1))/5000);
                Bin2 = ceil((str2num(C{5})-DomainStart(1))/5000);
                HiC(Bin1,Bin2) = HiC(Bin1,Bin2)+1;
                HiC(Bin2,Bin1) = HiC(Bin2,Bin1)+1;
            end
        end
        tline = fgetl(fid);
    end
    HiC = HiC([[1:5], [7:20]],:);
    HiC = HiC(:,[[1:5], [7:20]]);
    save('HiC_5kb.mat','HiC');
    fclose(fid);
end
figure(2)
imagesc(HiC)
colorbar
title('HiC interaction map');
axis square
ColorMap = load('RedBlue.txt');
L = size(ColorMap,1);
ColorMap = ColorMap(L:-1:1,:);
colormap(ColorMap/255);
PlotProp
savefig('HiC interaction map 5kb.fig')

%% Load in ChrList from each dataset, and remap the cell type based on the combined cell type identification
load('Clustering results final.mat')
% First plot the scd2 expression in different cell types
Scd2Exp = Matrix(:,114);
Ind = find(CellTypeID_new<8 & CellTypeID_new>0);
Scd2Exp_crop = Scd2Exp(Ind);
CellTypeID_new_crop = CellTypeID_new(Ind);
for i = 1:7
    Ind = find(CellTypeID_new_crop==i);
    Scd2ExpInGroup(i) = mean(Scd2Exp_crop(Ind));
    EnrichmentScd2ExpInGroup(i) = Scd2ExpInGroup(i)/mean(Scd2Exp_crop(find(CellTypeID_new_crop~=i)));
end
figure(1)
barh(1:7,EnrichmentScd2ExpInGroup,'LineWidth',2,'FaceColor',[0,158,115]/255,'BaseValue',1)
hold on
plot([1 1],[0 8],'-k','LineWidth',2)
hold off
xlabel('Fold enrichment of scd2 expression')
yticklabels({'Other','Endothelial cell','Hepatocyte','Macrophage','Megakaryocyte','Erythroblast','Proerythroblast'})
PlotProp
savefig(['figures_TraceAnalysisByClusters_SmallScale/Fold enrichment of scd2 expression.fig']);

Chr_all = [];
for ii = 1:length(FolderIDwithFineTracing)
    FolderID = FolderIDwithFineTracing(ii);
    % find cells in this folder (dataset)
    CellsInFolder = CellList(find([CellList.DatasetID]==FolderID));
    load([Folder{FolderID} '\FineTraceList.mat']);
    for jj = 0:NFOV(FolderID)-1
        I = zeros(ImageSize);
        I2 = zeros(ImageSize);
        % find cells in this FOV
        CellsInFOV = CellsInFolder(find([CellsInFolder.FOV] == jj));
        for i = 1:length(CellsInFOV)
            I(CellsInFOV(i).PixelList) = CellsInFOV(i).CellType;
            I2(CellsInFOV(i).PixelList) = i;
        end
        % find chromosomes in this FOV
        ChrInFOV = Chr(find([Chr.FOV] == jj));
        %find the corresponding cell type of this chromosome.
        if length(ChrInFOV)>0
            for i = 1:length(ChrInFOV)
                MeanY = ChrInFOV(i).MeanYWGA;
                MeanX = ChrInFOV(i).MeanXWGA;
                if MeanY>0 && MeanY<=ImageSize && MeanX>0 && MeanX<=ImageSize
                    ChrInFOV(i).CellType = I(MeanY,MeanX);
                    if I2(MeanY,MeanX)>0
                        ChrInFOV(i).RNACopyNumber = CellsInFOV(I2(MeanY,MeanX)).RNACopyNumber;
                    else
                        ChrInFOV(i).RNACopyNumber = [];
                    end
                else
                    ChrInFOV(i).CellType = 0;
                    ChrInFOV(i).RNACopyNumber = [];
                end
                ChrInFOV(i).DatasetID = FolderID;
            end
            Chr_all = [Chr_all ChrInFOV];
        end
    end
end
Chr = Chr_all(find([Chr_all.CellType]>0));
save('5kb_res_trace.mat','Chr');
display(['Total number of cells analyzed: ' num2str(length(CellList))])
display(['Total number of chromosomes analyzed: ' num2str(length(Chr))])
% compare overall mean spatial distance with Hi-C
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(Chr)
            if Chr(k).r(i) == 1 && Chr(k).r(j) == 1
                DisList = [DisList ((Chr(k).x(i)-Chr(k).x(j))^2+(Chr(k).y(i)-Chr(k).y(j))^2+(Chr(k).z(i)-Chr(k).z(j))^2)^0.5];
            end
        end
        Mean(i,j) = mean(DisList);
        Std(i,j) = std(DisList);
        SEM(i,j) = std(DisList)/(length(DisList))^0.5;
        NofData(i,j) = length(DisList);
        DisListAll{i,j} = DisList;
    end
end
figure(3)
imagesc(Mean)
colorbar
RedBlue
caxis([0 0.6])
PlotProp
axis square
savefig(['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance.fig']);

%% Plot the difference between the mean spatial distance matrices of Hepatocyte and other cell types.
Mean_Hep = [];
Mean_other = [];

ChrInGroup = Chr(find([Chr.CellType]==3));
N_hep = length(ChrInGroup);
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(ChrInGroup)
            if ChrInGroup(k).r(i) == 1 && ChrInGroup(k).r(j) == 1
                DisList = [DisList ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
            end
        end
        Mean_Hep(i,j) = mean(DisList);
    end
end

ChrInGroup = Chr(find([Chr.CellType]>=1 & [Chr.CellType]<=7 & [Chr.CellType]~=3));
N_nonhep = length(ChrInGroup);
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(ChrInGroup)
            if ChrInGroup(k).r(i) == 1 && ChrInGroup(k).r(j) == 1
                DisList = [DisList ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
            end
        end
        Mean_other(i,j) = mean(DisList);
    end
end   
figure(4)
imagesc(Mean_Hep)
colorbar
RedBlue
caxis([0 0.6])
PlotProp
axis square
title(['N = ' num2str(N_hep)])
savefig(['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance Hep.fig']);
figure(5)
imagesc(Mean_other)
colorbar
RedBlue
caxis([0 0.6])
PlotProp
axis square
title(['N = ' num2str(N_nonhep)])
savefig(['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance All except Hep.fig']);
figure(6)
imagesc(Mean_Hep-Mean_other)
colorbar
RedBlue
PlotProp
axis square
savefig(['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance Hep minus All others.fig']);
figure(7)
plot(1:TotalNumTADs,Mean_Hep(19,:)-Mean_other(19,:),'.-','MarkerSize',20,'LineWidth',2)
xlabel('ID number of 5-kb regions');
ylabel('Change of distance (um)')
PlotProp
savefig(['figures_TraceAnalysisByClusters_SmallScale/Change of distance to scd2 promoter in heptocytes.fig']);

% Plot the enrichment of contact probability between different loci with locus 19, in Hepatocyte/all other cell types.
ChrInGroup = Chr(find([Chr.CellType]==3));
for i = 1:TotalNumTADs
    j = 19;
    DisList = [];
    for k = 1:length(ChrInGroup)
        if ChrInGroup(k).r(i) == 1 && ChrInGroup(k).r(j) == 1
            DisList = [DisList ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
        end
    end
    ContactProb_Hep(i) = length(find(DisList<ThresholdOfContact))/length(DisList);
end

ChrInGroup = Chr(find([Chr.CellType]>=1 & [Chr.CellType]<=7 & [Chr.CellType]~=3));
for i = 1:TotalNumTADs
    j = 19;
    DisList = [];
    for k = 1:length(ChrInGroup)
        if ChrInGroup(k).r(i) == 1 && ChrInGroup(k).r(j) == 1
            DisList = [DisList ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
        end
    end
    ContactProb_other(i) = length(find(DisList<ThresholdOfContact))/length(DisList);
end  
figure(8)
plot(1:TotalNumTADs,ContactProb_Hep-ContactProb_other,'.-','MarkerSize',20,'LineWidth',2)
hold on
plot([0 20],[0 0],'--k','LineWidth',2)
hold off
xlabel('ID number of 5-kb regions');
ylabel('Change of contact prob')
PlotProp
savefig(['figures_TraceAnalysisByClusters_SmallScale/Change of contact prob with scd2 promoter in heptocytes.fig']);

figure(80)
plot(1:TotalNumTADs,ContactProb_Hep,'.-','MarkerSize',20,'LineWidth',2)
hold on
plot(1:TotalNumTADs,ContactProb_other,'.-','MarkerSize',20,'LineWidth',2)
plot([0 20],[0 0],'--k','LineWidth',2)
hold off
xlabel('ID number of 5-kb regions');
ylabel('Contact prob')
PlotProp
savefig(['figures_TraceAnalysisByClusters_SmallScale/Contact prob with scd2 promoter in heptocytes and non-hep.fig']);


% Ask if multiple contacts are coorporative with each other
ChrInGroup = Chr(find([Chr.CellType]==3));
i1 = 16;
i2 = 10;
j = 19;
Dis1 = [];
Dis2 = [];
for k = 1:length(ChrInGroup)
    if ChrInGroup(k).r(i1) == 1 && ChrInGroup(k).r(i2) == 1 &&ChrInGroup(k).r(j) == 1
        Dis1 = [Dis1 ((ChrInGroup(k).x(i1)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i1)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i1)-ChrInGroup(k).z(j))^2)^0.5];
        Dis2 = [Dis2 ((ChrInGroup(k).x(i2)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i2)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i2)-ChrInGroup(k).z(j))^2)^0.5];
    end
end
% contact 1 prob
p1 = length(find(Dis1<ThresholdOfContact))/length(Dis1);
p1_error = sqrt(p1*(1-p1)/length(Dis1));
% contact 1 probability given contact 2
p1g2 = length(find(Dis1<ThresholdOfContact & Dis2<ThresholdOfContact))/length(find(Dis2<ThresholdOfContact));
p1g2_error = sqrt(p1g2*(1-p1g2)/length(find(Dis2<ThresholdOfContact)));
% contact 1 probability given no contact 2
p1gno2 = length(find(Dis1<ThresholdOfContact & Dis2>=ThresholdOfContact))/length(find(Dis2>=ThresholdOfContact));
p1gno2_error = sqrt(p1gno2*(1-p1gno2)/length(find(Dis2>=ThresholdOfContact)));
figure(9)
barh([1 2 3], [p1 p1g2 p1gno2],'LineWidth',2)
hold on
% contact 2 prob
p2 = length(find(Dis2<ThresholdOfContact))/length(Dis2);
p2_error = sqrt(p2*(1-p2)/length(Dis2));
% contact 2 probability given contact 1
p2g1 = length(find(Dis1<ThresholdOfContact & Dis2<ThresholdOfContact))/length(find(Dis1<ThresholdOfContact));
p2g1_error = sqrt(p2g1*(1-p2g1)/length(find(Dis1<ThresholdOfContact)));
% contact 2 probability given no contact 1
p2gno1 = length(find(Dis2<ThresholdOfContact & Dis1>=ThresholdOfContact))/length(find(Dis1>=ThresholdOfContact));
p2gno1_error = sqrt(p2gno1*(1-p2gno1)/length(find(Dis1>=ThresholdOfContact)));
barh([4 5 6], [p2 p2g1 p2gno1],'LineWidth',2)
errorbar([p1 p1g2 p1gno2 p2 p2g1 p2gno1],[1:6],[p1_error p1g2_error p1gno2_error p2_error p2g1_error p2gno1_error],'.k','horizontal','LineWidth',2,'CapSize',10)
yticks(1:6)
yticklabels({'C1','C1 | C2','C1 | no-C2','C2','C2 | C1','C2 | no-C1'});
xlabel('Probability')
PlotProp
hold off
savefig(['figures_TraceAnalysisByClusters_SmallScale/Conditional contact probabilities.fig']);

% Ask if contacts are correlated with Scd2 copy number in the single hepatocyte cells.
% First analyze C1 with copy number
ChrInGroup = Chr(find([Chr.CellType]==3));
i = 16;
j = 19;
Dis = [];
Scd2CopyNumbers = [];
for k = 1:length(ChrInGroup)
    if ChrInGroup(k).r(i) == 1 &&ChrInGroup(k).r(j) == 1
        Dis = [Dis ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
        Scd2CopyNumbers = [Scd2CopyNumbers ChrInGroup(k).RNACopyNumber(114)];
    end
end
figure(10)
Scd2CopyNumbers_C1 = Scd2CopyNumbers(find(Dis<ThresholdOfContact));
Scd2CopyNumbers_noC1 = Scd2CopyNumbers(find(Dis>=ThresholdOfContact));
Mean1 = mean(Scd2CopyNumbers_C1);
Mean2 = mean(Scd2CopyNumbers_noC1);
SEM1 = std(Scd2CopyNumbers_C1)/sqrt(length(Scd2CopyNumbers_C1));
SEM2 = std(Scd2CopyNumbers_noC1)/sqrt(length(Scd2CopyNumbers_noC1));
bar([1 2],[Mean1 Mean2],'LineWidth',2);
hold on
errorbar([1 2],[Mean1 Mean2],[SEM1 SEM2],'.k','LineWidth',2,'CapSize',10)
p1 = ranksum(Scd2CopyNumbers_C1, Scd2CopyNumbers_noC1); % Wilcoxon rank sum test

% Second analyze C2 with copy number
ChrInGroup = Chr(find([Chr.CellType]==3));
i = 10;
j = 19;
Dis = [];
Scd2CopyNumbers = [];
for k = 1:length(ChrInGroup)
    if ChrInGroup(k).r(i) == 1 &&ChrInGroup(k).r(j) == 1
        Dis = [Dis ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
        Scd2CopyNumbers = [Scd2CopyNumbers ChrInGroup(k).RNACopyNumber(114)];
    end
end
figure(10)
Scd2CopyNumbers_C2 = Scd2CopyNumbers(find(Dis<ThresholdOfContact));
Scd2CopyNumbers_noC2 = Scd2CopyNumbers(find(Dis>=ThresholdOfContact));
Mean1 = mean(Scd2CopyNumbers_C2);
Mean2 = mean(Scd2CopyNumbers_noC2);
SEM1 = std(Scd2CopyNumbers_C2)/sqrt(length(Scd2CopyNumbers_C2));
SEM2 = std(Scd2CopyNumbers_noC2)/sqrt(length(Scd2CopyNumbers_noC2));
bar([3 4],[Mean1 Mean2],'LineWidth',2);
errorbar([3 4],[Mean1 Mean2],[SEM1 SEM2],'.k','LineWidth',2,'CapSize',10)
p2 = ranksum(Scd2CopyNumbers_C2, Scd2CopyNumbers_noC2); % Wilcoxon rank sum test
title(['p1 = ' num2str(p1) ', p2 = ' num2str(p2)])
ylabel('scd2 copy number');
xticks([1:4]);
xticklabels({'C1','no C1','C2','no C2'})
PlotProp
hold off
savefig(['figures_TraceAnalysisByClusters_SmallScale/scd2 expression in C vs no-C.fig']);
