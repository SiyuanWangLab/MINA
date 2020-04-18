clear all
close all
% Analyze the large scaling tracing, LAD and NAD by cell types

% List datasets to include:
Folder{1} = 'G:\Miao\190622';
NFOV(1) = 12; % number of fields of views
Folder{2} = 'G:\Miao\190706';
NFOV(2) = 12; % number of fields of views
Folder{3} = 'G:\Miao\190715';
NFOV(3) = 12; % number of fields of views
Folder{4} = 'G:\Miao\190725';
NFOV(4) = 14; % number of fields of views

ImageSize = 1536; % number of pxls
ChrName = 'chr19';
TotalNumTADs = 50;
DomainStartsAndEnds = 'DomainStartsAndEnds19.txt';
load('ChosenDomains19.mat');
DomainsToExclude = []; % list the domains to exclude, e.g. DomainsToExclude = [15, 2];
NormalizeTowardsMean = 0; % 1 means using the running average method to calcualte the expected distance between genomic loci
MaxNumMissedTAD = TotalNumTADs; % maximum number of missed TADs in polarization analysis

%% Load fetal liver hi-C data
if exist(['HiC.mat'])==2
    load('HiC.mat');
    NumBins = size(HiC,1);
    BinCenters = 40000*[1:NumBins] - 20000;
else
    % first define an empty HiC matrix based on the mESC Hi-C data
    HiC = load(mES_HiC);
    NumBins = size(HiC,1);
    BinCenters = 40000*[1:NumBins] - 20000;
    HiC = HiC*0;
    % load in fetal liver hic data
    fid = fopen(['GSM1718024_mouse1_INL_GCCAAT_L001_R.raw.txt'], 'r');
    tline = fgetl(fid);
    tline = fgetl(fid);
    while ischar(tline)
        C = strsplit(tline);
        if strcmp(C{1},ChrName) && strcmp(C{4},ChrName)
            Coord1 = ceil(str2num(C{2})/40000);
            Coord2 = ceil(str2num(C{5})/40000);
            HiC(Coord1,Coord2) = HiC(Coord1,Coord2)+1;
            HiC(Coord2,Coord1) = HiC(Coord2,Coord1)+1;
        end
        tline = fgetl(fid);
    end
    save('HiC.mat','HiC');    
end
% Load domain starts and ends
fid = fopen(DomainStartsAndEnds, 'r');
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
DomainStart = DomainStart(ChosenDomains);
DomainEnd = DomainEnd(ChosenDomains);
DomainCenter = DomainCenter(ChosenDomains);
for i = 1:length(DomainStart)
    BinIDs = [];
    for j = 1:NumBins
        if BinCenters(j)>DomainStart(i) && BinCenters(j)<DomainEnd(i)
            BinIDs = [BinIDs j];
        end
    end
    ChozenBinID(i,:) = [min(BinIDs) max(BinIDs)];
end
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        HiC_chozen(i, j) = sum(sum(HiC(ChozenBinID(i,1):ChozenBinID(i,2),ChozenBinID(j,1):ChozenBinID(j,2))));
        HiC_chozen(i, j) = HiC_chozen(i, j)/(ChozenBinID(i,2)-ChozenBinID(i,1)+1)/(ChozenBinID(j,2)-ChozenBinID(j,1)+1);
    end
end
figure(1)
for i = 1:length(HiC_chozen)
    HiC_chozen(i, i) = 0;
end
for i = 1:length(HiC_chozen)
    HiC_chozen(i, i) = max(max(HiC_chozen));
end
imagesc(HiC_chozen)
colorbar
title('HiC interaction map');
axis square
ColorMap = load('RedBlue.txt');
L = size(ColorMap,1);
ColorMap = ColorMap(L:-1:1,:);
ColorMap = [ColorMap; repmat(ColorMap(end,:), 3200,1)];
ColorMap = ColorMap(1:11:end,:);

colormap(ColorMap/255);
PlotProp
savefig('HiC interaction map.fig')

%% load in the known canonical gene table from USCS to calcualte the gene density for the domains
GeneCounts = zeros(size(DomainCenter));
fid = fopen(['knownCanonical.txt'], 'r');
tline = fgetl(fid);
while ischar(tline)
    C = strsplit(tline);
    if strcmp(C{1},ChrName)
        GenePosition = mean([str2num(C{2}) str2num(C{3})]);
        for i = 1:length(DomainCenter)
            if GenePosition>DomainStart(i) && GenePosition<DomainEnd(i)
                GeneCounts(i) = GeneCounts(i)+1;
            end
        end
    end
    tline = fgetl(fid);
end
GeneDensity = GeneCounts./(DomainEnd-DomainStart);
figure(2)
bar(GeneDensity, 'k');
xlabel('ID number of imaged TADs');
ylabel('Gene density')
savefig('Gene density.fig')

%% Load in ChrList from each dataset, and remap the cell type based on the combined cell type identification
load('Clustering results final.mat')
Chr_all = [];
for FolderID = 1:length(Folder)
    % find cells in this folder (dataset)
    CellsInFolder = CellList(find([CellList.DatasetID]==FolderID));
    load([Folder{FolderID} '\ChrList.mat']);
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
caxis([0 3])
title('Mean spatial distance');
RedBlue
PlotProp
axis square
title(['N = ' num2str(length(Chr))])
savefig(['Mean spatial distance_LargeScale.fig']);

GenomicDis = zeros(TotalNumTADs,TotalNumTADs);
GenomicDis_pxl = [];
HiC_pxl = [];
Dis_pxl = [];
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs  
        GenomicDis(i,j) = abs(DomainCenter(i)-DomainCenter(j));
        if i~= j  % without diagonal entries
            if isempty(find(DomainsToExclude == i)) && isempty(find(DomainsToExclude == j)) 
                GenomicDis_pxl = [GenomicDis_pxl GenomicDis(i,j)];
                HiC_pxl = [HiC_pxl HiC_chozen(i,j)];
                Dis_pxl = [Dis_pxl Mean(i,j)];
            end
        end
    end
end
figure(4)
loglog(Dis_pxl, 1./HiC_pxl, '.','MarkerSize', 10);
xlabel('Mean spatial distance (um)');
ylabel('Inversed HiC contact frequency');
PlotProp
hold on
axis square
Ind = find(HiC_pxl>0);
x = log(Dis_pxl(Ind))';
y = log(1./HiC_pxl(Ind))';
X = [ones(size(x)), x];
[b, bint] = regress(y,X);
display(['Slope = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])
R = corrcoef(x, y);
title(['R = ' num2str(R(2,1)) '; k = ' num2str(b(2)) '+-' num2str(bint(2,2)-b(2))])
display(['Correlation coefficient = ' num2str(R(2,1))])
x = Dis_pxl';
y_fit = exp(b(1))*x.^b(2);
loglog(x, y_fit, '-r','LineWidth', 4);
hold off

savefig(['HiC correlation coefficient_LargeScale.fig']);

%% Calculate mean spatial distance matrix, normalized matrix, correlation matrix and PCA for each cell type
mkdir figures_TraceAnalysisByClusters_LargeScale
for ii = 1:max(CellTypeID_new)
    Mean = [];
    Std = [];
    SEM = [];
    NofData = [];
    DisListAll = {};
    ChrInGroup = Chr(find([Chr.CellType]==ii));
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            DisList = [];
            for k = 1:length(ChrInGroup)
                if ChrInGroup(k).r(i) == 1 && ChrInGroup(k).r(j) == 1
                    DisList = [DisList ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
                end
            end
            Mean(i,j) = mean(DisList);
            Std(i,j) = std(DisList);
            SEM(i,j) = std(DisList)/(length(DisList))^0.5;
            NofData(i,j) = length(DisList);
            DisListAll{i,j} = DisList;
        end
    end
    NChrInGroup(ii) = length(ChrInGroup);

    figure(5)
    imagesc(Mean)
    colorbar
    caxis([0 3])
    title(['Mean spatial distance for Group ' num2str(ii) ', N=' num2str(length(ChrInGroup))]);
    ylabel('ID number of imaged TADs');
    xlabel('ID number of imaged TADs');
    RedBlue
    PlotProp
    axis square
    hold on % mark pixels that are NaN with dots
    [row, col] = find(isnan(Mean) == 1);
    plot(col,row,'y.','MarkerSize', 20);
    hold off
    saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Mean spatial distance for Group ' num2str(ii) '.jpg'])
    savefig(['figures_TraceAnalysisByClusters_LargeScale/Mean spatial distance for Group ' num2str(ii) '.fig']); 

    % calculate the genomic distances
    GenomicDis = zeros(TotalNumTADs,TotalNumTADs);
    GenomicDis_pxl = [];
    Dis_pxl = [];
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs  
            GenomicDis(i,j) = abs(DomainCenter(i)-DomainCenter(j));
            if i~= j  % without diagonal entries
                if isempty(find(DomainsToExclude == i)) && isempty(find(DomainsToExclude == j)) 
                    GenomicDis_pxl = [GenomicDis_pxl GenomicDis(i,j)];
                    Dis_pxl = [Dis_pxl Mean(i,j)];
                end
            end
        end
    end
    figure(6)
    plot(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
    xlabel('Genomic distance (bp)');
    ylabel('Mean spatial distance (um)');
    PlotProp
    hold on
    [x, IX] = sort(GenomicDis_pxl);
    y = Dis_pxl(IX);
    % fit to get new scaling factor
    Equ = 'b*x^s';
    StartPoint = [mean(y)/(mean(x))^(1/3), 1/3];
    % get rid of NaN values in y if there are any.
    Ind = find(y>0);
    y = y(Ind);
    x = x(Ind);
    [f, gof] = fit(x', y', Equ, 'Start', StartPoint);
    y_fit = f.b*x.^f.s;
    plot(x, y_fit, 'r-', 'LineWidth', 5);
    hold off
    ci = confint(f, 0.95);
    title(['Group ' num2str(ii) ', N=' num2str(length(ChrInGroup)) ', Scaling = ' num2str(f.s) '+-' num2str(ci(2,2)-f.s)])
    saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/scaling for Group ' num2str(ii) '.jpg'])
    savefig(['figures_TraceAnalysisByClusters_LargeScale/scaling for Group ' num2str(ii) '.fig']); 
    
    % calculate normalized distance matrix
    x = GenomicDis_pxl';
    y_fit = f.b*x.^f.s;
    Mean_adjust = Mean; 
    n = 0;
    if NormalizeTowardsMean == 0
        for i = 1:TotalNumTADs
            for j = 1:TotalNumTADs
                if i~=j
                    n = n+1;
                    Mean_adjust(i,j) = Mean(i,j)/y_fit(n);
                elseif i==j && ~isnan(Mean_adjust(i,j))
                    Mean_adjust(i,j) = 1;
                end
            end
        end
    else
        %calculate the running average of spatial distance,
        %use the 30 closest genomic distance entries for each data point.
        for i = 1:TotalNumTADs
            for j = 1:TotalNumTADs
                if i~=j
                    n = n+1;
                    GenomicDis_diff = abs(GenomicDis_pxl - GenomicDis_pxl(n));
                    [GenomicDis_sort, IX] = sort(GenomicDis_diff);
                    Dis_sort = Dis_pxl(IX);
                    Dis_sort = Dis_sort(find(isnan(Dis_sort)~=1));
                    Dis_sort_mean(n) = mean(Dis_sort(1:30));
                    Mean_adjust(i,j) = Mean(i,j)/Dis_sort_mean(n);
                else
                    Mean_adjust(i,j) = 1;
                end
            end
        end
        figure(6)
        hold on
        [GenomicDis_sort, IX] = sort(GenomicDis_pxl);
        plot(GenomicDis_sort, Dis_sort_mean(IX), 'k-', 'LineWidth', 5);
        hold off
        title(['Group ' num2str(ii) ', N=' num2str(length(ChrInGroup)) ', Scaling = ' num2str(f.s) '+-' num2str(ci(2,2)-f.s)])
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/scaling for Group ' num2str(ii) '.jpg'])
        savefig(['figures_TraceAnalysisByClusters_LargeScale/scaling for Group ' num2str(ii) '.fig']);
    end
    figure(7)
    imagesc(Mean_adjust)
    axis square
    RedBlue
    Range = max([max(max(Mean_adjust))-1 1-min(min(Mean_adjust))]);
    caxis([1-Range 1+Range])
    colorbar
    title(['Normalized distance, Group ' num2str(ii) ', N=' num2str(length(ChrInGroup))]);
    ylabel('ID number of imaged TADs');
    xlabel('ID number of imaged TADs');
    PlotProp
    hold on % mark pixels that are NaN with dots
    [row, col] = find(isnan(Mean_adjust) == 1);
    plot(col,row,'y.','MarkerSize', 20);
    hold off
    saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Normalized distance for Group ' num2str(ii) '.jpg'])
    savefig(['figures_TraceAnalysisByClusters_LargeScale/Normalized distance for Group ' num2str(ii) '.fig']);
    
    % calculate correlation matrix
    Pearson = [];
    for i = 1:TotalNumTADs
        for j = 1:TotalNumTADs
            X = [];
            Y = [];
            for k = 1:TotalNumTADs
                if k~=i && k~=j && ~isnan(Mean_adjust(k,i)) && ~isnan(Mean_adjust(k,j))
                    X = [X; Mean_adjust(k,i)];
                    Y = [Y; Mean_adjust(k,j)];
                end
            end
            if ~isempty(X)
                Pearson(i,j) = corr(X,Y, 'type', 'Pearson');
            else
                Pearson(i,j) = NaN;
            end
        end
    end
    
    figure(8)
    imagesc(Pearson)
    axis square

    BlueRed
    colorbar
    title(['Pearson correlation, Group ' num2str(ii) ', N=' num2str(length(ChrInGroup))])
    ylabel('ID number of imaged TADs');
    xlabel('ID number of imaged TADs');
    PlotProp
    hold on % mark pixels that are NaN with dots
    [row, col] = find(isnan(Pearson) == 1);
    plot(col,row,'y.','MarkerSize', 20);
    hold off
    saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Correlation matrix for Group ' num2str(ii) '.jpg'])
    savefig(['figures_TraceAnalysisByClusters_LargeScale/Correlation matrix for Group ' num2str(ii) '.fig']);
        
    % calculate the compartmentalization with PCA
    [coeff,score,latent,tsquared,explained,mu] = pca(Pearson);
    if ~isempty(coeff)
        coeff = coeff(:,1);
        % determine if compartment needs switching based on the correlation
        % with gene desntiy
        if corr(coeff, GeneDensity')<0
            coeff = -coeff;
        end
        CompartmentScores{ii} = coeff;
        figure(9)
        bar(coeff,'b');
        ylabel('Compartment score');
        xlabel('ID number of imaged TADs');
        hold on
        A = find(coeff>0);
        bar(A, coeff(A), 'r');
        hold off

        xlim([0 TotalNumTADs+1])

        PlotProp
        title(['PCA, Group ' num2str(ii) ', N=' num2str(length(ChrInGroup))])
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Compartment scores for Group ' num2str(ii) '.jpg'])
        savefig(['figures_TraceAnalysisByClusters_LargeScale/Compartment scores for Group ' num2str(ii) '.fig']);
        % save the compartment scores for this cell type in an array
        AllCoeff{ii} = coeff;
    end
    
    % calculate the polarization index of cell type and the randomization control
    if ~isempty(coeff) 
        coeffControl = randomizeCompartments(coeff); % generate randomized compartment assignment
        for i = 1:length(ChrInGroup)
            if length(find(ChrInGroup(i).r == 0))<=MaxNumMissedTAD % trace is long enough for polarization analysis
                try
                    ChrInGroup(i).PolarizationIndex = calculatePolarization(ChrInGroup(i),coeff);
                    ChrInGroup(i).PolarizationIndexControl = calculatePolarization(ChrInGroup(i),coeffControl);
                catch
                    ChrInGroup(i).PolarizationIndex = NaN;
                    ChrInGroup(i).PolarizationIndexControl = NaN;
                end
            else
                ChrInGroup(i).PolarizationIndex = NaN;
                ChrInGroup(i).PolarizationIndexControl = NaN;
            end
        end
        PI = [ChrInGroup.PolarizationIndex];
        PI = PI(~isnan(PI));
        PIcontrol = [ChrInGroup.PolarizationIndexControl];
        PIcontrol = PIcontrol(~isnan(PIcontrol));
        figure(10)
        if ~isempty(PI) && ~isempty(PIcontrol)
            dotboxv(PI,  PIcontrol);
            set(gca, 'LineWidth', 2, 'FontSize', 16, 'XTick',[1, 2.5],'XTickLabel', {'Obs', 'Ctr'});
            p = ranksum(PI, PIcontrol); % Wilcoxon rank sum test
            title(['Group ' num2str(ii) ', p = ' num2str(p) ', N=' num2str(length(PI))])
            ylabel('Polarization index');
            axis square
            saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Polarization index for Group ' num2str(ii) '.jpg'])
            savefig(['figures_TraceAnalysisByClusters_LargeScale/Polarization index for Group ' num2str(ii) '.fig']);
        end
        close(10)
    end
    
    % calulate the chromosome surface ratio of all TADs in this cell type (probability of TAD at surface of chromosme)
    AllSurfaceIndices{ii} = [];
    for i = 1:length(ChrInGroup)
        if length(find(ChrInGroup(i).r == 0))<=MaxNumMissedTAD % trace is long enough
            AllSurfaceIndices{ii} = [AllSurfaceIndices{ii}; calculateSurfaceIndex(ChrInGroup(i))];
        end
    end
    if ~isempty(AllSurfaceIndices{ii})
        for i = 1:TotalNumTADs
            List = AllSurfaceIndices{ii}(:,i);
            if length(find(isnan(List)~=1))~=0
                SurfaceRatio(i) = length(find(List == 1))/length(find(isnan(List)~=1));
                SurfaceRatioError(i) = 1.96*sqrt(SurfaceRatio(i)*(1-SurfaceRatio(i))/length(find(isnan(List)~=1)));
            else
                SurfaceRatio(i) = NaN;
                SurfaceRatioError(i) = NaN;
            end
        end
        % save the SurfaceRatio of this cell type in an array.
        AllSurfaceRatio{ii} = SurfaceRatio;
        AllSurfaceRatioError{ii} = SurfaceRatioError;
        % note the AllSurfaceIndices array saves the surface indices for each chromosome for all cell types
    end
    
    % calculate the lamina association ratio for all TADs in this cell
    % type (probability of TAD within certain threshold distance from lamina)
    AllLaminaAssociations = [];
    for i = 1:length(ChrInGroup)
        LaminaAssociations = zeros(1,TotalNumTADs);
        LaminaAssociations(find(ChrInGroup(i).r == 0)) = NaN;
        LaminaAssociations(find(ChrInGroup(i).r==1 & ChrInGroup(i).MinDisToLamina<0.2)) = 1;
        AllLaminaAssociations = [AllLaminaAssociations; LaminaAssociations];
    end
    for i = 1:TotalNumTADs
        List = AllLaminaAssociations(:,i);
        if length(find(isnan(List)~=1))~=0
            LaminaRatio(i) = length(find(List == 1))/length(find(isnan(List)~=1));
            LaminaRatioError(i) = 1.96*sqrt(LaminaRatio(i)*(1-LaminaRatio(i))/length(find(isnan(List)~=1)));
        else
            LaminaRatio(i) = NaN;
            LaminaRatioError(i) = NaN;
        end
    end
    % save the LaminaRatio of this cell type in an array.
    AllLaminaRatio{ii} = LaminaRatio;
    AllLaminaRatioError{ii} = LaminaRatioError;
    
    % calculate the nucleolar association ratio for all TADs in this cell
    % type (probability of TAD within certain threshold distance from nucleolus)
    AllNucleolarAssociations = [];
    for i = 1:length(ChrInGroup)
        NucleolarAssociations = zeros(1,TotalNumTADs);
        NucleolarAssociations(find(ChrInGroup(i).r == 0)) = NaN;
        NucleolarAssociations(find(ChrInGroup(i).r==1 & ChrInGroup(i).MinDisToNucleolar<0.2)) = 1;
        AllNucleolarAssociations = [AllNucleolarAssociations; NucleolarAssociations];
    end
    for i = 1:TotalNumTADs
        List = AllNucleolarAssociations(:,i);
        if length(find(isnan(List)~=1))~=0
            NucleolarRatio(i) = length(find(List == 1))/length(find(isnan(List)~=1));
            NucleolarRatioError(i) = 1.96*sqrt(NucleolarRatio(i)*(1-NucleolarRatio(i))/length(find(isnan(List)~=1)));
        else
            NucleolarRatio(i) = NaN;
            NucleolarRatioError(i) = NaN;
        end
    end
    % save the NucleolarRatio of this cell type in an array.
    AllNucleolarRatio{ii} = NucleolarRatio;
    AllNucleolarRatioError{ii} = NucleolarRatioError;

    % Here calcualte the mean RNA copy numbers within this group
    CellInGroup = CellList(find([CellList.CellType]==ii));
    RNACopyNumberMatrix = [];
    for i = 1:length(CellInGroup)
        RNACopyNumberMatrix = [RNACopyNumberMatrix CellInGroup(i).RNACopyNumber];
    end
    MeanRNACopyNumbers{ii} = mean(RNACopyNumberMatrix,2);
    
    % Save ChrInGrouop with polarization indices
    ChrInGroup_all{ii} = ChrInGroup;
end

%% Analyze if the change of A-B compartment score is largely associated with expression change of Chr19 genes.
load('GenesOfInterestWithOligos_FetalLiver_final');
GeneExpression_HigherScore = [];
GeneExpression_LowerScore = [];
HigherScore = [];
LowerScore = [];
for i = 1:TotalNumTADs
    for ii = 1:6
        if ~isempty(AllCoeff{ii})
            for jj = ii+1:7
                if ~isempty(AllCoeff{jj})
                    if AllCoeff{ii}(i)>AllCoeff{jj}(i)
                        HigherScoreCellTypeID = ii;
                        LowerScoreCellTypeID = jj;
                    else
                        HigherScoreCellTypeID = jj;
                        LowerScoreCellTypeID = ii;
                    end
                    TADid = ChosenDomains(i);
                    for j = 1:length(GenesOfInterestFinal)
                        if TADid == GenesOfInterestFinal(j).TADid % gene is in this TAD
                            GeneExpression_HigherScore = [GeneExpression_HigherScore; MeanRNACopyNumbers{HigherScoreCellTypeID}(j)];
                            GeneExpression_LowerScore = [GeneExpression_LowerScore; MeanRNACopyNumbers{LowerScoreCellTypeID}(j)];
                            HigherScore = [HigherScore; AllCoeff{HigherScoreCellTypeID}(i)];
                            LowerScore = [LowerScore; AllCoeff{LowerScoreCellTypeID}(i)];
                        end
                    end
                end
            end
        end
    end
end
% plot result
figure(1111)
ExpressionFoldChange = GeneExpression_HigherScore./GeneExpression_LowerScore;
ScoreIncrease = HigherScore-LowerScore;
Ind = find(ExpressionFoldChange>1);
plot(ExpressionFoldChange(Ind), ScoreIncrease(Ind), '.','MarkerSize', 20,'MarkerEdgeColor',[213,94,0]/255)
ExpressionFoldIncrease = ExpressionFoldChange(Ind);
ScoreChange = ScoreIncrease(Ind);
hold on
Ind = find(ExpressionFoldChange<=1);
plot(1./ExpressionFoldChange(Ind), -ScoreIncrease(Ind), '.','MarkerSize', 20,'MarkerEdgeColor',[0,158,115]/255)
ExpressionFoldIncrease = [ExpressionFoldIncrease; 1./ExpressionFoldChange(Ind)];
ScoreChange = [ScoreChange; -ScoreIncrease(Ind)];
plot([3 3], [min(ScoreChange)-0.05  max(ScoreChange)+0.05],'--k','LineWidth',2)
hold off
xlim([0.9 max(ExpressionFoldIncrease)+1])
ylim([min(ScoreChange)-0.05  max(ScoreChange)+0.05])
xlabel('Fold increase of gene expression')
ylabel('Compartment score change')
PlotProp
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase vs AB score change.jpg'])
savefig(['figures_TraceAnalysisByClusters_LargeScale/gene expression increase vs AB score change.fig']);
figure(1112)
Ind = find(ExpressionFoldIncrease>2);
pie([length(find(ScoreChange(Ind)>0)) length(find(ScoreChange(Ind)<0))],{'AB score increases', 'AB score decreases'});
title('Gene expression increase > 2-fold')
colormap([[213,94,0]/255; [0,158,115]/255])
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (2-fold or more).jpg'])
savefig(['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (2-fold or more).fig']);
figure(1113)
Ind = find(ExpressionFoldIncrease>3);
pie([length(find(ScoreChange(Ind)>0)) length(find(ScoreChange(Ind)<0))],{'AB score increases', 'AB score decreases'});
title('Gene expression increase > 3-fold')
colormap([[213,94,0]/255; [0,158,115]/255])
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (3-fold or more).jpg'])
savefig(['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (3-fold or more).fig']);
figure(1114)
Ind = find(ExpressionFoldIncrease>4);
pie([length(find(ScoreChange(Ind)>0)) length(find(ScoreChange(Ind)<0))],{'AB score increases', 'AB score decreases'});
title('Gene expression increase > 4-fold')
colormap([[213,94,0]/255; [0,158,115]/255])
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (4-fold or more).jpg'])
savefig(['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (4-fold or more).fig']);

%% plot the lamina association ratio with compartment score
mkdir figures_Lamina % this folder stores the lamina analyses
for ii = 1:max(CellTypeID_new)
    if ~isempty(AllCoeff{ii}) && ~isempty(AllLaminaRatio{ii})
        figure(70)
        subplot(2,1,1)
        bar(AllCoeff{ii},'b');
        ylabel('Compartment score');
        xlabel('ID number of imaged TADs');
        hold on
        A = find(AllCoeff{ii}>0);
        bar(A, AllCoeff{ii}(A), 'r');
        hold off
        xlim([0 TotalNumTADs+1])
        [R, P] = corr(AllCoeff{ii},AllLaminaRatio{ii}', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R = ' num2str(R) ', p = ' num2str(P)]);
        subplot(2,1,2)
        errorbar(1:TotalNumTADs,AllLaminaRatio{ii},AllLaminaRatioError{ii},'o');
        ylabel('Lamina ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        saveas(gcf, ['figures_Lamina/AB score and LaminaRatio for Group ' num2str(ii) '.jpg'])
        
        figure(72)
        errorbar(1:TotalNumTADs,AllLaminaRatio{ii},AllLaminaRatioError{ii},'o','MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 2);
        ylabel('Lamina ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        PlotProp
        savefig(['figures_Lamina/LaminaRatio for Group ' num2str(ii) '.fig'])
        
        figure(71)
        plot(AllCoeff{ii},AllLaminaRatio{ii}','.','MarkerSize', 20);
        hold on
        y = AllLaminaRatio{ii}';
        x = AllCoeff{ii};
        X = [ones(size(x)), x];
        b = regress(y,X);
        x_fit = [min(AllCoeff{ii})-0.05; max(AllCoeff{ii})+0.05];
        X_fit = [ones(size(x_fit)), x_fit];
        y_fit = X_fit*b;
        plot(x_fit, y_fit,'-','LineWidth',2)
        hold off
        xlabel('Compartment score');
        ylabel('Lamina ratio')
        [R, P] = corr(AllCoeff{ii},AllLaminaRatio{ii}', 'type', 'Pearson');
        R_LAD(ii) = R;
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        xlim([min(AllCoeff{ii})-0.05 max(AllCoeff{ii})+0.05])
        ylim([min(AllLaminaRatio{ii})-0.05 max(AllLaminaRatio{ii})+0.05])
        PlotProp
        savefig(['figures_Lamina/AB score and LaminaRatio correaltion for Group ' num2str(ii) '.fig']);
    
        if ii==3 || ii==6 || ii==7
            if ii == 3
                Color = [230,159,0]/255;
            elseif ii == 6
                Color = [204,121,167]/255;
            else
                Color = [0,114,178]/255;
            end
            figure(73)
            hold on
            plot(AllCoeff{ii},AllLaminaRatio{ii}','.','MarkerSize', 20,'MarkerEdgeColor',Color);
            y = AllLaminaRatio{ii}';
            x = AllCoeff{ii};
            X = [ones(size(x)), x];
            b = regress(y,X);
            x_fit = [min(AllCoeff{ii})-0.05; max(AllCoeff{ii})+0.05];
            X_fit = [ones(size(x_fit)), x_fit];
            y_fit = X_fit*b;
            plot(x_fit, y_fit,'-','LineWidth',2,'Color',Color);
            hold off
            xlabel('Compartment score');
            ylabel('Lamina ratio')
            PlotProp
            if ii == 7
                legend('Hep','Hep','Ery','Ery','Pro','Pro')
                savefig(['figures_Lamina/AB score and LaminaRatio correaltion (JOINT).fig']);
            end
        end
    end 
end
% plot the correlations for all cell types
figure(74)
barh(1:7,R_LAD(1:7),'LineWidth',2,'FaceColor',[86,180,233]/255)
xlabel('Correlation between compartment score and lamina ratio')
ylabel('Cell type')
yticklabels({'Other','Epithelial cell','Hepatocyte','Macrophage','Megakaryocyte','Erythroblast','Proerythroblast'})
PlotProp
savefig(['figures_Lamina/All lamina correlations.fig']);

%% plot the nucleolar ratio with compartment score 
mkdir figures_Nucleolus % this folder stores the nucleolus analyses
for ii = 1:max(CellTypeID_new)
    if ~isempty(AllCoeff{ii}) && ~isempty(AllNucleolarRatio{ii})
        figure(80)
        subplot(2,1,1)
        bar(AllCoeff{ii},'b');
        ylabel('Compartment score');
        xlabel('ID number of imaged TADs');
        hold on
        A = find(AllCoeff{ii}>0);
        bar(A, AllCoeff{ii}(A), 'r');
        hold off
        xlim([0 TotalNumTADs+1])
        [R, P] = corr(AllCoeff{ii},AllNucleolarRatio{ii}', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R = ' num2str(R) ', p = ' num2str(P)]);
        subplot(2,1,2)
        errorbar(1:TotalNumTADs,AllNucleolarRatio{ii},AllNucleolarRatioError{ii},'o');
        ylabel('Nucelolar ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        saveas(gcf, ['figures_Nucleolus/AB score and NucelolarRatio for Group ' num2str(ii) '.jpg'])
        
        figure(82)
        errorbar(1:TotalNumTADs,AllNucleolarRatio{ii},AllNucleolarRatioError{ii},'o','MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 2);
        ylabel('Nucelolar ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        PlotProp
        savefig(['figures_Nucleolus/NucelolarRatio for Group ' num2str(ii) '.fig'])
        
        % plot the correlation of AB score with nucleolar ratio for just the
        % A scores or the B scores
        figure(81)
        plot(AllCoeff{ii},AllNucleolarRatio{ii}','.','MarkerSize', 20);
        hold on
        y = AllNucleolarRatio{ii}';
        x = AllCoeff{ii};
        X = [ones(size(x)), x];
        b = regress(y,X);
        x_fit = [min(AllCoeff{ii})-0.05; max(AllCoeff{ii})+0.05];
        X_fit = [ones(size(x_fit)), x_fit];
        y_fit = X_fit*b;
        plot(x_fit, y_fit,'-','LineWidth',2)
        hold off
        xlabel('Compartment score');
        ylabel('Nucleolar ratio')
        [R, P] = corr(AllCoeff{ii},AllNucleolarRatio{ii}', 'type', 'Pearson');
        R_NAD(ii) = R;
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        xlim([min(AllCoeff{ii})-0.05 max(AllCoeff{ii})+0.05])
        ylim([0 max(AllNucleolarRatio{ii})+0.01])
        PlotProp
        savefig(['figures_Nucleolus/AB score and NucelolarRatio correaltion for Group ' num2str(ii) '.fig']);
    
        if ii==3 || ii==6 || ii==7
            if ii == 3
                Color = [230,159,0]/255;
            elseif ii == 6
                Color = [204,121,167]/255;
            else
                Color = [0,114,178]/255;
            end
            figure(83)
            hold on
            plot(AllCoeff{ii},AllNucleolarRatio{ii}','.','MarkerSize', 20,'MarkerEdgeColor',Color);
            y = AllNucleolarRatio{ii}';
            x = AllCoeff{ii};
            X = [ones(size(x)), x];
            b = regress(y,X);
            x_fit = [min(AllCoeff{ii})-0.05; max(AllCoeff{ii})+0.05];
            X_fit = [ones(size(x_fit)), x_fit];
            y_fit = X_fit*b;
            plot(x_fit, y_fit,'-','LineWidth',2,'Color',Color);
            hold off
            xlabel('Compartment score');
            ylabel('Nucleolar ratio')
            PlotProp
            if ii == 7
                legend('Hep','Hep','Ery','Ery','Pro','Pro')
                savefig(['figures_Nucleolus/AB score and NucelolarRatio correaltion (JOINT).fig']);
            end
        end
    end
end
% plot the correlations for all cell types
figure(84)
barh(1:7,R_NAD(1:7),'LineWidth',2,'FaceColor',[0,158,115]/255)
xlabel('Correlation between compartment score and nucleolar ratio')
ylabel('Cell type')
yticklabels({'Other','Epithelial cell','Hepatocyte','Macrophage','Megakaryocyte','Erythroblast','Proerythroblast'})
PlotProp
savefig(['figures_Nucleolus/All nucleolar correlations.fig']);

%% plot the Surface ratio, compartment score, and constitutive LADs along the chromosome for each cell type
% load in the constitutive lad coordiantes from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17051
% these are mm9 coordinates.
mkdir figures_SurfaceRatio % this folder stores the surface ratio analyses
fid = fopen('GSE17051_cLAD_regions.bed', 'r');
tline = fgetl(fid);
i = 0;
while ischar(tline)
    C = strsplit(tline);
    if strcmpi(C{1}, ChrName)
        i = i+1;
        LADStart(i) = str2num(C{2});
        LADEnd(i) = str2num(C{3});
    end
    tline = fgetl(fid);
end
fclose(fid);
% find traced TAD that overlap with LAD
LADid = [];
for i = 1:TotalNumTADs
    for j = 1:length(LADStart)
        if DomainStart(i) >= LADEnd(j) || DomainEnd(i) <= LADStart(j)
        else
            display(['TAD ' num2str(i) ' is overlaping with LAD.']);
            LADid = [LADid i];
            break
        end
    end
end

for ii = 1:max(CellTypeID_new)
    if ~isempty(AllCoeff{ii}) && ~isempty(AllSurfaceRatio{ii})
        figure(90)
        subplot(2,1,1)
        bar(AllCoeff{ii},'b');
        ylabel('Compartment score');
        xlabel('ID number of imaged TADs');
        hold on
        A = find(AllCoeff{ii}>0);
        bar(A, AllCoeff{ii}(A), 'r');
        bar(LADid, AllCoeff{ii}(LADid), 'k'); % LADs are plotted as black
        hold off
        xlim([0 TotalNumTADs+1])
        [R, P] = corr(AllCoeff{ii},AllSurfaceRatio{ii}', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R = ' num2str(R) ', p = ' num2str(P)]);
        subplot(2,1,2)
        errorbar(1:TotalNumTADs,AllSurfaceRatio{ii},AllSurfaceRatioError{ii},'o');
        ylabel('Surface ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        saveas(gcf, ['figures_SurfaceRatio/AB score and SR for Group ' num2str(ii) '.jpg'])
        
        % Plot the SR by itself
        figure(93)
        errorbar(1:TotalNumTADs,AllSurfaceRatio{ii},AllSurfaceRatioError{ii},'o','MarkerEdgeColor','red','MarkerFaceColor','red', 'LineWidth', 2);
        title(['Group' num2str(ii)]);
        ylabel('Surface ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        PlotProp
        savefig(['figures_SurfaceRatio/SR for Group ' num2str(ii) '.fig'])
        
        % plot the correlation of AB score with surface ratio for just the
        % A scores or the B scores
        figure(91)
        plot(AllCoeff{ii},AllSurfaceRatio{ii}','.b','MarkerSize', 20);
        A = find(AllCoeff{ii}>0);
        hold on
        plot(AllCoeff{ii}(A),AllSurfaceRatio{ii}(A)','.r','MarkerSize', 20);
        plot([0 0], [min(AllSurfaceRatio{ii})-0.05 max(AllSurfaceRatio{ii})+0.05],'--k','LineWidth',2)
        % generate line fits for the A and B regions
        y = AllSurfaceRatio{ii}(A)';
        x = AllCoeff{ii}(A);
        X = [ones(size(x)), x];
        b = regress(y,X);
        x_fit = [0; max(AllCoeff{ii})+0.05];
        X_fit = [ones(size(x_fit)), x_fit];
        y_fit = X_fit*b;
        plot(x_fit, y_fit,'-r','LineWidth',2)
        B = find(AllCoeff{ii}<0);
        y = AllSurfaceRatio{ii}(B)';
        x = AllCoeff{ii}(B);
        X = [ones(size(x)), x];
        b = regress(y,X);
        x_fit = [min(AllCoeff{ii})-0.05; 0];
        X_fit = [ones(size(x_fit)), x_fit];
        y_fit = X_fit*b;
        plot(x_fit, y_fit,'-b','LineWidth',2)
        hold off
        xlabel('Compartment score');
        ylabel('Surface ratio')
        [RA, PA] = corr(AllCoeff{ii}(A),AllSurfaceRatio{ii}(A)', 'type', 'Pearson');
        [RB, PB] = corr(AllCoeff{ii}(B),AllSurfaceRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii)]);
        PlotProp
        xlim([min(AllCoeff{ii})-0.05 max(AllCoeff{ii})+0.05])
        ylim([min(AllSurfaceRatio{ii})-0.05 max(AllSurfaceRatio{ii})+0.05])
        text(min(AllCoeff{ii}), min(AllSurfaceRatio{ii}),['RB = ' num2str(RB) ', pB = ' num2str(PB) ', RA = ' num2str(RA) ', pA = ' num2str(PA)])
        savefig(['figures_SurfaceRatio/AB score and SR correaltion for Group ' num2str(ii) '.fig']);
    end
end
%% Plot the PI vs the Lamina/Nucleolar contacts for different cell types
mkdir figures_PIandLaminaOrNucelolarContacts
for ii = 1:7
    ChrInGroup = ChrInGroup_all{ii};
    B = find(AllCoeff{ii}<0);
    
    PI = [ChrInGroup.PolarizationIndex];
    AllLaminaAssociations = [];
    LaminaAssociationCout = [];
    for i = 1:length(ChrInGroup)
        LaminaAssociations = zeros(1,TotalNumTADs);
        LaminaAssociations(find(ChrInGroup(i).r == 0)) = NaN;
        LaminaAssociations(find(ChrInGroup(i).r==1 & ChrInGroup(i).MinDisToLamina<0.2)) = 1;
        AllLaminaAssociations = [AllLaminaAssociations; LaminaAssociations];
        LaminaAssociations = LaminaAssociations(B);
        LaminaAssociationCout =[LaminaAssociationCout length(find(LaminaAssociations==1))];
    end
    AllNucleolarAssociations = [];
    NucleolarAssociationCout = [];
    for i = 1:length(ChrInGroup)
        NucleolarAssociations = zeros(1,TotalNumTADs);
        NucleolarAssociations(find(ChrInGroup(i).r == 0)) = NaN;
        NucleolarAssociations(find(ChrInGroup(i).r==1 & ChrInGroup(i).MinDisToNucleolar<0.2)) = 1;
        AllNucleolarAssociations = [AllNucleolarAssociations; NucleolarAssociations];
        NucleolarAssociations = NucleolarAssociations(B);
        NucleolarAssociationCout =[NucleolarAssociationCout length(find(NucleolarAssociations==1))];
    end
    Ind = find(~isnan(PI));
    PI = PI(Ind);
    LaminaAssociationCout = LaminaAssociationCout(Ind);
    NucleolarAssociationCout = NucleolarAssociationCout(Ind);
    figure(900)
    grouping = LaminaAssociationCout;
    grouping(find(grouping>10))=10;
    boxplot(PI,grouping,'Whisker',0,'Symbol','')
    [R, P] = corr(PI',LaminaAssociationCout', 'type', 'Pearson');
    title(['R = ' num2str(R) ', p = ' num2str(P)]);
    PlotProp
    ylabel('Polarization index')
    xlabel('Comparment B-lamina association count')
    savefig(['figures_PIandLaminaOrNucelolarContacts/PI and lamina contact for group' num2str(ii) '.fig'])
    
    figure(901)
    grouping = NucleolarAssociationCout;
    grouping(find(grouping>5))=5;
    boxplot(PI,grouping,'Whisker',0,'Symbol','')
    [R, P] = corr(PI',NucleolarAssociationCout', 'type', 'Pearson');
    title(['R = ' num2str(R) ', p = ' num2str(P)]);
    PlotProp
    ylabel('Polarization index')
    xlabel('Comparment B-nucleoli association count')
    savefig(['figures_PIandLaminaOrNucelolarContacts/PI and nucleoli contact for group' num2str(ii) '.fig'])
end


%% save the workspace
save('run6workspace.mat');
save('AllCoeff.mat','AllCoeff');
save('Chr.mat','Chr');

%% Plot the correlations of lamina/nucleolar association ratios with compartment scores on the same plot
figure(10000)
h = barh(1:7,[R_LAD(1:7)' R_NAD(1:7)'],'grouped','LineWidth',1)
h(1).FaceColor = [0,158,115]/255;
h(2).FaceColor = [213,94,0]/255;

yticklabels({'Other','Epithelial cell','Hepatocyte','Macrophage','Megakaryocyte','Erythroblast','Proerythroblast'})
PlotProp
savefig(['All lamina and nucleolar correlations.fig']);

%% Plot the PI vs the Lamina/Nucleolar contacts for different cell types
% consider 4 goups: with either contacts, with LA but no NA, with NA but no
% LA, and with both LA and NA.
mkdir figures_PIandLaminaOrNucelolarContacts
for ii = 1:7
    ChrInGroup = ChrInGroup_all{ii};
    B = find(AllCoeff{ii}<0);
    
    PI = [ChrInGroup.PolarizationIndex];
    AllLaminaAssociations = [];
    LaminaAssociationCout = [];
    for i = 1:length(ChrInGroup)
        LaminaAssociations = zeros(1,TotalNumTADs);
        LaminaAssociations(find(ChrInGroup(i).r == 0)) = NaN;
        LaminaAssociations(find(ChrInGroup(i).r==1 & ChrInGroup(i).MinDisToLamina<0.2)) = 1;
        AllLaminaAssociations = [AllLaminaAssociations; LaminaAssociations];
        LaminaAssociations = LaminaAssociations(B);
        LaminaAssociationCout =[LaminaAssociationCout length(find(LaminaAssociations==1))];
    end
    AllNucleolarAssociations = [];
    NucleolarAssociationCout = [];
    for i = 1:length(ChrInGroup)
        NucleolarAssociations = zeros(1,TotalNumTADs);
        NucleolarAssociations(find(ChrInGroup(i).r == 0)) = NaN;
        NucleolarAssociations(find(ChrInGroup(i).r==1 & ChrInGroup(i).MinDisToNucleolar<0.2)) = 1;
        AllNucleolarAssociations = [AllNucleolarAssociations; NucleolarAssociations];
        NucleolarAssociations = NucleolarAssociations(B);
        NucleolarAssociationCout =[NucleolarAssociationCout length(find(NucleolarAssociations==1))];
    end
    Ind = find(~isnan(PI));
    PI = PI(Ind);
    LaminaAssociationCout = LaminaAssociationCout(Ind);
    NucleolarAssociationCout = NucleolarAssociationCout(Ind);
    
    figure(1)
    dotboxh4(PI(find(LaminaAssociationCout>0 & NucleolarAssociationCout>0)),...
        PI(find(LaminaAssociationCout>0 & NucleolarAssociationCout==0)),...
        PI(find(LaminaAssociationCout==0 & NucleolarAssociationCout>0)),...
        PI(find(LaminaAssociationCout==0 & NucleolarAssociationCout==0)));
    set(gca, 'LineWidth', 1, 'FontSize', 16, 'YTick',[],'YTickLabel', {});
%     title(['Group ' num2str(ii) ', p = ' num2str(p)])
%     xlabel('Polarization index');
    savefig(['figures_PIandLaminaOrNucelolarContacts/PI and nucleolar and or lamina contact for group' num2str(ii) '.fig'])
    saveas(gcf,['figures_PIandLaminaOrNucelolarContacts/PI and nucleolar and or lamina contact for group' num2str(ii) '.png'])
    close(1)
end

