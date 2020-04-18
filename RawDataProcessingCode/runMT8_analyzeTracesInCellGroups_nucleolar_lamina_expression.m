% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

% updated on 190711 to analyze the correltion between gene expression and
% comparmentalization; also save chromosome list to facilitate joint 
% analysis of multiple datasets.

% updated on 190704 to extract the lamina profile and caculate the lamina
% association ratio

% updated on 190612 to include nucleolar abstraction and analysis of
% association with nucleolus

% updated on 190529 to include the consitutive LAD coordinates from Bas van
% Steensel's work and ask if they are in B compartments; also if the
% constitutive LADs have different surface ratio in comparison to all other
% TADs or other TADs in B.

% updated on 190526 to calculate the change in surface ratio for TAD that
% switched A-B compartment in different cell types, and also the change in
% gene expression in those TADs.

% updated on 190524 to search for genes correlated with AB scores and
% chromosome surface indices; also generate normalized matrix with the
% running average method to calcualte expected distances between genomic
% loci

% updated on 190521 to perform compartment polarization analysis
% updated on 190508 to add compartment anlaysis

clear all
close all
NuclearNumImage = 181; % Number of images in nuclear stack and its bead stack 
NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
StepSize = 0.2; %um
TotalNumTADs = 50;
ImageSize = 1536; % number of pxls
NFOV = 14; % number of fields of views
UmPerPxl = 0.108;
DomainStartsAndEnds = 'DomainStartsAndEnds19.txt';
load('ChosenDomains19.mat');
DomainsToExclude = []; % list the domains to exclude, e.g. DomainsToExclude = [15, 2];
ChrName = 'chr19';
MaxNumMissedTAD = 8; % maximum number of missed TADs in polarization analysis
NormalizeTowardsMean = 1; % 1 means using the running average method to calcualte the expected distance between genomic loci
AdaptiveThresholdingSensitivity = 0.5;
NuclearImagePath = 'sequential/Laser405_0_'; % set the path of the nuclear images

%% make new folder
mkdir figures_TraceAnalysisByClusters_LargeScale
mkdir figures_Nucleolus % this folder stores the nucleolus abstraction 
mkdir figures_Lamina % this folder stores the lamina abstraction
%% load the genomic coordiantes of domains
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
figure(1)
bar(GeneDensity, 'k');
xlabel('ID number of imaged TADs');
ylabel('Gene density')
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Gene density.jpg'])

%% load in chromosome traces
load('FinalClusteringResults.mat');
n = 0;
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
    
    % here extract nucleolar coordinates
    [ImageStack, InfoFile] = ReadZStack_MultiChannel(['sequential/STORM2_00_' FOVid],NumImage,FramesToWait,TotalNumChannels,1); % updated on 190717
    ImageMax = medfilt2(max(ImageStack, [],3)); 
    figure(1000)
    subplot(1,2,1)
    imagesc(ImageMax)
    colormap gray
    axis equal
    NormalizationFactor = max(max(ImageMax));
    T = adaptthresh(ImageMax/NormalizationFactor,AdaptiveThresholdingSensitivity);
    Nucleolus = zeros(size(ImageStack));
    load('DeltaZ.mat');
    load('tform.mat');
    for i = 1:size(ImageStack,3)
        CurrentImage = medfilt2(ImageStack(:,:,i));
        BW = imbinarize(CurrentImage/NormalizationFactor,T);
        % Here warp the Cy5 image into the Cy3 channel.
        BW = imtransform(BW, tform, 'XData', [1 ImageSize], 'Ydata', [1 ImageSize]);
        BW(find(BW<0.5))=0;
        BW(find(BW>=0.5))=1;
        Nucleolus(:,:,i) = BW;
    end
    subplot(1,2,2)
    imagesc(max(Nucleolus, [],3))
    colormap gray
    axis equal
    saveas(gcf, ['figures_Nucleolus/Nucelolus' FOVid '.jpg']);
    Ind = find(Nucleolus);
    [row,col,z] = ind2sub(size(Nucleolus),Ind); 
    z = z-DeltaZ; % cancel the color shift in Z
    NucelolarZ = z*StepSize;
    NucelolarX = col*UmPerPxl;
    NucelolarY = row*UmPerPxl;
    figure(1001)
    scatter3(NucelolarX,NucelolarY,NucelolarZ,'.');
    axis equal
    xlabel('x')
    ylabel('y')
    zlabel('z')
    axis ij
    savefig(['figures_Nucleolus/Nucelolus3D' FOVid '.fig']);

    % load drift parameters between the lamina images and the hyb0 images
    if exist(['LaminaDriftParams/DriftParams' FOVid '.mat'])==2
        load(['LaminaDriftParams/DriftParams' FOVid '.mat']);
        
        % here extract lamina coordinates
        [ImageStack, InfoFile] = ReadZStack([NuclearImagePath FOVid],NuclearNumImage,FramesToWait);
        I = mean(ImageStack,3);
        I = I/max(max(I));
        figure(2000)
        subplot(2,2,1)
        imagesc(I) % plot mean
        colormap gray
        axis equal
        T = adaptthresh(I);
        I = I./T;
        Min = quantile(I(:), 0.25);
        Max = quantile(I(:), 0.75);
        I = (I-Min)/(Max-Min);
        I(find(I<0)) = 0;
        I(find(I>1)) = 1;
        subplot(2,2,2)
        imagesc(I) % plot mean with adaptive background removal
        colormap gray
        axis equal
        se = strel('disk', 25);
        Ie = imerode(I,se);
        I = imreconstruct(Ie,I);
        subplot(2,2,3)
        imagesc(I) % plot the openned image
        colormap gray
        axis equal
        hy = fspecial('sobel');
        hx = hy';
        EdgeStack = [];
        for i  = 1:size(ImageStack,3)
            I = ImageStack(:,:,i);
            I = I/max(max(I));
            T = adaptthresh(I);
            I = I./T;
            Min = quantile(I(:), 0.25);
            Max = quantile(I(:), 0.75);
            I = (I-Min)/(Max-Min);
            I(find(I<0)) = 0;
            I(find(I>1)) = 1;
            Ie = imerode(I,se);
            I = imreconstruct(Ie,I);
            Iy = imfilter(double(I), hy, 'replicate');% calculate gradient image
            Ix = imfilter(double(I), hx, 'replicate');% calculate gradient image
            gradmag = sqrt(Ix.^2 + Iy.^2);% calculate gradient image
            gradmag_filtered = imgaussfilt(gradmag,5); % guassian filter
            bwedge = imbinarize(gradmag_filtered,'adaptive','Sensitivity', 0.1);% binarize 
            EdgeStack = cat(3, EdgeStack, bwedge);
        end
        subplot(2,2,4)
        imagesc(mean(EdgeStack,3))
        colormap gray
        axis equal
        saveas(gcf, ['figures_Lamina/Lamina' FOVid '.jpg']);
        Ind = find(EdgeStack);
        [row,col,z] = ind2sub(size(EdgeStack),Ind);
        % apply drift correction and convert to real distance
        LaminaZ = (z- Zdrift)*StepSize;
        LaminaX = (col- Xdrift)*UmPerPxl;
        LaminaY = (row- Ydrift)*UmPerPxl;
        figure(2001)
        scatter3(LaminaX,LaminaY,LaminaZ,'.');
        axis equal
        xlabel('x')
        ylabel('y')
        zlabel('z')
        axis ij
        savefig(['figures_Lamina/Lamina3D' FOVid '.fig']);
    end
    
    
    if exist(['Traces_LargeScale\TraceArrayRefined' FOVid '.mat'])==2 && exist(['tformsWGA/tformWGA_' FOVid '.mat'])==2
        load(['Traces_LargeScale\TraceArrayRefined' FOVid '.mat']);
        load(['tformsWGA/tformWGA_' FOVid '.mat']);
        
        %find the cell group types in this FOV.
        I = zeros(ImageSize);
        I2 = zeros(ImageSize);
        Ind = find([CellList.FOV] == jj);
        CellsInFOV = CellList(Ind);
        for i = 1:length(CellsInFOV)
            I(CellsInFOV(i).PixelList) = CellsInFOV(i).CellType;
            I2(CellsInFOV(i).PixelList) = i;
        end
        
        for i = 1:length(TraceArray)
            n = n+1;
            Chr(n).x = zeros(TotalNumTADs,1);
            Chr(n).y = zeros(TotalNumTADs,1);
            Chr(n).z = zeros(TotalNumTADs,1);
            Chr(n).r = zeros(TotalNumTADs,1);
            Chr(n).x(TraceArray{i}(:,end)) = TraceArray{i}(:,1);
            Chr(n).y(TraceArray{i}(:,end)) = TraceArray{i}(:,2);
            Chr(n).z(TraceArray{i}(:,end)) = TraceArray{i}(:,3);
            Chr(n).r(TraceArray{i}(:,end)) = 1;
            % calculate the minimum distance to nucleolar pixels
            Chr(n).MinDisToNucleolar = (pdist2([NucelolarX, NucelolarY,NucelolarZ],[Chr(n).x, Chr(n).y, Chr(n).z],'euclidean','Smallest',1))';
            
            % calculate the minimum distance to lamina pixels
            if exist(['LaminaDriftParams/DriftParams' FOVid '.mat'])==2
                Chr(n).MinDisToLamina = (pdist2([LaminaX, LaminaY,LaminaZ],[Chr(n).x, Chr(n).y, Chr(n).z],'euclidean','Smallest',1))';
            else
                Chr(n).MinDisToLamina = NaN;
            end
            
            x = TraceArray{i}(:,1)/UmPerPxl;
            y = TraceArray{i}(:,2)/UmPerPxl;
            % take into account the WGA image drift
            [X,Y] = transformPointsInverse(tform,x,y);
            
            %find the corresponding cell group type of this trace.
            MeanY = round(mean(Y));
            MeanX = round(mean(X));
            if MeanY>0 && MeanY<=ImageSize && MeanX>0 && MeanX<=ImageSize 
                Chr(n).CellType = I(MeanY,MeanX);
                if I2(MeanY,MeanX)>0
                    Chr(n).RNACopyNumber = CellsInFOV(I2(MeanY,MeanX)).RNACopyNumber;
                else
                    Chr(n).RNACopyNumber = [];   
                end
            else
                Chr(n).CellType = 0;
                Chr(n).RNACopyNumber = [];                
            end
            
            % this section was added on 190725 to facilitate the combined
            % analysis of mutliple folders (datasets)
            Chr(n).MeanYWGA = MeanY;
            Chr(n).MeanXWGA = MeanX;
            Chr(n).FOV = jj;            
        end
    end
end
display([num2str(length(find([Chr.CellType]~=0))) ' copies of Chr were analyzed.'])
save('ChrList.mat','Chr') %added on 190711 to facilitate joint analysis of multiple datasets.

%% Calculate mean spatial distance matrix, normalized matrix, correlation matrix and PCA for each cell type
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

    figure(2)
    imagesc(Mean)
    colorbar
    caxis([0 3])
    title(['Mean spatial distance for Group ' num2str(ii) ', N=' num2str(length(ChrInGroup))]);
    ylabel('ID number of imaged TADs');
    xlabel('ID number of imaged TADs');
    ColorMap = load('RedBlue.txt');
    colormap(ColorMap/255);
    PlotProp
    axis square
    hold on % mark pixels that are NaN with dots
    [row, col] = find(isnan(Mean) == 1);
    plot(col,row,'y.','MarkerSize', 20);
    hold off
    saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Mean spatial distance for Group ' num2str(ii) '.jpg'])
%     savefig(['figures_TraceAnalysisByClusters_LargeScale/Mean spatial distance for Group ' num2str(ii) '.fig']); 

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
    figure(3)
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
        figure(3)
        hold on
        [GenomicDis_sort, IX] = sort(GenomicDis_pxl);
        plot(GenomicDis_sort, Dis_sort_mean(IX), 'k-', 'LineWidth', 5);
        hold off
        title(['Group ' num2str(ii) ', N=' num2str(length(ChrInGroup)) ', Scaling = ' num2str(f.s) '+-' num2str(ci(2,2)-f.s)])
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/scaling for Group ' num2str(ii) '.jpg'])
    end
    figure(4)
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
    
    figure(5)
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
        figure(6)
        bar(coeff,'b');
        ylabel('Compartment score');
        xlabel('ID number of imaged TADs');
        hold on
        A = find(coeff>0);
        bar(A, coeff(A), 'r');
        hold off

        xlim([0 TotalNumTADs+1])
        axis square
        PlotProp
        title(['PCA, Group ' num2str(ii) ', N=' num2str(length(ChrInGroup))])
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Compartment scores for Group ' num2str(ii) '.jpg'])
        % save the compartment scores for this cell type in an array
        AllCoeff{ii} = coeff;
    end
    
    % calculate the polarization index of cell type and the randomization control
    if ~isempty(coeff) 
        coeffControl = randomizeCompartments(coeff); % generate randomized compartment assignment
        for i = 1:length(ChrInGroup)
            if length(find(ChrInGroup(i).r == 0))<=MaxNumMissedTAD % trace is long enough for polarization analysis
                ChrInGroup(i).PolarizationIndex = calculatePolarization(ChrInGroup(i),coeff);
                ChrInGroup(i).PolarizationIndexControl = calculatePolarization(ChrInGroup(i),coeffControl);
            else
                ChrInGroup(i).PolarizationIndex = NaN;
                ChrInGroup(i).PolarizationIndexControl = NaN;
            end
        end
        PI = [ChrInGroup.PolarizationIndex];
        PI = PI(~isnan(PI));
        PIcontrol = [ChrInGroup.PolarizationIndexControl];
        PIcontrol = PIcontrol(~isnan(PIcontrol));
        figure(7)
        if ~isempty(PI) && ~isempty(PIcontrol)
            dotboxv(PI,  PIcontrol);
            set(gca, 'LineWidth', 2, 'FontSize', 16, 'XTick',[1, 2.5],'XTickLabel', {'Obs', 'Ctr'});
            p = ranksum(PI, PIcontrol); % Wilcoxon rank sum test
            title(['Group ' num2str(ii) ', p = ' num2str(p) ', N=' num2str(length(PI))])
            ylabel('Polarization index');
            axis square
            saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/Polarization index for Group ' num2str(ii) '.jpg'])
        end
        close(7)
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
    % save the NucleolarRatio of this cell type in an array.
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
plot(ExpressionFoldChange(Ind), ScoreIncrease(Ind), '.')
ExpressionFoldIncrease = ExpressionFoldChange(Ind);
ScoreChange = ScoreIncrease(Ind);
hold on
Ind = find(ExpressionFoldChange<=1);
plot(1./ExpressionFoldChange(Ind), -ScoreIncrease(Ind), '.')
ExpressionFoldIncrease = [ExpressionFoldIncrease; 1./ExpressionFoldChange(Ind)];
ScoreChange = [ScoreChange; -ScoreIncrease(Ind)];
hold off
xlabel('Fold increase of gene expression')
ylabel('Compartment score change')
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase vs AB score change.jpg'])

figure(1112)
Ind = find(ExpressionFoldIncrease>2);
pie([length(find(ScoreChange(Ind)>0)) length(find(ScoreChange(Ind)<0))],{'AB score increases', 'AB score decreases'});
title('Gene expression increase > 2-fold')
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (2-fold or more).jpg'])
figure(1113)
Ind = find(ExpressionFoldIncrease>3);
pie([length(find(ScoreChange(Ind)>0)) length(find(ScoreChange(Ind)<0))],{'AB score increases', 'AB score decreases'});
title('Gene expression increase > 3-fold')
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (3-fold or more).jpg'])
figure(1114)
Ind = find(ExpressionFoldIncrease>4);
pie([length(find(ScoreChange(Ind)>0)) length(find(ScoreChange(Ind)<0))],{'AB score increases', 'AB score decreases'});
title('Gene expression increase > 4-fold')
saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/gene expression increase pie chart (4-fold or more).jpg'])



%% plot the lamina association ratio with compartment score or surface ratio
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
        
        % plot the correlation of AB score with lamina ratio for just the
        % A scores or the B scores
        figure(71)
        subplot(2,2,1)
        A = find(AllCoeff{ii}>0);
        plot(AllCoeff{ii}(A),AllLaminaRatio{ii}(A)','*');
        xlabel('Compartment score of A TAD');
        ylabel('Lamina ratio')
        [R, P] = corr(AllCoeff{ii}(A),AllLaminaRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(2,2,2)
        B = find(AllCoeff{ii}<0);
        plot(AllCoeff{ii}(B),AllLaminaRatio{ii}(B)','*');
        xlabel('Compartment score of B TAD');
        ylabel('Lamina ratio')
        [R, P] = corr(AllCoeff{ii}(B),AllLaminaRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        % plot the correlation between lamina ratio and gene density for
        % just the A scores or the B scores
        subplot(2,2,3)
        A = find(AllCoeff{ii}>0);
        plot(GeneDensity(A),AllLaminaRatio{ii}(A)','*');
        xlabel('Gene density in A TAD');
        ylabel('Lamina ratio');
        [R, P] = corr(GeneDensity(A)',AllLaminaRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(2,2,4)
        B = find(AllCoeff{ii}<0);
        plot(GeneDensity(B),AllLaminaRatio{ii}(B)','*');
        xlabel('Gene density in B TAD');
        ylabel('Lamina ratio');
        [R, P] = corr(GeneDensity(B)',AllLaminaRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        saveas(gcf, ['figures_Lamina/AB score LR and gene density for Group ' num2str(ii) '.jpg'])
    end

    
    if ~isempty(AllSurfaceRatio{ii}) && ~isempty(AllLaminaRatio{ii}) && ~isempty(AllCoeff{ii})
        figure(72)
        subplot(2,1,1)
        errorbar(1:TotalNumTADs,AllSurfaceRatio{ii},AllSurfaceRatioError{ii},'o');
        ylabel('Surface ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        [R, P] = corr(AllSurfaceRatio{ii}',AllLaminaRatio{ii}', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R = ' num2str(R) ', p = ' num2str(P)]);
        
        subplot(2,1,2)
        errorbar(1:TotalNumTADs,AllLaminaRatio{ii},AllLaminaRatioError{ii},'o');
        ylabel('Lamina ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        saveas(gcf, ['figures_Lamina/SurfaceRatio and LaminaRatio for Group ' num2str(ii) '.jpg'])
        
        % plot the correlation of surface ratio with Lamina ratio for just the
        % A scores or the B scores
        figure(73)
        subplot(1,2,1)
        A = find(AllCoeff{ii}>0);
        plot(AllSurfaceRatio{ii}(A),AllLaminaRatio{ii}(A)','*');
        xlabel('Surface ratio of A TAD');
        ylabel('Lamina ratio')
        [R, P] = corr(AllSurfaceRatio{ii}(A)',AllLaminaRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(1,2,2)
        B = find(AllCoeff{ii}<0);
        plot(AllSurfaceRatio{ii}(B),AllLaminaRatio{ii}(B)','*');
        xlabel('Surface ratio of B TAD');
        ylabel('Lamina ratio')
        [R, P] = corr(AllSurfaceRatio{ii}(B)',AllLaminaRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        saveas(gcf, ['figures_Lamina/SR and LR for Group ' num2str(ii) '.jpg'])
    end
end



%% plot the nucleolar ratio with compartment score or surface ratio
for ii = 1:max(CellTypeID_new)
    if ~isempty(AllCoeff{ii}) && ~isempty(AllNucleolarRatio{ii})
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
        [R, P] = corr(AllCoeff{ii},AllNucleolarRatio{ii}', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R = ' num2str(R) ', p = ' num2str(P)]);
        
        subplot(2,1,2)
        errorbar(1:TotalNumTADs,AllNucleolarRatio{ii},AllNucleolarRatioError{ii},'o');
        ylabel('Nucelolar ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        saveas(gcf, ['figures_Nucleolus/AB score and NucelolarRatio for Group ' num2str(ii) '.jpg'])
        
        % plot the correlation of AB score with nucleolar ratio for just the
        % A scores or the B scores
        figure(71)
        subplot(2,2,1)
        A = find(AllCoeff{ii}>0);
        plot(AllCoeff{ii}(A),AllNucleolarRatio{ii}(A)','*');
        xlabel('Compartment score of A TAD');
        ylabel('Nucleolar ratio')
        [R, P] = corr(AllCoeff{ii}(A),AllNucleolarRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(2,2,2)
        B = find(AllCoeff{ii}<0);
        plot(AllCoeff{ii}(B),AllNucleolarRatio{ii}(B)','*');
        xlabel('Compartment score of B TAD');
        ylabel('Nucleolar ratio')
        [R, P] = corr(AllCoeff{ii}(B),AllNucleolarRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        % plot the correlation between surface ratio and gene density for
        % just the A scores or the B scores
        subplot(2,2,3)
        A = find(AllCoeff{ii}>0);
        plot(GeneDensity(A),AllNucleolarRatio{ii}(A)','*');
        xlabel('Gene density in A TAD');
        ylabel('Nucleolar ratio');
        [R, P] = corr(GeneDensity(A)',AllNucleolarRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(2,2,4)
        B = find(AllCoeff{ii}<0);
        plot(GeneDensity(B),AllNucleolarRatio{ii}(B)','*');
        xlabel('Gene density in B TAD');
        ylabel('Nucleolar ratio');
        [R, P] = corr(GeneDensity(B)',AllNucleolarRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        saveas(gcf, ['figures_Nucleolus/AB score NR and gene density for Group ' num2str(ii) '.jpg'])
    end

    
    if ~isempty(AllSurfaceRatio{ii}) && ~isempty(AllNucleolarRatio{ii}) && ~isempty(AllCoeff{ii})
        figure(72)
        subplot(2,1,1)
        errorbar(1:TotalNumTADs,AllSurfaceRatio{ii},AllSurfaceRatioError{ii},'o');
        ylabel('Surface ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        [R, P] = corr(AllSurfaceRatio{ii}',AllNucleolarRatio{ii}', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R = ' num2str(R) ', p = ' num2str(P)]);
        
        subplot(2,1,2)
        errorbar(1:TotalNumTADs,AllNucleolarRatio{ii},AllNucleolarRatioError{ii},'o');
        ylabel('Nucelolar ratio (95% CI)');
        xlabel('ID number of imaged TADs');
        xlim([0 TotalNumTADs+1])
        saveas(gcf, ['figures_Nucleolus/SurfaceRatio and NucelolarRatio for Group ' num2str(ii) '.jpg'])
        
        % plot the correlation of surface ratio with nucleolar ratio for just the
        % A scores or the B scores
        figure(73)
        subplot(1,2,1)
        A = find(AllCoeff{ii}>0);
        plot(AllSurfaceRatio{ii}(A),AllNucleolarRatio{ii}(A)','*');
        xlabel('Surface ratio of A TAD');
        ylabel('Nucleolar ratio')
        [R, P] = corr(AllSurfaceRatio{ii}(A)',AllNucleolarRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(1,2,2)
        B = find(AllCoeff{ii}<0);
        plot(AllSurfaceRatio{ii}(B),AllNucleolarRatio{ii}(B)','*');
        xlabel('Surface ratio of B TAD');
        ylabel('Nucleolar ratio')
        [R, P] = corr(AllSurfaceRatio{ii}(B)',AllNucleolarRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        saveas(gcf, ['figures_Nucleolus/SR and NR for Group ' num2str(ii) '.jpg'])
    end
end


%% plot the Surface ratio, compartment score, and constitutive LADs along the chromosome for each cell type
% load in the constitutive lad coordiantes from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE17051
% these are mm9 coordinates.
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
        figure(8)
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
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/AB score and SR for Group ' num2str(ii) '.jpg'])
        
        % plot the correlation of AB score with surface ratio for just the
        % A scores or the B scores
        figure(80)
        subplot(2,2,1)
        A = find(AllCoeff{ii}>0);
        plot(AllCoeff{ii}(A),AllSurfaceRatio{ii}(A)','*');
        xlabel('Compartment score of A TAD');
        ylabel('Surface ratio')
        [R, P] = corr(AllCoeff{ii}(A),AllSurfaceRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(2,2,2)
        B = find(AllCoeff{ii}<0);
        plot(AllCoeff{ii}(B),AllSurfaceRatio{ii}(B)','*');
        xlabel('Compartment score of B TAD');
        ylabel('Surface ratio')
        [R, P] = corr(AllCoeff{ii}(B),AllSurfaceRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        % plot the correlation between surface ratio and gene density for
        % just the A scores or the B scores
        subplot(2,2,3)
        A = find(AllCoeff{ii}>0);
        plot(GeneDensity(A),AllSurfaceRatio{ii}(A)','*');
        xlabel('Gene density in A TAD');
        ylabel('Surface ratio');
        [R, P] = corr(GeneDensity(A)',AllSurfaceRatio{ii}(A)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        
        subplot(2,2,4)
        B = find(AllCoeff{ii}<0);
        plot(GeneDensity(B),AllSurfaceRatio{ii}(B)','*');
        xlabel('Gene density in B TAD');
        ylabel('Surface ratio');
        [R, P] = corr(GeneDensity(B)',AllSurfaceRatio{ii}(B)', 'type', 'Pearson');
        title(['Group' num2str(ii) ', R: ' num2str(R) ', p: ' num2str(P)]);
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/AB score SR and gene density for Group ' num2str(ii) '.jpg'])
        
        % plot the surface ratio for TADs that overlap constitutive LADs and for all other
        % TADs.
        NonLADid = [];
        for i = 1:TotalNumTADs
            if length(find(LADid == i))==0
                NonLADid = [NonLADid i];
            end
        end
        figure(81)
        subplot(1,2,1)
        dotboxv(AllSurfaceRatio{ii}(LADid),  AllSurfaceRatio{ii}(NonLADid));
        set(gca, 'LineWidth', 2, 'FontSize', 16, 'XTick',[1, 2.5],'XTickLabel', {'cLAD', 'other TAD'});
        p = ranksum(AllSurfaceRatio{ii}(LADid), AllSurfaceRatio{ii}(NonLADid)); % Wilcoxon rank sum test
        title(['Group ' num2str(ii) ', p: ' num2str(p)])
        ylabel('Surface ratio');

        % plot the surface ratio for TADs that overlap constitutive LADs
        % and for all other TADs in the B compartment.
        BNonLADid = [];
        for i = 1:TotalNumTADs
            if length(find(LADid == i))==0 && AllCoeff{ii}(i)<0
                BNonLADid = [BNonLADid i];
            end
        end
        subplot(1,2,2)
        dotboxv(AllSurfaceRatio{ii}(LADid),  AllSurfaceRatio{ii}(BNonLADid));
        set(gca, 'LineWidth', 2, 'FontSize', 16, 'XTick',[1, 2.5],'XTickLabel', {'cLAD', 'other TAD in B'});
        p = ranksum(AllSurfaceRatio{ii}(LADid), AllSurfaceRatio{ii}(BNonLADid)); % Wilcoxon rank sum test
        title(['Group ' num2str(ii) ', p: ' num2str(p)])
        ylabel('Surface ratio');
        saveas(gcf, ['figures_TraceAnalysisByClusters_LargeScale/LAD surface ratio for Group ' num2str(ii) '.jpg'])
        close(81)
    end
end
















