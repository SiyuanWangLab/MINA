% updated on 190508 to add compartment anlaysis

clear all
close all
TotalNumTADs = 19;
ImageSize = 1536; % number of pxls
NFOV = 14; % number of fields of views
UmPerPxl = 0.108;
DomainStartsAndEnds = 'DomainStartsAndEnds19_5kb_19loci.txt';
load('ChosenDomains19_5kb_19loci.mat');
DomainsToExclude = []; % list the domains to exclude, e.g. DomainsToExclude = [15, 2];
ChrName = 'chr19';


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


%%
load('FinalClusteringResults.mat');
mkdir figures_TraceAnalysisByClusters_SmallScale

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
    if exist(['Traces_SmallScale/TraceArrayRefined' FOVid '.mat'])==2 && exist(['tformsWGA/tformWGA_' FOVid '.mat'])==2
        load(['Traces_SmallScale/TraceArrayRefined' FOVid '.mat']);
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
save('FineTraceList.mat','Chr') %added on 190725 to facilitate joint analysis of multiple datasets.
% plot mean spatial distance matrix for cells in each group
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

    figure(1)
    imagesc(Mean)
    colorbar
    caxis([0 0.6])
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
    saveas(gcf, ['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance for Group ' num2str(ii) '.jpg'])
    savefig(['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance for Group ' num2str(ii) '.fig']); 

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
    figure(2)
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
    saveas(gcf, ['figures_TraceAnalysisByClusters_SmallScale/scaling for Group ' num2str(ii) '.jpg'])
    
    % calculate normalized distance matrix
    x = GenomicDis_pxl';
    y_fit = f.b*x.^f.s;
    Mean_adjust = Mean; 
    n = 0;
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
    figure(3)
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
    saveas(gcf, ['figures_TraceAnalysisByClusters_SmallScale/Normalized distance for Group ' num2str(ii) '.jpg'])
    
end

%% plot mean spatial distance matrix for all non-hep cells
Mean = [];
Std = [];
SEM = [];
NofData = [];
DisListAll = {};
ChrInGroup = [];
for ii = 1:2
    ChrInGroup = [ChrInGroup, Chr(find([Chr.CellType]==ii))];
end
for ii = 4:7 
    ChrInGroup = [ChrInGroup, Chr(find([Chr.CellType]==ii))];
end
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
figure(1)
imagesc(Mean)
colorbar
caxis([0 0.6]) 
title(['Mean spatial distance for non-hepatocytes, N=' num2str(length(ChrInGroup))]);
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
saveas(gcf, ['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance for non-hepatocytes.jpg'])
savefig(['figures_TraceAnalysisByClusters_SmallScale/Mean spatial distance for non-hepatocytes.fig']);

%% Plot the enrichment of contact probability between different loci with locus 19, in Hepatocyte-all other cell types.
ThresholdOfContact = 0.15; %um, distance shorter than this is considered contact
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






