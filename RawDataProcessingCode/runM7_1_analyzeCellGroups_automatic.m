clear all
close all
% updated on 190505 to perform tsne and clustering analysis on only the
% first 53 marker genes.

% updated on 190503 to generate the gene projected tsne plot here rather
% than in runM6

ImageSize = 1536; % number of pxls
load('CodeBookSubPool3_190602.mat');
MinSize = 0.25e4; % Set min cell size cutoff
MaxSize = 2e4; % Set max cell size cutoff
MinRNACount = 10; % Set min total RNA count
NFOV = 14; % number of fields of views

NicheRadius = 200; % define the radius of niche for enrichment analysis
Kvalue = 60; % k value in Louvain clustering
CroppedGeneID = [1:55];  % crop these genes in the gene list for tsne and cluster analyses
%%
% load in the genes of interest list 
load('GenesOfInterestWithOligos_FetalLiver_final.mat');

load('SingleCellAnalysisResults.mat');
% get rid of cells on the edge of each FOV
Ind = find([CellList.OnEdge] == 0);
CellList = CellList(Ind);
% get rid of cells that are too large or too small
for i = 1:length(CellList)
    Sizes(i) = length(CellList(i).PixelList);
end
Ind = find(Sizes>=MinSize & Sizes<=MaxSize);
CellList = CellList(Ind);
% get rid of cells having too few RNA counts
Ind = find([CellList.TotalRNACopyNumber]>=MinRNACount);
CellList = CellList(Ind);

% tsne analysis
k = length(CellList); % number of cells to be analyzed
N = length(CellList(1).RNACopyNumber); % number of genes to be analyzed 
Matrix = zeros(k, N);
for i = 1:k
    for j = 1:N
        Matrix(i,j) = CellList(i).RNACopyNumber(j);
    end
end

% this section was added on 190505 to perform tsne and clustering with only
% the cropped genes
Matrix_crop = Matrix(:,CroppedGeneID);
Matrix_crop_sumCol = sum(Matrix_crop,2);
Ind = find(Matrix_crop_sumCol>0);
Matrix_crop = Matrix_crop(Ind,:);
Matrix = Matrix(Ind,:);
CellList = CellList(Ind);
display(['Number of cells analyzed in Tsne = ' num2str(length(CellList))]);

% perform tsne analysis
rng default
Y = tsne(Matrix_crop,'Algorithm','exact','NumPCAComponents',50,'Perplexity', 60,'Exaggeration',20,'Distance','cosine');
figure(3)
gscatter(Y(:,1),Y(:,2))
xlabel('tsne1');
ylabel('tsne2');
title('tsne')
saveas(gcf, 'tsne analysis.jpg')

% project expressin levels of individual genes onto tsne plots
mkdir figures_tsne
xrange1 = min(Y(:,1))-(max(Y(:,1))-min(Y(:,1)))/20;
xrange2 = max(Y(:,1))+(max(Y(:,1))-min(Y(:,1)))/20;
yrange1 = min(Y(:,2))-(max(Y(:,2))-min(Y(:,2)))/20;
yrange2 = max(Y(:,2))+(max(Y(:,2))-min(Y(:,2)))/20;
for i = 1:size(Matrix,2)
    figure(31)
    scatter(Y(:,1),Y(:,2),10,Matrix(:,i),'filled')
    colormap jet
    colorbar
    xlim([xrange1, xrange2])
    ylim([yrange1, yrange2])
    xlabel('tsne1');
    ylabel('tsne2');
    title(['Gene' num2str(i) ': ' Codebook(i).GeneShortName])
    saveas(gcf, ['figures_tsne\' 'Gene' num2str(i) '_' Codebook(i).GeneShortName '.jpg'])
end

%% Louvain clustering
% build weight matrix
Idx = knnsearch(Matrix_crop, Matrix_crop, 'K' , Kvalue,'Distance','cosine');
W = zeros(size(Y,1),size(Y,1));
for i = 1:size(Idx,1)
    for j = 1:size(Idx,2)
        W(i,Idx(i,j))=1;
    end
end
W2 = zeros(size(Y,1),size(Y,1));
for i = 1:size(Idx,1)
    for j = 1:size(Idx,2)
        W2(i,Idx(i,j))=1-pdist(W([i,Idx(i,j)],:),'jaccard');
    end
end
COMTY = cluster_jl(W2);

[M, Level] = max(COMTY.MOD);
CellTypeID = COMTY.COM{Level(1)};
%%
ColorMap = colormap(jet);
ColorMap = ColorMap([1, ceil((1:(max(CellTypeID)-1))/(max(CellTypeID)-1)*length(ColorMap))],:);
figure(4)
gscatter(Y(:,1),Y(:,2),CellTypeID,ColorMap);
title('tsne')
legend('Location','bestoutside')
saveas(gcf, 'tsne analysis with clusters.jpg')

Matrix2 = [];
for i = 1:max(CellTypeID)
    Idx = find(CellTypeID == i);
    Matrix2 = [Matrix2; mean(Matrix(Idx,:),1)];
end

% normalize gene expression for each cell
for i = 1:size(Matrix,1)
    Matrix_norm(i,:) = Matrix(i,:)/sum(Matrix(i,:));
end
% calculate the average profile for all cells
AverageProfile = mean(Matrix_norm, 1);
% build the profile for each cluster
ProfileForClusters = [];
for i = 1:max(CellTypeID)
    Idx = find(CellTypeID == i);
    ProfileForClusters = [ProfileForClusters; mean(Matrix_norm(Idx,:),1)];
end
for i = 1:size(ProfileForClusters,1)
    Log2Fold(i,:) = log2(ProfileForClusters(i,:)./AverageProfile);
end
figure(5)
imagesc(Log2Fold)
hold on
ColorMap2 = load('RedBlue.txt');
ColorMap2 = flipud(ColorMap2);
colormap(ColorMap2/255);
colorbar
xlabel('Gene ID');
ylabel('Cluster ID');
title('Log2 fold change of Gene expression in each cluster')

Lim = max(max(Log2Fold));
caxis([-Lim Lim])
% draw lines to divide the 11 cell types
for i = 1:size(Log2Fold,2)
    if GenesOfInterestFinal(i).Cluster == GenesOfInterestFinal(i+1).Cluster-1
        plot([i+0.5 i+0.5], [0.5 max(CellTypeID)+0.5],'k');
    end
    if GenesOfInterestFinal(i+1).Cluster<GenesOfInterestFinal(i).Cluster
        plot([i+0.5 i+0.5], [0.5 max(CellTypeID)+0.5],'k');
        break
    end
end
hold off
saveas(gcf, 'Gene expression for cell clusters.jpg')
savefig(['Gene expression for cell clusters.fig'])

CellTypeID_new = CellTypeID;
Matrix2_new = Matrix2;

%% write xls file
A{1,1} = 'ClusterID';
for i = 1:size(Matrix2_new, 1)
    A{1,i+1} = num2str(i);
end
for i = 1:length(Codebook)
    A{i+1,1} = Codebook(i).GeneShortName;
    for j = 1: size(Matrix2_new, 1)
        A{i+1,j+1} = Matrix2_new(j,i);
    end
end
if exist('GeneExpressionOfClusters.xlsx') == 2
    delete GeneExpressionOfClusters.xlsx
end
xlswrite('GeneExpressionOfClusters.xlsx',A);

%%
for i = 1:length(CellList)
    CellList(i).CellType = CellTypeID_new(i);
end
mkdir figuresCellType
for i = 0:NFOV-1
    A = zeros(ImageSize,ImageSize);
    for j = 1:length(CellList)
        if CellList(j).FOV == i
            A(CellList(j).PixelList) = CellList(j).CellType;
        end
    end
    Lrgb = label2rgb(A, ColorMap, 'w');
    figure(8)
    imshow(Lrgb)
    title('Cell type plot')
    saveas(gcf, ['figuresCellType/figure8_' num2str(i) '.jpg'])
end

%% analyze the enrichment of pairs of cell types in neighborhood
NicheEnrichment = zeros(max(CellTypeID_new),max(CellTypeID_new));
for i = 1:length(CellList)
    for j = 1:length(CellList)
        if i~=j && CellList(i).FOV==CellList(j).FOV && (sum((CellList(i).Center-CellList(j).Center).^2))^0.5<NicheRadius
            Type1 = CellList(i).CellType;
            Type2 = CellList(j).CellType;
            if Type1 ~= Type2
            	NicheEnrichment(Type1,Type2) = NicheEnrichment(Type1,Type2)+1;
                NicheEnrichment(Type2,Type1) = NicheEnrichment(Type2,Type1)+1;
            else
                NicheEnrichment(Type1,Type2) = NicheEnrichment(Type1,Type2)+1;
            end
        end            
    end
end
Sum = 0;
for i = 1:max(CellTypeID_new)
    for j = i:max(CellTypeID_new)
        Sum = Sum + NicheEnrichment(i,j);
    end
end
NicheEnrichment = NicheEnrichment/Sum;
% caluclate the percentage of each cell type
for i = 1:max(CellTypeID_new)
    CellTypePercentage(i) = length(find([CellList.CellType]==i));
end
CellTypePercentage = CellTypePercentage/sum(CellTypePercentage);
for i = 1:max(CellTypeID_new)
    for j = 1:max(CellTypeID_new)
        if i~=j
            NicheEnrichment(i,j) = NicheEnrichment(i,j)/(2*CellTypePercentage(i)*CellTypePercentage(j));
        else
            NicheEnrichment(i,j) = NicheEnrichment(i,j)/(CellTypePercentage(i)*CellTypePercentage(j));
        end
    end
end
figure(9)
imagesc(log2(NicheEnrichment))

ColorMap2 = load('RedBlue.txt');
ColorMap2 = flipud(ColorMap2);
colormap(ColorMap2/255);
colorbar
caxis([-3 3])
xlabel('Cluster ID');
ylabel('Cluster ID');
title('Log2 enrichment of pairs of cell types in neighborhood')
saveas(gcf, 'Enrichment of pairs of cell types in neighborhood.jpg')
save('FinalClusteringResults.mat','CellList', 'CellTypeID','Matrix2','CellTypeID_new','Matrix2_new','NicheEnrichment');
