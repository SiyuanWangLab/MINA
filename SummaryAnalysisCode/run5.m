clear all
close all
% Manually merge clusters

MergeClusters = {[12 17], 13, [1 6 11], 8, 14, [2 7 10], [3 9], [4 5 15 16 18]}; 

%% merge clusters based on manual input
load('Clustering results original.mat');
CellTypeID_new = CellTypeID;
for i = 1:length(MergeClusters)
    for j = 1:length(MergeClusters{i})
        Ind = find(CellTypeID == MergeClusters{i}(j));
        CellTypeID_new(Ind) = i;
    end
end
figure(6)
ColorMap = colormap(jet);
ColorMap = ColorMap([1, ceil((1:(max(CellTypeID_new)-1))/(max(CellTypeID_new)-1)*length(ColorMap))],:);
gscatter(Y(:,1),Y(:,2),CellTypeID_new,ColorMap);
xlabel('tsne1')
ylabel('tsne2')
legend({'Other','Epithelial cell','Hepatocyte','Macrophage','Megakaryocyte','Erythroblast','Proerythroblast','Doublet'},'Location','bestoutside')
savefig('tsne with clusters (less clusters).fig')

% normalize gene expression for each cell
for i = 1:size(Matrix,1)
    Matrix_norm(i,:) = Matrix(i,:)/sum(Matrix(i,:));
end
% calculate the average profile for all cells
AverageProfile = mean(Matrix_norm, 1);
% build the profile for each cluster
Log2Fold = [];
ProfileForClusters = [];
for i = 1:max(CellTypeID_new)
    Idx = find(CellTypeID_new == i);
    ProfileForClusters = [ProfileForClusters; mean(Matrix_norm(Idx,:),1)];
end
for i = 1:size(ProfileForClusters,1)
    Log2Fold(i,:) = log2(ProfileForClusters(i,:)./AverageProfile);
end
% Reorder the Log2Fold matrix so that the order the marker genes are:
% Clusters 3 - 4 - 5 - 6 - 7 - 2 - 1
Log2Fold_new = [];
Log2Fold_new = [Log2Fold_new Log2Fold(:,11:15)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,48:54)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,16:20)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,47)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,21:25)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,36:38)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,39:43)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,6:10)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,1:5)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,26:30)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,44:46)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,31:35)];
Log2Fold_new = [Log2Fold_new Log2Fold(:,55)];
% Reorder the Log2Fold matrix so that the order of the cell types are:
% Hepatocyte - Macrophage - Megakaryocyte - Erythroblast - Proerythroblast
% - Epithelial cell - Other - Doublet
Log2Fold_new2 = [Log2Fold_new(3:7,:);Log2Fold_new(2,:);Log2Fold_new(1,:);Log2Fold_new(8,:)];

figure(7)
imagesc(Log2Fold_new2)
ColorMap2 = load('RedBlue.txt');
ColorMap2 = flipud(ColorMap2);
colormap(ColorMap2/255);
colorbar
xlabel('Marker gene ID');
ylabel('Cell type');
yticklabels({'Hepatocyte','Macrophage','Megakaryocyte','Erythroblast','Proerythroblast','Epithelial cell','Other','Doublet'})
title('Log2 enrichement of Gene expression')
Lim = max(max(Log2Fold_new2));
caxis([-Lim Lim])
savefig(['Gene expression for clusters (less clusters and reordered).fig'])
for i = 1:length(CellList)
    CellList(i).CellType = CellTypeID_new(i);
end
save('Clustering results final.mat','CellList','Y','CellTypeID_new','Matrix');
