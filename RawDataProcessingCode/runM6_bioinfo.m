clear all
close all
% updated on 190503 to move the tsne plotting to runM7

load('CodeBookSubPool3_190602.mat');
MinSize = 0.25e4; % Set min cell size cutoff
MaxSize = 2e4; % Set max cell size cutoff
MinRNACount = 10; % Set min total RNA count

%%
load('SingleCellAnalysisResults.mat');
% get rid of cells on the edge of each FOV
Ind = find([CellList.OnEdge] == 0);
CellList = CellList(Ind);

% get rid of cells that are too large or too small
for i = 1:length(CellList)
    Sizes(i) = length(CellList(i).PixelList);
end
figure(1)
counts = hist(Sizes, [0.125:0.25:6.875]*1e4);
hist(Sizes, [0.125:0.25:6.875]*1e4);
xlabel('Cell area (pixels)');
ylabel('Cell count');
hold on
plot([MinSize MinSize], [0 max(counts)*1.1], 'r');
plot([MaxSize MaxSize], [0 max(counts)*1.1], 'r');
ylim([0 max(counts)*1.1])
hold off
title('Cell size cutoff')
saveas(gcf, 'Cell size cutoff.jpg')
Ind = find(Sizes>=MinSize & Sizes<=MaxSize);
CellList = CellList(Ind);

% get rid of cells having too few RNA counts
figure(2)
counts = hist([CellList.TotalRNACopyNumber],[5:10:max([CellList.TotalRNACopyNumber])]);
hist([CellList.TotalRNACopyNumber],[5:10:max([CellList.TotalRNACopyNumber])]);
Ind = find([CellList.TotalRNACopyNumber]>=MinRNACount);
CellList = CellList(Ind);
hold on
plot([MinRNACount MinRNACount], [0 max(counts)*1.1], 'r');
hold off
ylim([0 max(counts)*1.1])
xlabel('Total RNA copy number per cell');
ylabel('Cell count');
title('Cellular RNA count cutoff')
saveas(gcf, 'Cellular RNA count cutoff.jpg')
