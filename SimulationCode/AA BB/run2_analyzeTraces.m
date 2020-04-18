clear all
close all
% input scores
load('G:\Miao\Summary\AllCoeff.mat')
Scores = AllCoeff{3};
figure(1)
bar(Scores);
xlabel('ID number of imaged TADs')
ylabel('Input scores')
PlotProp
saveas(gcf, ['input scores.jpg'])
savefig('input scores.fig')

for n = 1:100
    load(['Traces\Trace' num2str(n) '.mat']);
    Chr(n).x = X;
    Chr(n).y = Y;
    Chr(n).z = Z;
	Chr(n).r = ones(size(X));
end

TotalNumTADs = length(X);
Mean = [];
Std = [];
SEM = [];
NofData = [];
DisListAll = {};
ChrInGroup = Chr;
for i = 1:TotalNumTADs
    for j = 1:TotalNumTADs
        DisList = [];
        for k = 1:length(ChrInGroup)
            DisList = [DisList ((ChrInGroup(k).x(i)-ChrInGroup(k).x(j))^2+(ChrInGroup(k).y(i)-ChrInGroup(k).y(j))^2+(ChrInGroup(k).z(i)-ChrInGroup(k).z(j))^2)^0.5];
        end
        Mean(i,j) = mean(DisList);
        Std(i,j) = std(DisList);
        SEM(i,j) = std(DisList)/(length(DisList))^0.5;
        NofData(i,j) = length(DisList);
        DisListAll{i,j} = DisList;
    end
end
NChrInGroup = length(ChrInGroup);

figure(2)
imagesc(Mean)
colorbar

title(['Mean spatial distance, N=' num2str(length(ChrInGroup))]);
ylabel('ID number of imaged TADs');
xlabel('ID number of imaged TADs');
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp
axis square
saveas(gcf, ['Mean spatial distance.jpg'])
savefig(['Mean spatial distance.fig']);

%% calculate the polarization index of cell type and the randomization control
coeff = Scores;
if ~isempty(coeff)
    coeffControl = randomizeCompartments(coeff); % generate randomized compartment assignment
    for i = 1:length(ChrInGroup)
        ChrInGroup(i).PolarizationIndex = calculatePolarization(ChrInGroup(i),coeff);
        ChrInGroup(i).PolarizationIndexControl = calculatePolarization(ChrInGroup(i),coeffControl);
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
        title(['p = ' num2str(p) ', N=' num2str(length(PI))])
        ylabel('Polarization index');
        axis square
        saveas(gcf, ['Polarization index.jpg'])
        savefig(['Polarization index.fig']);
    end
    close(7)
end

% calulate the chromosome surface ratio of all TADs in this cell type (probability of TAD at surface of chromosme)
AllSurfaceIndices = [];
for i = 1:length(ChrInGroup)
    AllSurfaceIndices = [AllSurfaceIndices; calculateSurfaceIndex(ChrInGroup(i))];
end
if ~isempty(AllSurfaceIndices)
    for i = 1:TotalNumTADs
        List = AllSurfaceIndices(:,i);
        if length(find(isnan(List)~=1))~=0
            SurfaceRatio(i) = length(find(List == 1))/length(find(isnan(List)~=1));
            SurfaceRatioError(i) = 1.96*sqrt(SurfaceRatio(i)*(1-SurfaceRatio(i))/length(find(isnan(List)~=1)));
        else
            SurfaceRatio(i) = NaN;
            SurfaceRatioError(i) = NaN;
        end
    end
    % save the SurfaceRatio of this cell type in an array.
    AllSurfaceRatio = SurfaceRatio;
    AllSurfaceRatioError = SurfaceRatioError;
    % note the AllSurfaceIndices array saves the surface indices for each chromosome for all cell types
end
%% plot the correlation of AB score with surface ratio 
figure(80)
plot(coeff,AllSurfaceRatio','.b','MarkerSize', 20);
A = find(coeff>0);
hold on
plot(coeff(A),AllSurfaceRatio(A)','.r','MarkerSize', 20);
plot([0 0], [min(AllSurfaceRatio)-0.05 max(AllSurfaceRatio)+0.05],'--k','LineWidth',2)
% generate line fits for the A and B regions
y = AllSurfaceRatio(A)';
x = coeff(A);
X = [ones(size(x)), x];
b = regress(y,X);
x_fit = [0; max(coeff)+0.05];
X_fit = [ones(size(x_fit)), x_fit];
y_fit = X_fit*b;
plot(x_fit, y_fit,'-r','LineWidth',2)
B = find(coeff<0);
y = AllSurfaceRatio(B)';
x = coeff(B);
X = [ones(size(x)), x];
b = regress(y,X);
x_fit = [min(coeff)-0.05; 0];
X_fit = [ones(size(x_fit)), x_fit];
y_fit = X_fit*b;
plot(x_fit, y_fit,'-b','LineWidth',2)
hold off
xlabel('Compartment score');
ylabel('Surface ratio')
[RA, PA] = corr(coeff(A),AllSurfaceRatio(A)', 'type', 'Pearson');
[RB, PB] = corr(coeff(B),AllSurfaceRatio(B)', 'type', 'Pearson');
PlotProp
xlim([min(coeff)-0.05 max(coeff)+0.05])
ylim([min(AllSurfaceRatio)-0.05 max(AllSurfaceRatio)+0.05])
text(min(coeff), min(AllSurfaceRatio),['RB = ' num2str(RB) ', pB = ' num2str(PB) ', RA = ' num2str(RA) ', pA = ' num2str(PA)])
savefig(['AB score and SR.fig']);
saveas(gcf, ['AB score and SR.jpg'])


