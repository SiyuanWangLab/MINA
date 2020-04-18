clear all
close all
TotalNumTADs = 19;
DomainsToExclude = []; % list the domains to exclude, e.g. DomainsToExclude = [15, 2];
ChrName = '19';
mES_HiC = 'mES_raw\uij\uij.chr19'; % use this file to find chromosome total bin number
DomainStartsAndEnds = 'DomainStartsAndEnds19_5kb_19loci.txt';
load('ChosenDomains19_5kb_19loci.mat');
NFOV = 14; % number of fields of views

%%
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

%% analyze mean spatial distance matrix
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
    if exist(['Traces_SmallScale\TraceArrayRefined' FOVid '.mat'])==2
        load(['Traces_SmallScale\TraceArrayRefined' FOVid '.mat']);
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
        end
    end
end
display([num2str(n) ' copies of Chr' ChrName ' were analyzed.'])
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

figure(2)
imagesc(Mean)
colorbar
title('Mean spatial distance');
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp
axis square
savefig(['Mean spatial distance.fig']);

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

                Dis_pxl = [Dis_pxl Mean(i,j)];
            end
        end
    end
end


figure(4)
imagesc(GenomicDis)
colorbar
title('Pairwise genomic distance');
axis square
ColorMap = load('RedBlue.txt');
colormap(ColorMap/255);
PlotProp

figure(5)
plot(GenomicDis_pxl, Dis_pxl, '.', 'MarkerSize', 20);
xlabel('Genomic distance (bp)');
ylabel('Mean spatial distance (um)');
PlotProp
hold on
[x, IX] = sort(GenomicDis_pxl);
y = Dis_pxl(IX);
% fit to the fractal globule model
Equ = 'b*x^(1/3)';
StartPoint = [mean(y)/(mean(x))^(1/3)];
[f, gof] = fit(x', y', Equ, 'Start', StartPoint);
y_fit = f.b*x.^(1/3);
plot(x, y_fit, 'g-', 'LineWidth', 5);
% fit to get new scaling factor
Equ = 'b*x^s';
StartPoint = [mean(y)/(mean(x))^(1/3), 1/3];
[f, gof] = fit(x', y', Equ, 'Start', StartPoint);
y_fit = f.b*x.^f.s;
plot(x, y_fit, 'r-', 'LineWidth', 5);
hold off
ci = confint(f, 0.95);
title(['Slope = ' num2str(f.s) '+-' num2str(ci(2,2)-f.s)])
savefig(['SpatialVsGenomicDis.fig']);
