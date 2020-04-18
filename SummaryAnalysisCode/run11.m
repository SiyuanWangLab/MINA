clear all 
close all
%% plot the entire regions with scd2 gene and enhancers.
% use mm10 (GRCm38.p6 corrdinates from Ensembl)
% scd 2 gene: mm10: 44,293,674-44,306,864

DomainStartsAndEnds = 'DomainStartsAndEnds19_5kb_19loci.txt';
load('ChosenDomains19_5kb_19loci.mat');

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
DomainSize = DomainEnd-DomainStart+1;
%%
figure(1)
for i = 1:length(DomainStart)
    area([DomainStart(i) DomainCenter(i) DomainEnd(i)], [0 DomainSize(i)/2 0], 'EdgeColor','b','FaceColor','b');

    hold on
end

xlabel('Genomic position');


% plot scd2 gene
plot([44293674 44306864], [0.6e4 0.6e4], 'k','LineWidth',2);
% plot scd2 exons
plot([44293674 44294411], [0.6e4 0.6e4], 'k','LineWidth',8);
plot([44294904 44295186], [0.6e4 0.6e4], 'k','LineWidth',8);
plot([44298040 44298170], [0.6e4 0.6e4], 'k','LineWidth',8);
plot([44299609 44299814], [0.6e4 0.6e4], 'k','LineWidth',8);
plot([44301212 44301444], [0.6e4 0.6e4], 'k','LineWidth',8);
plot([44303001 44306864], [0.6e4 0.6e4], 'k','LineWidth',8);
% plot enhancers
plot([44251873 44252601], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44253001 44254875], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44253745 44255200], [0.6e4 0.6e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44257108 44257961], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44258201 44258800], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44259001 44259200], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44259601 44260200], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44260401 44260600], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);
plot([44273801 44274400], [0.4e4 0.4e4], 'Color', [230 159 0]/255,'LineWidth',8);

% plot promoters
plot([44285002 44304599], [0.4e4 0.4e4], 'Color', [1 0.5 0.5],'LineWidth',8);
plot([44292400 44299801], [0.4e4 0.4e4], 'r','LineWidth',8);


% plot promoter flanking region
plot([44263802 44268402], [0.4e4 0.4e4], 'Color', [1 0.5 0.5],'LineWidth',8);
plot([44244801 44250199], [0.4e4 0.4e4], 'Color', [1 0.5 0.5],'LineWidth',8);
plot([44237001 44238600], [0.4e4 0.4e4], 'Color', [1 0.5 0.5],'LineWidth',8);


% plot AC125101.1-201 
plot([44225796 44246570], [0.6e4 0.6e4], 'k','LineWidth',2);
% plot AC125101.1-201  exons
plot([44246570	44246481], [0.6e4 0.6e4], 'k','LineWidth',8);
plot([44227973	44225796], [0.6e4 0.6e4], 'k','LineWidth',8);

ylim([-1000 8000])
xlim([DomainStart(6)  44306864])
yticklabels({''})
yticks([])
box off
PlotProp


