clear all
close all
% plot everage energy
Energy_all = [];
for n = 1:100
    load(['Traces\Trace' num2str(n) '.mat']);
    Energy_all = [Energy_all; Energy];
end

figure(1)
errorbar(1:size(Energy_all,2), mean(Energy_all,1), std(Energy_all)/(size(Energy_all,1))^0.5,'.','MarkerEdgeColor','red','CapSize',0);
xlabel('Simulation steps')
ylabel('Energy')
PlotProp
saveas(gcf, ['Energy.jpg'])
savefig('Energy.fig')

[h,p] = ttest2(Energy_all(:,end),Energy_all(:,end-10000))