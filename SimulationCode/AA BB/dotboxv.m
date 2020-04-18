function yy = dotboxv(group1,group2)
SpreadWidth = 0.5;
ErrorBarWidth = 0.6;
for i = 1:length(group1)
    x1(i) = SpreadWidth*(rand-0.5)+1;
end
for i = 1:length(group2)
    x2(i) = SpreadWidth*(rand-0.5)+2.5;
end

X2 = 2.5;
X1 = 1;
BoxWidth = 0.8;
y1 = quantile(group1, [0.25 0.5 0.75]);
y2 = quantile(group2, [0.25 0.5 0.75]);
rectangle('Position',[X1-BoxWidth/2, y1(1), BoxWidth,y1(3)-y1(1)], 'FaceColor', [153, 204, 255]/255,'EdgeColor',[0 128 255]/255);
rectangle('Position',[X2-BoxWidth/2, y2(1), BoxWidth,y2(3)-y2(1)], 'FaceColor', [153, 204, 255]/255,'EdgeColor',[0 128 255]/255);
hold on

scatter_patches(x1, group1, 10, 'k','EdgeColor', 'none');
% hold on
scatter_patches(x2, group2, 10, 'k','EdgeColor', 'none');
alpha(0.35);
ylim([0 1.02])
xlim([0 3])

x1 = 1;
x2 = 2.5;
BoxWidth = 0.8;
y1 = quantile(group1, [0.25 0.5 0.75]);
y2 = quantile(group2, [0.25 0.5 0.75]);
for i = 1:3
    plot([x1-BoxWidth/2 x1+BoxWidth/2], [y1(i) y1(i)],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
    plot([x2-BoxWidth/2 x2+BoxWidth/2], [y2(i) y2(i)],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
    if i == 2
        plot([x1-BoxWidth/2 x1+BoxWidth/2], [y1(i) y1(i)],'-r','LineWidth', 2);
        plot([x2-BoxWidth/2 x2+BoxWidth/2], [y2(i) y2(i)],'-r','LineWidth', 2);
    end
end
plot([x1-BoxWidth/2 x1-BoxWidth/2], [y1(1) y1(3)],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([x1+BoxWidth/2 x1+BoxWidth/2], [y1(1) y1(3)],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([x2-BoxWidth/2 x2-BoxWidth/2], [y2(1) y2(3)],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([x2+BoxWidth/2 x2+BoxWidth/2], [y2(1) y2(3)],'-k','LineWidth', 2, 'Color', [0 128 255]/255);


hold off
xlim([0.2 3.2])
ylim([0 1.1])
yy = 1;