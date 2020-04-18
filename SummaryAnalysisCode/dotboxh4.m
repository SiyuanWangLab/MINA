function yy = dotboxh4(group1,group2,group3,group4)
SpreadWidth = 0.5;
ErrorBarWidth = 0.6;
for i = 1:length(group1)
    x1(i) = SpreadWidth*(rand-0.5)+1;
end
for i = 1:length(group2)
    x2(i) = SpreadWidth*(rand-0.5)+2.5;
end
for i = 1:length(group3)
    x3(i) = SpreadWidth*(rand-0.5)+4;
end
for i = 1:length(group4)
    x4(i) = SpreadWidth*(rand-0.5)+5.5;
end

X1 = 1;
X2 = 2.5;
X3 = 4;
X4 = 5.5;
BoxWidth = 0.8;
y1 = quantile(group1, [0.25 0.5 0.75]);
y2 = quantile(group2, [0.25 0.5 0.75]);
y3 = quantile(group3, [0.25 0.5 0.75]);
y4 = quantile(group4, [0.25 0.5 0.75]);
rectangle('Position',[y1(1), X1-BoxWidth/2, y1(3)-y1(1),BoxWidth], 'FaceColor', [153, 204, 255]/255,'EdgeColor',[0 128 255]/255);
rectangle('Position',[y2(1), X2-BoxWidth/2, y2(3)-y2(1),BoxWidth], 'FaceColor', [153, 204, 255]/255,'EdgeColor',[0 128 255]/255);
rectangle('Position',[y3(1), X3-BoxWidth/2, y3(3)-y3(1),BoxWidth], 'FaceColor', [153, 204, 255]/255,'EdgeColor',[0 128 255]/255);
rectangle('Position',[y4(1), X4-BoxWidth/2, y4(3)-y4(1),BoxWidth], 'FaceColor', [153, 204, 255]/255,'EdgeColor',[0 128 255]/255);
hold on

scatter_patches(group1,x1, 1, 'k','EdgeColor', 'none');
hold on
scatter_patches(group2,x2, 1, 'k','EdgeColor', 'none');
hold on
scatter_patches(group3,x3, 1, 'k','EdgeColor', 'none');
hold on
scatter_patches(group4,x4, 1, 'k','EdgeColor', 'none');
alpha(0.35);
xlim([0 1.02])
ylim([0 3])

x1 = 1;
x2 = 2.5;
x3 = 4;
x4 = 5.5;
BoxWidth = 0.8;
y1 = quantile(group1, [0.25 0.5 0.75]);
y2 = quantile(group2, [0.25 0.5 0.75]);
y3 = quantile(group3, [0.25 0.5 0.75]);
y4 = quantile(group4, [0.25 0.5 0.75]);
for i = 1:3
    plot([y1(i) y1(i)],[x1-BoxWidth/2 x1+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
    plot([y2(i) y2(i)],[x2-BoxWidth/2 x2+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
    plot([y3(i) y3(i)],[x3-BoxWidth/2 x3+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
    plot([y4(i) y4(i)],[x4-BoxWidth/2 x4+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
    if i == 2
        plot([y1(i) y1(i)],[x1-BoxWidth/2 x1+BoxWidth/2], '-r','LineWidth', 2);
        plot([y2(i) y2(i)],[x2-BoxWidth/2 x2+BoxWidth/2], '-r','LineWidth', 2);
        plot([y3(i) y3(i)],[x3-BoxWidth/2 x3+BoxWidth/2], '-r','LineWidth', 2);
        plot([y4(i) y4(i)],[x4-BoxWidth/2 x4+BoxWidth/2], '-r','LineWidth', 2);
    end
end
plot([y1(1) y1(3)],[x1-BoxWidth/2 x1-BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y1(1) y1(3)],[x1+BoxWidth/2 x1+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y2(1) y2(3)],[x2-BoxWidth/2 x2-BoxWidth/2],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y2(1) y2(3)],[x2+BoxWidth/2 x2+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y3(1) y3(3)],[x3-BoxWidth/2 x3-BoxWidth/2],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y3(1) y3(3)],[x3+BoxWidth/2 x3+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y4(1) y4(3)],[x4-BoxWidth/2 x4-BoxWidth/2],'-k','LineWidth', 2, 'Color', [0 128 255]/255);
plot([y4(1) y4(3)],[x4+BoxWidth/2 x4+BoxWidth/2], '-k','LineWidth', 2, 'Color', [0 128 255]/255);

hold off
ylim([0.2 6.2])
xlim([0 1.1])
yy = 1;