clear all
close all
%% load in real chromosome traces
TotalNumTADs = 50;
for n = 1:100
    load(['Traces\Trace' num2str(n) '.mat']);
    Chr(n).x = X;
    Chr(n).y = Y;
    Chr(n).z = Z;
	Chr(n).r = ones(size(X));
end

%% plot one of the chromosomes
I = 5; % chosen chromosome ID
ShowIndex = 0; % show index number of TADs if this value is 1
cmap = colormap(jet);
cmap = cmap(ceil((size(cmap,1)/TotalNumTADs)*(1:TotalNumTADs)),:);
% change color map to AB compartments
load('G:\Miao\Summary\AllCoeff.mat')
Scores = AllCoeff{3};
cmap_AB = zeros(size(cmap));
A = find(Scores>0);
B = find(Scores<0);
cmap_AB(B,3) = 1;
cmap_AB(A,1) = 1;
cmap = cmap_AB;

SphereSize = 0.2;
LineThickness = 0.04;

%%
figure(1)
b = bar(1:TotalNumTADs, ones(1,TotalNumTADs), 1 , 'FaceColor','flat');
for i = 1:TotalNumTADs
    b.CData(i,:) = cmap(i,:);
end
xlabel('ID number of imaged loci')
yticks([])
title('Color map');
set(gca, 'LineWidth', 1, 'FontSize', 12);
savefig('Color map.fig');
saveas(gcf, ['Color map.jpg'])

figure(2)
Ind = find(Chr(I).r == 1);
x = Chr(I).x(Ind);
y = Chr(I).y(Ind);
z = Chr(I).z(Ind);

% generate original sphere
[X, Y, Z]= sphere;
X = X*SphereSize;
Y = Y*SphereSize;
Z = Z*SphereSize;

for i = 1:length(Ind)
    % plot spheres
    surf(X+x(i),Y+y(i),Z+z(i),'EdgeColor','none','FaceColor',cmap(Ind(i),:));
    hold on
end
camlight
lighting gouraud
ax = gca;
% ax.Color = 'k'; % make background black
axis equal
box off
xlabel('x')
ylabel('y')
zlabel('z')

% generate original cylinder
[X,Y,Z] = cylinder;
X = X*LineThickness;
Y = Y*LineThickness;

% use a cubic spline interpolation to generate a smooth curve that links
% the points
values = csapi(1:length(x), [x y z]', 1:0.05:length(x));
for i = 1:length(Ind)-1
    range = ((i-1)*20+1):(i*20+1);
    for j = 1:20
        x1 = values(1,range(j));
        y1 = values(2,range(j));
        z1 = values(3,range(j));
        x2 = values(1,range(j+1));
        y2 = values(2,range(j+1));
        z2 = values(3,range(j+1));
        L = ((x1-x2)^2+(y1-y2)^2+(z1-z2)^2)^0.5;
        Lxy = ((x1-x2)^2+(y1-y2)^2)^0.5;
        AngleFromZ = acos((z2-z1)/L);
        if y2-y1>0
            AngleFromX = acos((x2-x1)/Lxy);
        else
            AngleFromX = -acos((x2-x1)/Lxy);
        end
        Ry = [cos(AngleFromZ) 0 sin(AngleFromZ); 0 1 0; -sin(AngleFromZ) 0 cos(AngleFromZ)];
        Rz = [cos(AngleFromX) -sin(AngleFromX) 0; sin(AngleFromX) cos(AngleFromX) 0; 0 0 1];
        XYZ1 = Rz*Ry*[X(1,:);Y(1,:);Z(1,:)];
        XYZ2 = Rz*Ry*[X(2,:);Y(2,:);Z(2,:)*L];
        Xnew = [XYZ1(1,:); XYZ2(1,:)]+x1;
        Ynew = [XYZ1(2,:); XYZ2(2,:)]+y1;
        Znew = [XYZ1(3,:); XYZ2(3,:)]+z1;
%         surf(Xnew,Ynew,Znew,'EdgeColor','none','FaceColor',cmap(Ind(i),:))
        surf(Xnew,Ynew,Znew,'EdgeColor','none','FaceColor',[0.99 0.99 0.99])
        hold on
    end
end
hold off
% show loci index
if ShowIndex == 1
    text(x,y,z,[repmat('     ',length(Ind),1), num2str((Ind))],'Color','w');
end
set(gca, 'LineWidth', 1, 'FontSize', 12);
savefig('Single trace.fig');
saveas(gcf, ['Single trace.jpg'])
