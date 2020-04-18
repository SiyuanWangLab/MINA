function  [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity, Xwidth, Ywidth, Zwidth] = fitMultipleFoci(ImageStack,LocalMaxThresh,MaxFociNumber)
% updated on 190314.

% updated on 190221 to fit at most MaxiFociNumber of brightest beads/foci

% updated on 181203 to include median filtering before finding potential
% center postion, and then use the unfiltered image to fit the exact
% center. The reason for the change is that the previous strategy of
% fitting the center with filtered image may cause a minor shift in the
% fitted center position.

% Assuming the foci in the image stack are well isolated, fit the 3D
% positions of all the foci.

%% threshold to find local mixum in the image stack
display('start extended maxima transformation')
ImageMax = medfilt2(max(ImageStack, [],3)); % updated on 181203

%% do background subtraction
background = imopen(ImageMax, strel('disk', 4));
ImageMax = ImageMax-background;
ImageMax(find(ImageMax<0)) = 0;

BW = imextendedmax(ImageMax,LocalMaxThresh);

BW(1:7,:) = 0;
BW(end-6:end,:) = 0;
BW(:,1:7) = 0;
BW(:,end-6:end) = 0;

%% identify connected ares as the individual foci
display('start connected component analysis')
CC = bwconncomp(BW, 8);
S = regionprops(CC, 'Centroid', 'Area');
Ind = find([S.Area]<20);
S = S(Ind);

% limit the total nubmer of objects to <=MaxFociNumber
if length(S)>MaxFociNumber
    for i = 1:length(S)
        X0 = round(S(i).Centroid(1));
        Y0 = round(S(i).Centroid(2));
        CenterIntensityList(i) = ImageMax(Y0, X0);
    end
    [B, Ind] = sort(CenterIntensityList,'descend');
    S = S(Ind(1:MaxFociNumber));
end

for i = 1:length(S)
    [M, Ind] = max(permute(ImageStack(round(S(i).Centroid(2)), round(S(i).Centroid(1)),:), [3 1 2]));
    S(i).Centroid(3) = Ind;
end

display('foci center identified; start fitting')

%% fit the 3D positions of each focus
j = 0;

for i = 1:length(S)
    X0 = round(S(i).Centroid(1));
    Y0 = round(S(i).Centroid(2));
    I = round(S(i).Centroid(3));
    CenterIntensity = ImageStack(Y0, X0, I);
    
    Data1 = mean(ImageStack(Y0,X0-7:X0+7,I),1);
    Data2 = mean(ImageStack(Y0-7:Y0+7,X0,I),2);
    
    GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
    StartPoint1 = [CenterIntensity X0 1 0];
    StartPoint2 = [CenterIntensity Y0 1 0];
    
    [f1, gof1] = fit((X0-7:X0+7)', Data1', GaussEqu, 'Start', StartPoint1);
    [f2, gof2] = fit((Y0-7:Y0+7)', Data2, GaussEqu, 'Start', StartPoint2);
    
    Zprofile = permute(ImageStack(Y0, X0, :), [3 1 2]);
    if I-3<1
        Irange = 1:7;
    elseif I+3> length(Zprofile)
        Irange = length(Zprofile)-6:length(Zprofile);
    else Irange = (I-3:I+3);
    end
    
    Data3 = Zprofile(Irange);
    StartPoint3 = [CenterIntensity I 1 0];
    [f3, gof3] = fit(Irange', Data3, GaussEqu, 'Start', StartPoint3);
    
    if f1.b>X0-7 && f1.b<X0+7 && f2.b>Y0-7 && f2.b<Y0+7 && f3.b>Irange(1) && f3.b<Irange(end) %&& f1.c<2 && f2.c<2
        j = j+1;
        Xfit(j) = f1.b;
        Yfit(j) = f2.b;
        Zfit(j) = f3.b;
        Xgof(j) = gof1;
        Ygof(j) = gof2;
        Zgof(j) = gof3;

        Intensity(j) = (f1.a*f2.a*f3.a)^(1/3);
        Xwidth(j) = f1.c;
        Ywidth(j) = f2.c;
        Zwidth(j) = f3.c;
    end
end
display('done fitting')

if j==0
    Xfit = [];
    Yfit = [];
    Zfit = [];
    Xgof = [];
    Ygof = [];
    Zgof = [];
    Intensity = [];
    Xwidth = [];
    Ywidth = [];
    Zwidth = [];
end

