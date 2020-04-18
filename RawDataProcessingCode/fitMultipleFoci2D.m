function  [Xfit, Yfit, Xgof, Ygof, Intensity] = fitMultipleFoci2D(Image,LocalMaxThresh)
% updated on 190106 to fit at most 100 brightest beads; also get rid of hot
% pixel by medfilt2

% Assuming the foci in the image stack are well isolated, fit the 2D
% positions of all the foci.

%threshold to find local maximum in the image stack
NoRows = size(Image, 1);
NoColumns = size(Image, 2);

I_medfilt = medfilt2(Image);
background_medfilt = imopen(I_medfilt, strel('disk', 5));
I_medfilt = I_medfilt-background_medfilt;
BW = imregionalmax(I_medfilt);

background = imopen(Image, strel('disk', 5));
I = Image-background;
BW = imextendedmax(BW.*I,LocalMaxThresh);
BW(1:7,:) = 0;
BW(end-6:end,:) = 0;
BW(:,1:7) = 0;
BW(:,end-6:end) = 0;

%identify connected ares as the individual foci
CC = bwconncomp(BW,8);
S = regionprops(CC, 'Centroid','Area');

% reduce the number of foci to the 100 brightest
if length(S)>100
    for i = 1:length(S)
        X0 = round(S(i).Centroid(1));
        Y0 = round(S(i).Centroid(2));
        CenterIntensity(i) = I(Y0, X0);
    end
    [B, Ind] = sort(CenterIntensity,'descend');
    S = S(Ind(1:100));
end

%fit the 2D positions of each focus
j = 0;

for i = 1:length(S)
    X0 = round(S(i).Centroid(1));
    Y0 = round(S(i).Centroid(2));
    CenterIntensity = Image(Y0, X0);
    
    Data1 = mean(Image(Y0,X0-7:X0+7),1);
    Data2 = mean(Image(Y0-7:Y0+7,X0),2);
    
    GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
    StartPoint1 = [CenterIntensity X0 1 0];
    StartPoint2 = [CenterIntensity Y0 1 0];
    
    [f1, gof1] = fit((X0-7:X0+7)', Data1', GaussEqu, 'Start', StartPoint1);
    [f2, gof2] = fit((Y0-7:Y0+7)', Data2, GaussEqu, 'Start', StartPoint2);
    
    if f1.b>15 && f1.b<NoColumns-14 && f2.b>15 && f2.b<NoRows-14 
        j = j+1;
        Xfit(j) = f1.b;
        Yfit(j) = f2.b;
        Xgof(j) = gof1;
        Ygof(j) = gof2;
        Intensity(j) = CenterIntensity;
    end
end
if j==0
    Xfit = [];
    Yfit = [];
    Xgof = [];
    Ygof = [];
    Intensity = [];
end

