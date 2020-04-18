function [Xfit, Yfit, Zfit] = fitFoci(ImageStack, roiList, Index, Index2)
% updated on 181203 to include median filtering before finding potential
% center postion, and then use the unfiltered image to fit the exact
% center. The reason for the change is that the previous strategy of
% fitting the center with filtered image may cause a minor shift in the
% fitted center position.

x1 = round(roiList.rect(1));
x2 = round(roiList.rect(1)+roiList.rect(3));
y1 = round(roiList.rect(2));
y2 = round(roiList.rect(2)+roiList.rect(4));

CropedImage = ImageStack(y1:y2, x1:x2,:);
CropedImage_medfilt = CropedImage; % updated on 181203
for i = 1:size(CropedImage,3) 
    CropedImage_medfilt(:,:,i) = medfilt2(CropedImage(:,:,i)); % updated on 181203
    MaxIntensAboveBG(i) = max(max(CropedImage_medfilt(:,:,i)))-mean2(CropedImage_medfilt(:,:,i)); % updated on 181203
end
[C,I] = max(MaxIntensAboveBG);

figure(1000*Index2)
subplot(8,10,Index*2-1)
imagesc(CropedImage(:,:,I));
title(Index)
colormap gray
axis equal
hold on
% imwrite(uint16(CropedImage(:,:,I)),strcat('Fig',num2str(HybID), '.tiff'),'tiff');

MaxIntensity = max(max(CropedImage_medfilt(:,:,I))); % updated on 181203
[Y0, X0] = find(CropedImage_medfilt(:,:,I) == MaxIntensity); % updated on 181203
Y0 = Y0(1);
X0 = X0(1);

Data1 = mean(ImageStack(Y0+y1-1,X0-3+x1-1:X0+3+x1-1,I),1);
Data2 = mean(ImageStack(Y0-3+y1-1:Y0+3+y1-1,X0+x1-1,I),2);

GaussEqu = 'a*exp(-(x-b)^2/2/c^2)+d';
StartPoint1 = [MaxIntensity X0 1 0];
StartPoint2 = [MaxIntensity Y0 1 0];

f1 = fit((X0-3:X0+3)', Data1', GaussEqu, 'Start', StartPoint1);
f2 = fit((Y0-3:Y0+3)', Data2, GaussEqu, 'Start', StartPoint2);

plot(f1.b, f2.b, 'x');
hold off

Zprofile = permute(CropedImage(Y0, X0, :), [3 1 2]);
if I-3<1
    Irange = 1:7;
elseif I+3> length(Zprofile)
    Irange = length(Zprofile)-6:length(Zprofile);
else Irange = (I-3:I+3);
end

Data3 = Zprofile(Irange);
StartPoint3 = [MaxIntensity I 1 0];
f3 = fit(Irange', Data3, GaussEqu, 'Start', StartPoint3);

figure(1000*Index2)
subplot(8,10,Index*2)
plot(f3, Irange', Data3);
legend('off');

Xfit = f1.b+x1-1;
Yfit = f2.b+y1-1;
Zfit = f3.b;