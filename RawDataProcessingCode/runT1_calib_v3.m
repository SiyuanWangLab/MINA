% updated on 190314.

clear all
close all

% calibrate in 3D the two channels. Select 20 beads.
% generate the warping matrix from Cy5 channel (right) to the Cy3 channel (left).
NumImage = 131;
FramesToWait = 5; % frames to wait at each height
AdjrsquareThreshold = 0.9;
ImageSize = 1536; % number of pxls
% use this section for the case without dual view
FileName_left = ['calib/movie_0001_561'];
FileName_right = ['calib/movie_0001_647'];

[ImageStack_left, InfoFile] = ReadZStack(FileName_left,NumImage,FramesToWait);
[ImageStack_right, InfoFile] = ReadZStack(FileName_right,NumImage,FramesToWait);

ImageMeanLeft = mean(ImageStack_left,3);
ImageMeanRight = mean(ImageStack_right,3);

% %% use this section for the case with dual view
% FileName = ['calib/movie_0001'];
% [ImageStack, InfoFile] = ReadZStack(FileName,NumImage,FramesToWait);
% for i = 1:size(ImageStack, 3)
%     ImageStack_left(:,:,i) = ImageStack(:, 1:ImageSize, i);
%     ImageStack_right(:,:,i) = ImageStack(:, (ImageSize+1):(ImageSize*2), i);
% end
% ImageMeanLeft = mean(ImageStack_left,3);
% ImageMeanRight = mean(ImageStack_right,3);

%% The bead intensity is normalized to their own maximum. Thus it doesn't
% matter if the left and right channels have different signal intensity.
RGB(:,:,1) = ImageMeanLeft/max(max(ImageMeanLeft));
RGB(:,:,2) = ImageMeanRight/max(max(ImageMeanRight));
RGB(:,:,3) = zeros(ImageSize,ImageSize);
RGB(find(RGB>1)) = 1;
figure(200)
imagesc(RGB)
axis equal

WhetherROI = questdlg('Do you want to select ROIs ?'); %ask question
if(strcmp(WhetherROI, 'Yes'))
    selectROI;    
    save('CalibRoiList.mat','roiList');
else
    load('CalibRoiList.mat');
    for i = 1:length(roiList)
        hold on
        x1=roiList(i).rect(1);
        y1=roiList(i).rect(2);
        x2=x1+roiList(i).rect(3);
        y2=y1;
        x3=x2;
        y3=y1+roiList(i).rect(4);
        x4=x1;
        y4=y3;
        xi=[x1, x2, x3, x4, x1];
        yi=[y1, y2, y3, y4, y1];
        plot(xi, yi, 'r-');
        text(x1, y1, num2str(roiList(i).index), 'Color', 'Blue', ...
            'FontWeight', 'Bold', 'FontSize',24); %write numbers on the image
    end
end
lmol_match =[];
rmol_match =[];
for i = 1:length(roiList)
    [XfitL, YfitL, ZfitL, XgofL, YgofL, ZgofL] = fitFoci_gof(ImageStack_left, roiList(i), 2*i-1, 1);
    [XfitR, YfitR, ZfitR, XgofR, YgofR, ZgofR] = fitFoci_gof(ImageStack_right, roiList(i), 2*i, 1);
    
    if XgofL.adjrsquare>AdjrsquareThreshold && YgofL.adjrsquare>AdjrsquareThreshold && ZgofL.adjrsquare>AdjrsquareThreshold ...
            && XgofR.adjrsquare>AdjrsquareThreshold && YgofR.adjrsquare>AdjrsquareThreshold && ZgofR.adjrsquare>AdjrsquareThreshold
        lmol_match = cat(1, lmol_match, [XfitL, YfitL, ZfitL]);
        rmol_match = cat(1, rmol_match, [XfitR, YfitR, ZfitR]);
    else
        figure(1000)
        subplot(16,20,(2*i-1)*2-1)
        title('throw away')
        subplot(16,20,(2*i-1)*2)
        title('throw away')
        subplot(16,20,(2*i)*2-1)
        title('throw away')
        subplot(16,20,(2*i)*2)
        title('throw away')
    end
end
tform = cp2tform(rmol_match(:,1:2),lmol_match(:,1:2),'polynomial',3);
save('tform.mat','tform');
DeltaZ = mean(rmol_match(:,3)-lmol_match(:,3))
save('DeltaZ.mat','DeltaZ');
%%
TransImg = imtransform(RGB(:,:,2), tform, 'XData', [1 ImageSize], 'Ydata', [1 ImageSize]);
RGB(:,:,2) = TransImg;
figure(500)
imagesc(RGB)
axis equal
%%
mkdir('calibfigs')
figure(200)
savefig(['calibfigs/figs200.fig']);
figure(1000)
savefig(['calibfigs/figs1000.fig']);
figure(500)
savefig(['calibfigs/figs500.fig']);


