% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

% updated on 190329 
clear all
close all
ImageSize = 1536; % number of pxls
StepSize = 1394; %pxl; 150/0.1076
XSteps = 4;
YSteps = 4;
NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
NFOV = 14; % number of fields of views
FileNameForWGA = 'sequential/STORM2_31_';
WGAChannel = 2; % channel index in the composite z stack
%%
StartPxl = (ImageSize-StepSize)/2+1; %(1536-1394)/2+1=72
EndPxl = ImageSize-StartPxl+1; % ImageSize-71;
TotalImage = zeros(StepSize*YSteps, StepSize*XSteps);
X = 1-StepSize;
Y = 1;

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
    FileName1 = [FileNameForWGA FOVid];
    [MovieFP, InfoFile] = ReadZStack_MultiChannel(FileName1,NumImage,FramesToWait,TotalNumChannels,WGAChannel); % updated on 190717
    Image1 = mean(MovieFP,3);
    % use adaptive thresholding to get more uniform intensity across filed of view
    Image1 = double(Image1);
    Image1 = Image1 - min(min(Image1));
    Image1 = Image1/max(max(Image1));
    T = adaptthresh(Image1, 0.1,'NeighborhoodSize',41);
    Image1 = Image1./T;
    Min = quantile(Image1(:), 0.1);
    Max = quantile(Image1(:), 0.99);
    Image1 = (Image1-Min)/(Max-Min)*255;
    Image1(find(Image1<0)) = 0;
    Image1(find(Image1>255)) = 255;

    % adjust image oreintation so that they can be correctly stiched
    % together
    Image1(:,1:end) = Image1(:,end:-1:1);
    Image1 = Image1';
    if mod(jj,2*XSteps) == 0 || mod(jj,2*XSteps) == XSteps
        X = X + StepSize;
    elseif mod(jj,2*XSteps) >= 1 && mod(jj,2*XSteps) <= XSteps-1 
        Y = Y + StepSize;
    else
        Y = Y - StepSize;
    end
    Image1 = Image1(StartPxl:EndPxl,StartPxl:EndPxl);

    TotalImage(Y:Y+StepSize-1, X:X+StepSize-1) = Image1;
end

TotalImage = uint8(TotalImage);
figure(1)
imagesc(TotalImage);
axis equal
colormap gray
imwrite(TotalImage,'TotalImage.png')
save('TotalImage.mat','TotalImage');

