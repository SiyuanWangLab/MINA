clear all
close all
% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

% updated on 190220 for the combined analysis the RNA MERFISH and DNA
% tracing. The RNA MERFISH uses 16 rounds of readout hybs and all readout
% hybs are in one (Cy7) fluorescent channel.

% updated on 190106 to fit only the brighted focus in the croped area

% updated on 190308 because hyb0 has no RNA signals and hyb1-16 are real
% signals.

FOVnumber = 0; % which filed of view to look at
NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
LocalMaxThresh = 30; % brightness threshold for bead identification
NFOV = 14; % number of fields of views
ImageSize = 1536; % number of pxls
RoundsOfHybs = 16;% total number of hybs

%%
jj = FOVnumber;   

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
if RoundsOfHybs<=10
    FileName1 = ['sequential/STORM2_0_' FOVid];
else
    FileName1 = ['sequential/STORM2_00_' FOVid];
end
%%
[MovieFP, InfoFile] = ReadZStack_MultiChannel(FileName1,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
Image1 = mean(MovieFP,3);
figure(1)
imagesc(Image1);
axis equal
hold on
rect = getrect;
Image1Crop = Image1(rect(2):(rect(2)+rect(4)), rect(1):(rect(1)+rect(3)));
[Xfit1, Yfit1, Xgof1, Ygof1, Intensity1] = fitMultipleFoci2D(Image1Crop,LocalMaxThresh);
[M,I] = max(Intensity1);
Xfit1=Xfit1(I);
Yfit1=Yfit1(I);

Xfit1 = Xfit1+rect(1)-1;
Yfit1 = Yfit1+rect(2)-1;
plot(Xfit1, Yfit1, 'xr');
hold off

ShiftList = [];
for ii = 1:RoundsOfHybs-1
    if RoundsOfHybs<=10
        RoundID = [num2str(ii)];
    else
        if ii<10
            RoundID = ['0' num2str(ii)];
        else
            RoundID = [num2str(ii)];
        end
    end
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
    FileName2 = ['sequential/STORM2_' RoundID '_' FOVid];
    
    [MovieFP, InfoFile] = ReadZStack_MultiChannel(FileName2,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
    Image2 = mean(MovieFP,3);
    figure(ii+1)
    imagesc(Image2);
    axis equal
    hold on
    rect = getrect;
    Image2Crop = Image2(rect(2):(rect(2)+rect(4)), rect(1):(rect(1)+rect(3)));
    [Xfit2, Yfit2, Xgof2, Ygof2, Intensity2] = fitMultipleFoci2D(Image2Crop,LocalMaxThresh);
    [M,I] = max(Intensity2);
    Xfit2=Xfit2(I);
    Yfit2=Yfit2(I);
    Xfit2 = Xfit2+rect(1)-1;
    Yfit2 = Yfit2+rect(2)-1;
    plot(Xfit2, Yfit2, 'xr');
    hold off
    ShiftList = [ShiftList; [Xfit2-Xfit1, Yfit2-Yfit1]];
end

save('ShiftList.mat','ShiftList')