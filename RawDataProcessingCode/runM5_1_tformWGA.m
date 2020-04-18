% created on 190329 to measure the shift between the WGA images and the
% hyb0 images.

clear all
close all

NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
InitialLocalMaxThresh = 100; % brightness threshold for bead identification
NFOV = 14; % number of fields of views
ImageSize = 1536; % number of pxls
RoundsOfHybs = 16;% total number of hybs

WGAbeadFileHeader = 'sequential/STORM2_31_';

%%
mkdir tformsWGA

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
    if RoundsOfHybs<=10
        FileName1 = ['sequential/STORM2_0_' FOVid];
    else
        FileName1 = ['sequential/STORM2_00_' FOVid];
    end
    [MovieFP1, InfoFile] = ReadZStack_MultiChannel(FileName1,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
    Image1 = mean(MovieFP1,3);
    LocalMaxThresh = InitialLocalMaxThresh;

    [Xfit1, Yfit1, Xgof1, Ygof1, Intensity1] = fitMultipleFoci2D(Image1,LocalMaxThresh);
	NumLandMarks = length(Xfit1);
    while NumLandMarks<10 && LocalMaxThresh>10
        LocalMaxThresh = LocalMaxThresh-10;
        [Xfit1, Yfit1, Xgof1, Ygof1, Intensity1] = fitMultipleFoci2D(Image1,LocalMaxThresh);
        NumLandMarks = length(Xfit1);
    end
    display([num2str(NumLandMarks) ' beads identified.']);

    % find the minimum distance to each bead (from all other beads): MinDisArray
	MinDisMatrix = ones(NumLandMarks,NumLandMarks)*ImageSize;
	for i = 1:NumLandMarks
		for j = 1:NumLandMarks
            if i~=j
    			MinDisMatrix(i,j) = ((Xfit1(i)-Xfit1(j))^2+(Yfit1(i)-Yfit1(j))^2)^0.5;
            end
		end
	end
	MinDisArray = min(MinDisMatrix);

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
    FileName2 = [WGAbeadFileHeader FOVid];
        
    [MovieFP2, InfoFile] = ReadZStack_MultiChannel(FileName2,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
    Image2 = mean(MovieFP2,3);
    LocalMaxThresh = InitialLocalMaxThresh;
    [Xfit2, Yfit2, Xgof2, Ygof2, Intensity2] = fitMultipleFoci2D(Image2,LocalMaxThresh);
    NumLandMarks2 = length(Xfit2);
    while NumLandMarks2<10 && LocalMaxThresh>10
        LocalMaxThresh = LocalMaxThresh-10;
        [Xfit2, Yfit2, Xgof2, Ygof2, Intensity2] = fitMultipleFoci2D(Image2,LocalMaxThresh);
        NumLandMarks2 = length(Xfit2);
    end

    display([num2str(NumLandMarks2) ' beads identified.']);
    % find the minimum distance to each bead (from all other beads): MinDisArray
    MinDisMatrix2 = ones(NumLandMarks2,NumLandMarks2)*ImageSize;
    for i = 1:NumLandMarks2
        for j = 1:NumLandMarks2
            if i~=j
                MinDisMatrix2(i,j) = ((Xfit2(i)-Xfit2(j))^2+(Yfit2(i)-Yfit2(j))^2)^0.5;
            end
        end
    end
    MinDisArray2 = min(MinDisMatrix2);
        
    % try to match the beads identified in later image to the beads
    % identified in the first image
    Xfit1_match = [];
    Yfit1_match = [];
    Xfit2_match = [];
    Yfit2_match = [];
    N = 0;
    for i = 1:NumLandMarks
        [m, Ind] = min(((Xfit2-Xfit1(i)).^2+(Yfit2-Yfit1(i)).^2).^0.5);
        if ~isempty(m)
            if m<MinDisArray(i)/2 && m<MinDisArray2(Ind)/2 && m<60
                N = N+1;
                Xfit1_match(N) = Xfit1(i);
                Yfit1_match(N) = Yfit1(i);
                Xfit2_match(N) = Xfit2(Ind);
                Yfit2_match(N) = Yfit2(Ind);
            end
        end
    end
    display([num2str(N) ' beads initially matched.']);
    
    if N>=2
        MeanX = mean(Xfit1_match-Xfit2_match);
        MeanY = mean(Yfit1_match-Yfit2_match);
        A = [1, 0, 0; 0, 1, 0; MeanX, MeanY, 1];
        tform = affine2d(A);
        
        [x_new, y_new] = transformPointsForward(tform, Xfit2_match',Yfit2_match');
        
        MeanShiftX = median(Xfit1_match-x_new');
        MeanShiftY = median(Yfit1_match-y_new');
        Ind = find(abs(Xfit1_match-x_new'-MeanShiftX)<1 & abs(Yfit1_match-y_new'-MeanShiftY)<1);
        display([num2str(length(Ind)) ' beads finally matched.']);
        if length(Ind)>=2
            Xfit1_match = Xfit1_match(Ind);
            Yfit1_match = Yfit1_match(Ind);
            Xfit2_match = Xfit2_match(Ind);
            Yfit2_match = Yfit2_match(Ind);
            MeanX = mean(Xfit1_match-Xfit2_match);
            MeanY = mean(Yfit1_match-Yfit2_match);
            A = [1, 0, 0; 0, 1, 0; MeanX, MeanY, 1];
            tform = affine2d(A);
            [x_new, y_new] = transformPointsForward(tform, Xfit2_match',Yfit2_match');
            
            if sum(abs(Xfit1_match-x_new')<2 & abs(Yfit1_match-y_new')<2) == length(Ind)
                if jj<10
                    FOVid = ['0' num2str(jj)];
                else
                    FOVid = [num2str(jj)];                
                end
                save(['tformsWGA/tformWGA_' FOVid '.mat'],'tform');
            end
        end
    end
end
