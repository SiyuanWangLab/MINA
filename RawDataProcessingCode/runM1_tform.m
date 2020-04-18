clear all
close all
% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

% updated on 190220 for the combined analysis the RNA MERFISH and DNA
% tracing. The RNA MERFISH uses 16 rounds of readout hybs and all readout
% hybs are in one (Cy7) fluorescent channel.

% modified on 190106 to include adaptive thresholding for bead
% identification.

NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
InitialLocalMaxThresh = 100; % brightness threshold for bead identification
NFOV = 14; % number of fields of views
ImageSize = 1536; % number of pxls
RoundsOfHybs = 16;% total number of hybs

%%
mkdir tforms
if exist('ShiftList.mat')==2
    load('ShiftList.mat');
else
    ShiftList = zeros(RoundsOfHybs-1,2);
end

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
    [MovieFP, InfoFile] = ReadZStack_MultiChannel(FileName1,NumImage,FramesToWait,TotalNumChannels,3); % updated on 190717
    Image1 = mean(MovieFP,3);
    LocalMaxThresh = InitialLocalMaxThresh;
    [Xfit1, Yfit1, Xgof1, Ygof1, Intensity1] = fitMultipleFoci2D(Image1,LocalMaxThresh);
	NumLandMarks = length(Xfit1);
    while NumLandMarks<10 && LocalMaxThresh>15
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
        LocalMaxThresh = InitialLocalMaxThresh;
        [Xfit2, Yfit2, Xgof2, Ygof2, Intensity2] = fitMultipleFoci2D(Image2,LocalMaxThresh);
        NumLandMarks2 = length(Xfit2);
        while NumLandMarks2<10 && LocalMaxThresh>15
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
        
%         figure(1)
%         imagesc(Image1);
%         axis equal
%         hold on
%         plot(Xfit1, Yfit1, 'xr')
%         hold off
        
%         figure(2)
%         imagesc(Image2);
%         axis equal
%         hold on
%         plot(Xfit2, Yfit2, 'xr')
%         hold off
        
        % try to match the beads identified in later image to the beads
        % identified in the first image
        Xfit1_match = [];
        Yfit1_match = [];
        Xfit2_match = [];
        Yfit2_match = [];
        N = 0;
        for i = 1:NumLandMarks
            Xfit1_new = Xfit1(i) + ShiftList(ii,1);
            Yfit1_new = Yfit1(i) + ShiftList(ii,2);
            [m, Ind] = min(((Xfit2-Xfit1_new).^2+(Yfit2-Yfit1_new).^2).^0.5);
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
        
%         figure(3)
%         plot(Xfit1_match,Yfit1_match,'x');
%         hold on
%         plot(Xfit2_match,Yfit2_match,'o');
%         hold off
%         axis square
%         xlim([1 size(Image1,2)])
%         ylim([1 size(Image1,1)])
%         axis ij
        
        if N>=2
            MeanX = mean(Xfit1_match-Xfit2_match);
            MeanY = mean(Yfit1_match-Yfit2_match);
            A = [1, 0, 0; 0, 1, 0; MeanX, MeanY, 1];
            tform = affine2d(A);
            
%             figure(4)
%             imagesc(TransImg)
%             axis equal       
            
            [x_new, y_new] = transformPointsForward(tform, Xfit2_match',Yfit2_match');
            
%             figure(5)
%             plot(Xfit1_match,Yfit1_match,'x');
%             hold on
%             plot(x_new,y_new,'o');
%             hold off
%             axis square
%             xlim([1 size(Image1,2)])
%             ylim([1 size(Image1,1)])
%             axis ij
            
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
                
%                 figure(6)
%                 plot(Xfit1_match,Yfit1_match,'x');
%                 hold on
%                 plot(x_new,y_new,'o');
%                 hold off
%                 axis square
%                 xlim([1 size(Image1,2)])
%                 ylim([1 size(Image1,1)])
%                 axis ij

                if sum(abs(Xfit1_match-x_new')<2 & abs(Yfit1_match-y_new')<2) == length(Ind)
                    save(['tforms/tform_' num2str(ii) '_' num2str(jj) '.mat'],'tform');
                else
                    break
                end
                
            else 
                break
            end
            
        else
            break
        end
    end
end
