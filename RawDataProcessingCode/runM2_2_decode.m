% updated on 190220 for the combined analysis the RNA MERFISH and DNA
% tracing. The RNA MERFISH uses 16 rounds of readout hybs and all readout
% hybs are in one (Cy7) fluorescent channel.

% updated 181204 to read in the correct images at each imaging height -
% excluding the frame during piezo stage movement.
clear all
close all
NumSteps = 11;
FramesToWait = 15; % 25 here means average every um of tissue images
ImageSize = 1536;
NFOV = 14; % number of fields of views
RoundsOfHybs = 16;
ParameterRange = 2:7;

% here define the bits
FileNameForBit{1} = 'sequential/STORM1_00_';
TformForBit{1} = '';
FileNameForBit{2} = 'sequential/STORM1_08_';
TformForBit{2} = 'tforms/tform_8_';
FileNameForBit{3} = 'sequential/STORM1_09_';
TformForBit{3} = 'tforms/tform_9_';
FileNameForBit{4} = 'sequential/STORM1_01_';
TformForBit{4} = 'tforms/tform_1_';
FileNameForBit{5} = 'sequential/STORM1_02_';
TformForBit{5} = 'tforms/tform_2_';
FileNameForBit{6} = 'sequential/STORM1_10_';
TformForBit{6} = 'tforms/tform_10_';
FileNameForBit{7} = 'sequential/STORM1_03_';
TformForBit{7} = 'tforms/tform_3_';
FileNameForBit{8} = 'sequential/STORM1_11_';
TformForBit{8} = 'tforms/tform_11_';
FileNameForBit{9} = 'sequential/STORM1_12_';
TformForBit{9} = 'tforms/tform_12_';
FileNameForBit{10} = 'sequential/STORM1_04_';
TformForBit{10} = 'tforms/tform_4_';
FileNameForBit{11} = 'sequential/STORM1_05_';
TformForBit{11} = 'tforms/tform_5_';
FileNameForBit{12} = 'sequential/STORM1_13_';
TformForBit{12} = 'tforms/tform_13_';
FileNameForBit{13} = 'sequential/STORM1_06_';
TformForBit{13} = 'tforms/tform_6_';
FileNameForBit{14} = 'sequential/STORM1_14_';
TformForBit{14} = 'tforms/tform_14_';
FileNameForBit{15} = 'sequential/STORM1_15_';
TformForBit{15} = 'tforms/tform_15_';
FileNameForBit{16} = 'sequential/STORM1_07_';
TformForBit{16} = 'tforms/tform_7_';

%%
SE = strel('square',3); % size of RNA foci in the identified image matrix
load('AllCodes.mat')
NumCodesAll = length(GoodCodes);
%%
for iii = ParameterRange
    load('Thresholds.mat');
    Thresholds = Thresholds/Thresholds(1)*(iii*20);
    mkdir(['results' num2str(iii)]);
    
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

        N = 0;
        for ii = 1:RoundsOfHybs-1
            if exist(['tforms/tform_' num2str(ii) '_' num2str(jj) '.mat'])
                N = N+1;
            end
        end
        if N == RoundsOfHybs-1
            for kk = 1:NumSteps
                StartFrame = (kk-1)*FramesToWait+3; %updated 181204 from +2 to +3

                EndFrame = (kk-1)*FramesToWait+FramesToWait+1; 
                ImageStack = zeros(ImageSize,ImageSize,16);
                for ii = 1:16
                    FileName1 = [FileNameForBit{ii} FOVid];
                    Thresh = Thresholds(ii);
                    [MovieFP, InfoFile] = ReadDax([FileName1, '.dax'],'startFrame', StartFrame, 'endFrame', EndFrame);
                    Image1 = mean(MovieFP,3);

                    if length(TformForBit{ii})>0 
                        load([TformForBit{ii} num2str(jj) '.mat']);
                        MeanX =tform.T(3,1);
                        MeanY =tform.T(3,2);
                        Image1 = imtranslate(Image1, [MeanX MeanY]);
                    end

                    background = imopen(Image1, strel('disk', 5));

                    Image2 = Image1-background;

                    RM = imregionalmax(Image2);

                    RM2 = imextendedmax(RM.*Image2,Thresh);

                    LPFImage = imdilate(RM2,SE);

                    ImageStack(:,:,ii) = LPFImage;         
                end
            
                UnitImageStack = ImageStack;

                for ii = 1:NumCodesAll
                    load(['CodeMatrices\CodeMatrix' num2str(ii) '.mat']);
                    Dis = (sum((UnitImageStack-CodeMatrix).^2,3)).^0.5;
                    Ind = find(Dis<=1);

                    % now count the molecules
                    BW = zeros(ImageSize,ImageSize);
                    BW(Ind) = 1;
                    CC(ii) = bwconncomp(BW,8);
                    
                    Ind = find(Dis==0);
                    BW = zeros(ImageSize,ImageSize);
                    BW(Ind) = 1;
                    CC_perfect(ii) = bwconncomp(BW,8);
                end
                save(['results' num2str(iii) '/CC_' FOVid '_' num2str(kk) '.mat'],'CC');
                save(['results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat'],'CC_perfect');
            end
        end
    end
end




