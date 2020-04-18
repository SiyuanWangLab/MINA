clear all
close all
% updated on 190220 for the combined analysis the RNA MERFISH and DNA
% tracing. The RNA MERFISH uses 16 rounds of readout hybs and all readout
% hybs are in one (Cy7) fluorescent channel.

EndFrame = 26; % analyze the first um of tissue sample
StartFrame =2;
NFOV = 14; % number of fields of views used for finding thresholds
ThreshOriginal = 100;
RoundsOfHybs = 16;

% here define the bits
FileNameForBit{1} = 'sequential/STORM1_00_';
FileNameForBit{2} = 'sequential/STORM1_08_';
FileNameForBit{3} = 'sequential/STORM1_09_';
FileNameForBit{4} = 'sequential/STORM1_01_';
FileNameForBit{5} = 'sequential/STORM1_02_';
FileNameForBit{6} = 'sequential/STORM1_10_';
FileNameForBit{7} = 'sequential/STORM1_03_';
FileNameForBit{8} = 'sequential/STORM1_11_';
FileNameForBit{9} = 'sequential/STORM1_12_';
FileNameForBit{10} = 'sequential/STORM1_04_';
FileNameForBit{11} = 'sequential/STORM1_05_';
FileNameForBit{12} = 'sequential/STORM1_13_';
FileNameForBit{13} = 'sequential/STORM1_06_';
FileNameForBit{14} = 'sequential/STORM1_14_';
FileNameForBit{15} = 'sequential/STORM1_15_';
FileNameForBit{16} = 'sequential/STORM1_07_';

% load bulk sequencing results
load('CodeBookSubPool3_190602.mat');

%%
N = length(Codebook); % number of genes

FPKMs = zeros(N,1);
fid = fopen(['genes.fpkm_tracking'], 'r');
tline = fgetl(fid);
tline = fgetl(fid);
i = 0;
while ischar(tline)
    i = i+1;
    C = strsplit(tline);
    Gene(i).Name = str2num(C{4}(8:end));
    Gene(i).ShortName = C{5};
    Gene(i).FPKM = str2num(C{10});
    tline = fgetl(fid);
end
fclose(fid);

display('Gene names and FPKM loaded.');

for i = 1:N
    flag = 0;
    for j = 1:length(Gene)
        if strcmpi(Codebook(i).GeneShortName,Gene(j).ShortName)
            FPKMs(i) = FPKMs(i) + Gene(j).FPKM;
            flag = flag+1;
        end
    end
    if flag == 0
        display(['No FPKM found for ' Codebook(i).GeneShortName]);

    elseif flag > 1
        display(['More than 1 FPKM found for ' Codebook(i).GeneShortName]);

    end
end

%% calculate the expected foci number ratio in each image
ExpectedFociRatio = zeros(16,1);
for i = 1:N
    for j = 1:16
        if Codebook(i).Code(j) == '1';
            ExpectedFociRatio(j) = ExpectedFociRatio(j)+FPKMs(i);
        end
    end
end
ExpectedFociRatio =ExpectedFociRatio/ExpectedFociRatio(1);
% figure(1)
% bar(ExpectedFociRatio)

%% count foci in the first round of imaging
FociCount = 0;
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
        FileName1 = [FileNameForBit{1} FOVid];
        Thresh = ThreshOriginal;
        [MovieFP, InfoFile] = ReadDax([FileName1, '.dax'],'startFrame', StartFrame, 'endFrame', EndFrame);
        Image1 = mean(MovieFP,3);

        background = imopen(Image1, strel('disk', 5));
        Image2 = Image1-background;
        RM = imregionalmax(Image2);
        RM2 = imextendedmax(RM.*Image2,Thresh);
        I = RM2.*Image2;
        I2 = I(find(I>0));
        FociCount = FociCount+ length(I2);
    end
end
ExpectedFociNumber =  ExpectedFociRatio*FociCount;
Thresholds(1) = Thresh;

%% try to find best thresholds for Imaging Round 2-16, so that the observed foci number match the expected
for ii = 2:16
    if ExpectedFociNumber(ii) == 0
        Thresholds(ii) = nan;
    else

        Thresh = ThreshOriginal;

        FociCount = 0;
        AllFOVImages = [];
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
            for kk = 1:RoundsOfHybs-1
                if exist(['tforms/tform_' num2str(kk) '_' num2str(jj) '.mat'])
                    N = N+1;
                end
            end
            if N == RoundsOfHybs-1        
                FileName1 = [FileNameForBit{ii} FOVid];
                [MovieFP, InfoFile] = ReadDax([FileName1, '.dax'],'startFrame', StartFrame, 'endFrame', EndFrame);
                Image1 = mean(MovieFP,3);

                background = imopen(Image1, strel('disk', 5));
                Image2 = Image1-background;
                AllFOVImages = cat(3, AllFOVImages, Image2);
                RM = imregionalmax(Image2);
                RM2 = imextendedmax(RM.*Image2,Thresh);
                I = RM2.*Image2;
                I2 = I(find(I>0));
                FociCount = FociCount+ length(I2);
            end
        end
        if FociCount>ExpectedFociNumber(ii)
            DeltaThresh = 10;
            flag = 0;
        elseif FociCount<ExpectedFociNumber(ii)
            DeltaThresh = -10;
            flag = 0;
        else
            flag = 1;
            Thresholds(ii) = Thresh;
        end
        while flag == 0
            PreviousFociCount = FociCount;
            FociCount = 0;
            PreviousThresh = Thresh;
            Thresh = Thresh + DeltaThresh
            for jj = 1:size(AllFOVImages,3)
                Image2 = AllFOVImages(:,:,jj);
                RM = imregionalmax(Image2);
                RM2 = imextendedmax(RM.*Image2,Thresh);
                I = RM2.*Image2;
                I2 = I(find(I>0));
                FociCount = FociCount+ length(I2);
            end
            if (PreviousFociCount>ExpectedFociNumber(ii) && FociCount<=ExpectedFociNumber(ii)) || ...
                    (PreviousFociCount<ExpectedFociNumber(ii) && FociCount>=ExpectedFociNumber(ii))
                flag = 1;
                if abs(PreviousFociCount-ExpectedFociNumber(ii)) < abs(FociCount-ExpectedFociNumber(ii))
                    Thresholds(ii) = PreviousThresh;
                else
                    Thresholds(ii) = Thresh;
                end
            end
        end
    end
end
for ii = 2:16
    if isnan(Thresholds(ii))

        Thresholds(ii) = mean(Thresholds(find(Thresholds>0)));
    end
end

save('Thresholds.mat', 'Thresholds');