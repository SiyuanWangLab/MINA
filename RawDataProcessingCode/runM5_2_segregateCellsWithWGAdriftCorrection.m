% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

clear all
close all

ImageSize = 1536;
NFOV = 14; % number of fields of views
load('CodeBookSubPool3_190602.mat');
NumSteps = 11;
NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
FileNameForWGA = 'sequential/STORM2_31_';
WGAChannel = 2; % channel index in the composite z stack

%%
SE = strel('disk',3);% update on 190430
N = length(Codebook); % number of genes
load('BestParameterSet.mat');
iii = BestParameterSet % best MERFISH analysis result folder name
mkdir figures
n=0; %total cell count
m = 0; % total molecule count 
MolOutput{1,1} = 'GeneID';
MolOutput{1,2} = 'MoleculeX';
MolOutput{1,3} = 'MoleculeY';
MolOutput{1,4} = 'CellID';
MolOutput{1,5} = 'FieldOfView';
MolOutput{1,6} = 'CellCenterX';
MolOutput{1,7} = 'CellCenterY';

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
    
    FileName = [FileNameForWGA FOVid];
    [MovieFP, InfoFile] = ReadZStack_MultiChannel(FileName,NumImage,FramesToWait,TotalNumChannels,WGAChannel); % updated on 190717
    I = mean(MovieFP,3);

    
    I = I/max(max(I));
    figure(1)
    imshow(I)
    title('original WGA image')
    saveas(gcf, ['figures/figure1_' FOVid '.jpg'])
    % use adaptive thresholding to get more uniform intensity across filed of view
    I = double(I);
    I = I - min(min(I));
    I = I/max(max(I));
    T = adaptthresh(I, 0.1,'NeighborhoodSize',41);
    I = I./T;
    Min = quantile(I(:), 0.1);
    Max = quantile(I(:), 0.99);
    I = (I-Min)/(Max-Min);
    I(find(I<0)) = 0;
    I(find(I>1)) = 1;
    figure(2)
    imshow(I)
    title('processed WGA image')
    saveas(gcf, ['figures/figure2_' FOVid '.jpg'])
    
    % try to close the image to get more connected WGA boundaries
    se = strel('disk', 15);
    Iclose = imclose(I, se);
    % watershed based on processed WGA image

    L = watershed(Iclose);

    Lrgb = label2rgb(L, 'jet', 'w', 'shuffle');
    figure(3)
    imshow(Lrgb)
    title('Colored watershed label matrix')
    saveas(gcf, ['figures/figure3_' FOVid '.jpg'])    
    
    Lrgb_edge = mean(Lrgb,3);
    Lrgb_edge(find(Lrgb_edge<255)) = 0;
    Lrgb_edge = cat(3, Lrgb_edge, zeros(ImageSize));
    Lrgb_edge = cat(3, Lrgb_edge, zeros(ImageSize));
    figure(4)
    imshow(I/max(max(I)))
    hold on
    himage = imshow(Lrgb_edge);
    himage.AlphaData = 0.3;
    hold off
    title('Segregation result projected onto WGA image')
    saveas(gcf, ['figures/figure4_' FOVid '.jpg'])    
    
    %% find molecular positions
    for i = 1:N
        MolPositions{i} = [];
    end
    for kk = 1:NumSteps
        if exist(['results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat'])
            load(['results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat']);

            for i = 1:N
              stats = regionprops(CC_perfect(i),'centroid');
                
                NewPosition = round(cat(1,stats.Centroid));

                MolPositions{i} = [MolPositions{i}; NewPosition];
            end
        end
    end
    % cancel the drift: move the molecules to the WGA image frame
    if exist(['tformsWGA/tformWGA_' FOVid '.mat'])
        load(['tformsWGA/tformWGA_' FOVid '.mat']);
        for i = 1:N
            if ~isempty(MolPositions{i})
                [u,v] = transformPointsInverse(tform, MolPositions{i}(:,1),MolPositions{i}(:,2));
                MolPositions{i}(:,1) = u;
                MolPositions{i}(:,2) = v;
            end
        end
    else
        for i = 1:N
            MolPositions{i} = [];
        end
    end
    
    A = ones(ImageSize,ImageSize);
    A(ImageSize,ImageSize)=0;
    figure(5)
    imshow(A)
    axis equal
    colormap gray
    hold on
    for i = 1:N
        if ~isempty(MolPositions{i})
            plot(MolPositions{i}(:,1), MolPositions{i}(:,2), '.')
        end
    end
    hold off
    title('Perfect match; all molecules')

    saveas(gcf, ['figures/figure5_' FOVid '.jpg'])    
    %% generate cell list
    NumCells = max(max(L));
    for i = 1:NumCells
        n = n+1;
        CellList(n).CellID = n;
        CellList(n).PixelList = find(L==i);
        SingleCellImg = zeros(ImageSize,ImageSize);
        SingleCellImg(CellList(n).PixelList) = 1;
        stats = regionprops(SingleCellImg, 'centroid');
        SingleCellImg = imerode(SingleCellImg,SE); % update on 190430
        
        CellList(n).FOV = str2num(FOVid);
        CellList(n).Center = stats.Centroid;
        CellList(n).RNACopyNumber = zeros(N,1);
        for j = 1:N
            for k = 1:size(MolPositions{j},1)
                MolY = round(MolPositions{j}(k,2));
                MolX = round(MolPositions{j}(k,1));
                if MolY>0 && MolY<=ImageSize && MolX>0 && MolX<=ImageSize 
                    if SingleCellImg(MolY,MolX)==1
                        CellList(n).RNACopyNumber(j) =CellList(n).RNACopyNumber(j)+1;
                        m = m+1;

                        MolList(m).GeneID = j;
                        MolList(m).MoleculeX = MolX;
                        MolList(m).MoleculeY = MolY;
                        MolList(m).CellID = n;
                        MolList(m).FieldOfView = str2num(FOVid);
                        MolList(m).CellCenterX = CellList(n).Center(1);
                        MolList(m).CellCenterY = CellList(n).Center(2);

                        MolOutput{m+1,1} = num2str(j);
                        MolOutput{m+1,2} = num2str(MolX);
                        MolOutput{m+1,3} = num2str(MolY);
                        MolOutput{m+1,4} = num2str(n);
                        MolOutput{m+1,5} = FOVid;
                        MolOutput{m+1,6} = num2str(CellList(n).Center(1));
                        MolOutput{m+1,7} = num2str(CellList(n).Center(2));
                    end
                end
            end
        end
        CellList(n).TotalRNACopyNumber = sum(CellList(n).RNACopyNumber);
    end
end
%% 
% identify cells on the edges of FOV
EdgeImage = zeros(ImageSize,ImageSize);
EdgeImage(1:end,1) = 1;
EdgeImage(1:end,end) = 1;
EdgeImage(1,1:end) = 1;
EdgeImage(end,1:end) = 1;
for i = 1:n
    SingleCellImg = zeros(ImageSize,ImageSize);
    SingleCellImg(CellList(i).PixelList) = 1;
    if sum(sum(EdgeImage.*SingleCellImg))== 0 
        CellList(i).OnEdge = 0;
    else
        CellList(i).OnEdge = 1;
    end
end

%%
save('SingleCellAnalysisResults.mat','CellList','MolList')
XlsOutput{1,1} = 'CellID';
XlsOutput{1,2} = 'FieldOfView';
XlsOutput{1,3} = 'CellCenterX';
XlsOutput{1,4} = 'CellCenterY';
for i = 1:N
    XlsOutput{1,i+4} = Codebook(i).GeneShortName;
end
for i = 1:length(CellList)
    XlsOutput{i+1,1} = CellList(i).CellID;
    XlsOutput{i+1,2} = CellList(i).FOV;
    XlsOutput{i+1,3} = num2str(CellList(i).Center(1));
    XlsOutput{i+1,4} = num2str(CellList(i).Center(2));
    for j = 1:N
        XlsOutput{i+1,j+4} = num2str(CellList(i).RNACopyNumber(j));
    end
end
if exist('CellList.xlsx') == 2
    delete CellList.xlsx
end
if exist('MolList.xlsx') == 2
    delete MolList.xlsx
end
xlswrite('CellList.xlsx', XlsOutput);
xlswrite('MolList.xlsx', MolOutput);
        
        

