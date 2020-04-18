clear all
close all
% Count the single cell copy numbers, combining multiple datasets

% List datasets to include:
Folder{1} = 'G:\Miao\190622';
NFOV(1) = 12; % number of fields of views
Folder{2} = 'G:\Miao\190706';
NFOV(2) = 12; % number of fields of views
Folder{3} = 'G:\Miao\190715';
NFOV(3) = 12; % number of fields of views
Folder{4} = 'G:\Miao\190725';
NFOV(4) = 14; % number of fields of views

ImageSize = 1536; % number of pixels in x or y dimension of the image
NumSteps = 11; % number of z steps (z sections) analyzed in MERFISH
load('CodeBookSubPool3_190602.mat');
load('BestParameterSets.mat');

%%
N = length(Codebook); % number of genes
n=0; % total cell count
m=0; % total molecule count 
SE = strel('disk',3);% update on 190430
CellListAll = [];
for FolderID = 1:length(Folder)
    iii = BestParameterSet(FolderID);
    load([Folder{FolderID} '\SingleCellAnalysisResults.mat']);
    for jj = 0:NFOV(FolderID)-1
        if NFOV(FolderID)<=10
            FOVid = num2str(jj);
        elseif NFOV(FolderID)>10 && NFOV(FolderID)<=100
            if jj<10
                FOVid = ['0' num2str(jj)];
            else
                FOVid = [num2str(jj)];
            end
        elseif NFOV(FolderID)>100
            if jj<10
                FOVid = ['00' num2str(jj)];
            elseif jj<100
                FOVid = ['0' num2str(jj)];
            else
                FOVid = [num2str(jj)];
            end
        end
        % find molecular positions
        for i = 1:N
            MolPositions{i} = [];
        end
        for kk = 1:NumSteps
            if exist([Folder{FolderID} '/results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat'])
                load([Folder{FolderID} '/results' num2str(iii) '/CC_perfect_' FOVid '_' num2str(kk) '.mat']);
                for i = 1:N
                    stats = regionprops(CC_perfect(i),'centroid');
                    NewPosition = round(cat(1,stats.Centroid));
                    MolPositions{i} = [MolPositions{i}; NewPosition];
                end
            end
        end
        % cancel the drift: move the molecules to the WGA image frame
        if exist([Folder{FolderID} '/tformsWGA/tformWGA_' FOVid '.mat'])
            load([Folder{FolderID} '/tformsWGA/tformWGA_' FOVid '.mat']);
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
        % update cell list
        Ind = find([CellList.FOV] == str2num(FOVid));
        for i = 1:length(Ind)
            n = n+1;
            CellList(Ind(i)).CellID = n;
            CellList(Ind(i)).DatasetID = FolderID;
            SingleCellImg = zeros(ImageSize,ImageSize);
            SingleCellImg(CellList(Ind(i)).PixelList) = 1;
            SingleCellImg = imerode(SingleCellImg,SE); % update on 190430
            CellList(Ind(i)).RNACopyNumber = zeros(N,1);
            for j = 1:N
                for k = 1:size(MolPositions{j},1)
                    MolY = round(MolPositions{j}(k,2));
                    MolX = round(MolPositions{j}(k,1));
                    if MolY>0 && MolY<=ImageSize && MolX>0 && MolX<=ImageSize 
                        if SingleCellImg(MolY,MolX)==1
                            CellList(Ind(i)).RNACopyNumber(j) =CellList(Ind(i)).RNACopyNumber(j)+1;
                            m = m+1;
                        end
                    end
                end
            end
            CellList(Ind(i)).TotalRNACopyNumber = sum(CellList(Ind(i)).RNACopyNumber);
        end
    end
    CellListAll = [CellListAll CellList];
end
save('SingleCellAnalysisResults.mat','CellListAll')
display(['Total cell number: ' num2str(n)]);
display(['Total number of molecules in cells: ' num2str(m)]);
