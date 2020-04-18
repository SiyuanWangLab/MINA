function [ImageStack, InfoFile] = ReadZStack(FileName, NumImage, StepInterval)
% updated 181204 so that it read in the correct images at each z height -
% exlude the frames during piezo stage movement
[MovieFP, InfoFile] = ReadDax([FileName, '.dax'],'startFrame', 1, 'endFrame', NumImage);
for i = 1:floor(NumImage/StepInterval)
    ImageStack(:,:,i) = mean(MovieFP(:,:, (i-1)*StepInterval+3:i*StepInterval+1),3);
end
