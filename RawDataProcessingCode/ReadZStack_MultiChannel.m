function [ImageStack, InfoFile] = ReadZStack_MultiChannel(FileName, NumImage, StepInterval,TotalNumChannels, ChannelIndex)
% updated 181204 so that it read in the correct images at each z height -
% exlude the frames during piezo stage movement
[MovieFP, InfoFile] = ReadDax([FileName, '.dax'],'startFrame', 1, 'endFrame', NumImage);
for i = 1:floor(NumImage/(StepInterval*TotalNumChannels))
    ImageStack(:,:,i) = mean(MovieFP(:,:, ((i-1)*TotalNumChannels+ChannelIndex-1)*StepInterval+3:((i-1)*TotalNumChannels+ChannelIndex)*StepInterval+1),3);
end
