function [ImageStack, InfoFile] = ReadZStack_MultiChannel(FileName, NumImage, StepInterval,TotalNumChannels, ChannelIndex)
% corrected on 201020 to exclude a frame in the averaging. In that frame
% the shutter is opening in the middle of the frame. I was previously
% averaging e.g. frames 3-6, 8-11, 13-16. Now I am averaging frames 4-6,
% 9-11, 14-16.

% updated 181204 so that it read in the correct images at each z height -
% exlude the frames during piezo stage movement
[MovieFP, InfoFile] = ReadDax([FileName, '.dax'],'startFrame', 1, 'endFrame', NumImage);
for i = 1:floor(NumImage/(StepInterval*TotalNumChannels))
    ImageStack(:,:,i) = mean(MovieFP(:,:, ((i-1)*TotalNumChannels+ChannelIndex-1)*StepInterval+4:((i-1)*TotalNumChannels+ChannelIndex)*StepInterval+1),3);
end
