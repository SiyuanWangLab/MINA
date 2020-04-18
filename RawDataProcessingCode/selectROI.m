%% ask the user to select ROIs


roiList =[];%used to record the ROIs

h=msgbox('Please select ROI/ROIs (rectangles)', 'Help', 'help', 'modal'); %this is the dialog 
uiwait(h); %wait till the user to click ok button

roi.rect = getrect; %select a rectangle on the image.

% hold on the image and plot the selected rectangle on the image.
hold on
x1=roi.rect(1);
y1=roi.rect(2);
x2=x1+roi.rect(3);
y2=y1;
x3=x2;
y3=y1+roi.rect(4);
x4=x1;
y4=y3;
xi=[x1, x2, x3, x4, x1];
yi=[y1, y2, y3, y4, y1];
hold on
plot(xi, yi, 'r-');
roi.index = 1;
nroi=1;
text(x1, y1, num2str(roi.index), 'Color', 'Blue', ...
    'FontWeight', 'Bold', 'FontSize',24); %write numbers on the image

roiList = cat(1, roiList, roi);

%ask user whether more ROIs are needed
moreROI = questdlg('Do you want to select more ROIs ?'); %ask question
%keep looping till the user select 'No'. 
while(strcmp(moreROI, 'Yes'))
    roi.rect = getrect;
    x1=roi.rect(1);
    y1=roi.rect(2);
    x2=x1+roi.rect(3);
    y2=y1;
    x3=x2;
    y3=y1+roi.rect(4);
    x4=x1;
    y4=y3;
    xi=[x1, x2, x3, x4, x1];
    yi=[y1, y2, y3, y4, y1];
    hold on
    plot(xi, yi, 'r-');
    roi.index = roi.index + 1;
    nroi= nroi + 1;
    text(x1, y1, num2str(roi.index), 'Color', 'Blue', ...
    'FontWeight', 'Bold', 'FontSize',24);
    roiList = cat(1, roiList, roi);
    moreROI = questdlg('Do you want to select more ROIs ?');
end