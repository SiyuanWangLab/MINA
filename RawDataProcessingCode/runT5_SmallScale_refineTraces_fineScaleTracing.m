% updated on 191106 to generate a rectangular neighborhood arround each
% chromosome territory for the refiting.

% updated on 190717 to read multi-channel z stacks. The 488 bead images are
% Channel 3 in the 3-channel z stacks

% updated on 190314.

% This program refine the traces identified by runT4_linkTraces: For each 
% round of imaging, it looks for missing foci among the traces, and try to
% re-fit the missing foci in the chromosome territory defined by the fitted
% foci in each trace. 

% updated on 181203 to replace ReadZStack_medfilt with ReadZStack, 
% because the new fitFoci functions perform the median filtering inside.

clear all
close all
NFOV = 14; % number of fields of views
NumFoci = 19; % total number of secondary hybs
NumImage = 541; % Number of images in each dax 
FramesToWait = 5; % frames to wait at each height for each channel
TotalNumChannels = 3; % total number of channels in the multi-channel z stack
ImageSize = 1536; % number of pxls
WidthThreshMax = 4; % this threshold removes fitted "foci" that are too wide. Previous 4
WidthThreshMin = 0.5; % this threshold removes fitted "foci" that are too narrow.
AdjrsquareThreshold = 0.7;
UmPerPxl = 0.108;
StepSize = 0.2; %um
Hyb0IsBit1 = 0; % change this to 1 if hyb0 is bit1.

%%
FinalFociCount = zeros(NumFoci,1); % final count of foci in the traces
se = strel('disk',3);
% find foci that are missing in each identified chromosome territory and refit them
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
    if exist(['Traces_SmallScale\TraceArray' FOVid '.mat'])==2
        load(['Traces_SmallScale\TraceArray' FOVid '.mat']);
        load(['TracingDriftParams/DriftParams' FOVid '.mat']);        
        for i = 1:NumFoci
            if i<=20
                if Hyb0IsBit1 == 0
                    if i<10
                        FileName = ['sequential/STORM2_0' num2str(i) '_' FOVid];
                    else
                        FileName = ['sequential/STORM2_' num2str(i) '_' FOVid];
                    end
                elseif Hyb0IsBit1 == 1
                    if i-1<10
                        FileName = ['sequential/STORM2_0' num2str(i-1) '_' FOVid];
                    else
                        FileName = ['sequential/STORM2_' num2str(i-1) '_' FOVid];
                    end
                end
                Yd = Ydrift(i);
                Xd = Xdrift(i);
                Zd = Zdrift(i);
            else
                j = i;
                if Hyb0IsBit1 == 0
                    if j<10
                        FileName = ['sequential/STORM2_0' num2str(j) '_' FOVid];
                    else
                        FileName = ['sequential/STORM2_' num2str(j) '_' FOVid];
                    end
                elseif Hyb0IsBit1 == 1
                    if j-1<10
                        FileName = ['sequential/STORM2_0' num2str(j-1) '_' FOVid];
                    else
                        FileName = ['sequential/STORM2_' num2str(j-1) '_' FOVid];
                    end
                end
                Yd = Ydrift(j);
                Xd = Xdrift(j);
                Zd = Zdrift(j);
            end
            [ImageStack, InfoFile] = ReadZStack_MultiChannel(FileName,NumImage,FramesToWait,TotalNumChannels,2); % updated on 190717

            % determine if the traces are missing foci in this image
            m = 0;
            jList = [];% corrected on 190610
            for j = 1:length(TraceArray)
                if isempty(find(TraceArray{j}(:,end)==i)) % indeed missing, then refit
                    Ichr = zeros(ImageSize,ImageSize);
                    % define chromosome territory
                    r = round(TraceArray{j}(:,2)/UmPerPxl+Yd);
                    c = round(TraceArray{j}(:,1)/UmPerPxl+Xd);
                    Ind = find(r>=1 & r<=ImageSize & c>=1 & c<=ImageSize);
                    r = r(Ind);
                    c = c(Ind);
                    if i<=NumFoci

                        z = round(TraceArray{j}(:,3)/StepSize+Zd);
                    else
                        z = round(TraceArray{j}(:,3)/StepSize+Zd);
                    end
                    z = z(Ind);
                    
                    % this following part was modified on 191106 to
                    % generate a rectangular chromosome neighborhood 
                    Ichr(min(r):max(r),min(c):max(c)) = 1;


                    Ichr(1:7,:) = 0;
                    Ichr(:,1:7) = 0;
                    Ichr(end-6:end,:) = 0;
                    Ichr(:,end-6:end) = 0;
                    if length(find(Ichr == 1))>0 % corrected on 190607
                        m = m+1;
                        jList = [jList j];% corrected on 190610
                        IchrIndex{m} = find(Ichr == 1);
                        Zrange{m}(1) = min(max(1, min(z)-2),36); % corrected on 190606
                        Zrange{m}(2) = max(1,min(36, max(z)+2)); % corrected on 190606
                    end
                end
            end
            display('finished preparing Ichr')
            % fit focus in this chromosome territory
            if m>0
                [Xfit, Yfit, Zfit, Xgof, Ygof, Zgof, Intensity, Xwidth, Ywidth, Zwidth] = ...
                    fitFociInChr(ImageStack,IchrIndex,Zrange);
                display('finsihed fitting foci in Ichr')
                clear IchrIndex
                clear Zrange % corrected on 190606
            end
            n = 0;
            for j = 1:length(TraceArray)
                if length(find(jList == j))>0 % corrected on 190610
                    n = n+1;
                    if ~isnan(Xfit(n)) && Xgof(n).adjrsquare>AdjrsquareThreshold ...
                            && Ygof(n).adjrsquare>AdjrsquareThreshold ...
                            && Zgof(n).adjrsquare>AdjrsquareThreshold ...
                            && Xwidth(n)<WidthThreshMax && Ywidth(n)<WidthThreshMax ...
                            && Xwidth(n)>WidthThreshMin && Ywidth(n)>WidthThreshMin
                        if i<=NumFoci
                            TraceArray{j}(end+1,:) = [(Xfit(n)-Xd)*UmPerPxl, ...
                                (Yfit(n)-Yd)*UmPerPxl, (Zfit(n)-Zd)*StepSize, ...
                                Intensity(n), Xwidth(n), Ywidth(n), Zwidth(n), Xgof(n).adjrsquare, ...
                                Ygof(n).adjrsquare, Zgof(n).adjrsquare, i];
                        else
                            TraceArray{j}(end+1,:) = [(Xfit(n)-Xd)*UmPerPxl, ...
                                (Yfit(n)-Yd)*UmPerPxl, (Zfit(n)-Zd)*StepSize, ...
                                Intensity(n), Xwidth(n), Ywidth(n), Zwidth(n), Xgof(n).adjrsquare, ...
                                Ygof(n).adjrsquare, Zgof(n).adjrsquare, i];
                        end
                    end
                end
            end
            display(['Tried to fit ' num2str(n) ' new foci for FOV' FOVid ' hyb' num2str(i)])
        end
        % sort the order of foci in traces
        figure(200)
        for j = 1:length(TraceArray)
            [B, Ind] = sort(TraceArray{j}(:,end));
            TraceArray{j} = TraceArray{j}(Ind,:);
            % plot new traces in each FOV
            plot(TraceArray{j}(:,1), TraceArray{j}(:,2),'.-');
            hold on
            FinalFociCount(TraceArray{j}(:,end)) = FinalFociCount(TraceArray{j}(:,end))+1;
        end
        hold off
        axis equal
        axis ij
        xlabel('x (um)');
        ylabel('y (um)');
        savefig(['Traces_SmallScale\TracesRefined_' FOVid '.fig']);
        save(['Traces_SmallScale\TraceArrayRefined' FOVid '.mat'],'TraceArray');
    end
end
figure(4)
bar(FinalFociCount)
xlabel('Hyb number')
ylabel('Foci count after refinement')
savefig(['TracingFociCount_final_SmallScale.fig']);
    

