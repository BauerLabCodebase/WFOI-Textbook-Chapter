poolobj = gcp('nocreate'); % If no pool, do not create new one.
numcores = feature('numcores');
if isempty(poolobj)
    parpool('local',numcores-1);
end
%% Get Run and System Info
% excelFile="C:\Users\xiaodanwang\Documents\GitHub\BauerLab\MATLAB\examples\Code Modification\exampleTiffOIS+Gcamp.xlsx";
% excelRows=[15];  % Rows from Excell Database
%clear all;close all;

excelFile="X:\XW\CodeModification\CodeMeeting.xlsx";
excelRows=[4];
istransform = 0;

fRange_ISA = [0.01 0.08];
fRange_delta = [0.5 4];
runsInfo = parseRuns(excelFile,excelRows);
[row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!
runNum=start_ind_mouse';

%% Get Masks and WL
previousDate = [];
for runInd=runNum %for each mouse
    
    runInfo=runsInfo(runInd);
    currentDate = runInfo.recDate;
    
    if  ~contains(runInfo.system,'EastOIS2')
        
        if ~exist(runInfo.saveMaskFile,'file')
            raw = readtiff(runInfo.rawFile{1});
            raw = reshape(raw,size(raw,1),size(raw,2),runInfo.numCh,[]);
            if ~exist(runInfo.saveFolder)
                mkdir(runInfo.saveFolder);
            end
            if istransform
                [isbrain,xform_isbrain,I,WL,xform_WL] = getMask(raw,runInfo.darkFramesInd,runInfo.invalidFramesInd,runInfo.rgbInd);
                save(runInfo.saveMaskFile,'xform_WL','isbrain','I','xform_isbrain','WL','-v7.3');
            else
                [isbrain,WL] = getMask_noLandmark(raw,runInfo.darkFramesInd,runInfo.invalidFramesInd,runInfo.rgbInd);               
                save(runInfo.saveMaskFile,'isbrain','WL','-v7.3');
            end
        end
        
    else     %elseif?? strcmp('EastOIS2',runInfo.system)
        
        if ~exist(runInfo.saveMaskFile,'file')
            load(fullfile(runInfo.rawFile{1}))
            %%% Check Dropped Frame. If there is, assume it to be the 80th
            %%% dark frame
            sessionInfo = sysInfo(runInfo.system);
            if sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end)),'all')/sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end-1)),'all')>5
                numCh = size(raw_unregistered,3);
                raw_unregistered = reshape(raw_unregistered,128,128,[]);
                raw_unregistered(:,:,2:end) = raw_unregistered(:,:,1:end-1);
                raw_unregistered = reshape(raw_unregistered,128,128,numCh,[]);
            end
            
            firstFrame_cam1  = squeeze(raw_unregistered(:,:,sessionInfo.rgb(2),runInfo.darkFramesInd(end)+1));
            
            clear raw_unregistered
            load(fullfile(runInfo.rawFile{2}))
            if sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end)),'all')/sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end-1)),'all')>5
                numCh = size(raw_unregistered,3);
                raw_unregistered = reshape(raw_unregistered,128,128,[]);
                raw_unregistered(:,:,2:end) = raw_unregistered(:,:,1:end-1);
                raw_unregistered = reshape(raw_unregistered,128,128,numCh,[]);
            end
            
            firstFrame_cam2  = squeeze(raw_unregistered(:,:,sessionInfo.rgb(1),runInfo.darkFramesInd(end)+1));
            
            clear raw_unregistered
            
            if ~strcmp(previousDate,currentDate)
                [mytform,fixed_cam1,registered_cam2] = getTransformation(firstFrame_cam1,firstFrame_cam2);
                if ~exist(runInfo.saveFolder)
                    mkdir(runInfo.saveFolder)
                end
                save(fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-tform.mat')),'mytform')
            end
            load(fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-tform.mat')),'mytform')
            previousDate = currentDate;
            fixed = firstFrame_cam1./max(max(firstFrame_cam1));
            unregistered = firstFrame_cam2./max(max(firstFrame_cam2));
            registered = imwarp(unregistered, mytform,'OutputView',imref2d(size(unregistered)));
            %Create White Light Image
            WL = zeros(128,128,3);
            WL(:,:,1) = registered;
            WL(:,:,2) = fixed;
            WL(:,:,3) = fixed;
            if istransform
                [isbrain,xform_isbrain,I,seedcenter,WLcrop,xform_WLcrop,xform_WL] = getLandMarksandMask_Zyla(WL);
                save(runInfo.saveMaskFile,'xform_WL','isbrain','I','xform_isbrain','WL','-v7.3');
            else
                mask = roipoly(WL);
                isbrain = single(uint8(mask));
                save(runInfo.saveMaskFile,'isbrain','WL','-v7.3');
            end
            close all
            
        end
    end
end

%% Process Data
% Process multiple OIS files with the following naming convention:

% filedate: YYMMDD
% subjname: any alphanumberic convention is ok, e.g. "Mouse1"
% suffix: filename run number, e.g. "fc1"
% In Andor Solis, files are saved as "YYMMDD-MouseName-suffixN.tif"

% "directory" is where processed data will be stored. If the
% folder does not exist, it is created, and if directory is not
% specified, data will be saved in the current folder:

runNum = numel(runsInfo);

for runInd = 1:runNum
    clearvars -except runsInfo runNum runInd istransform
    runInfo=runsInfo(runInd);
    
    %If is this already processed
    if exist(runInfo.saveHbFile,'file')
        disp([runInfo.saveHbFile,' Already processed'])
        load(runInfo.saveHbFile,'xform_datahb')
        load(runInfo.saveMaskFile,'xform_isbrain')
        if ~isempty(runInfo.fluorChInd)
            load(runInfo.saveFluorFile,'xform_datafluorCorr')
        end
        if ~isempty(runInfo.FADChInd)
            load(runInfo.saveFADFile,'xform_dataFADCorr')
        end
    else
        load(runInfo.saveMaskFile); %anmol-load savemaskfile (I don't think there's a way that the saveMaskFile will not exist if the code reaches this point)
        
        disp(['Processing ', runInfo.saveHbFile ' on ', runInfo.system])
        sessionInfo = sysInfo(runInfo.system);
        
        %Processing Raw Data
        if ~contains(runInfo.system,'EastOIS2')
            raw = readtiff(runInfo.rawFile{1});
        else
            load(fullfile(runInfo.rawFile{1}))
            if sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end)),'all')/sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end-1)),'all')>5
                numCh = size(raw_unregistered,3);
                raw_unregistered = reshape(raw_unregistered,128,128,[]);
                raw_unregistered(:,:,2:end) = raw_unregistered(:,:,1:end-1);
                raw_unregistered = reshape(raw_unregistered,128,128,numCh,[]);
            end
            binnedRaw_cam1 = raw_unregistered(:,:,sessionInfo.cam1Chan,:);
            clear raw_unregistered
            load(fullfile(runInfo.rawFile{2}))
            if sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end)),'all')/sum(raw_unregistered(:,:,1,runInfo.darkFramesInd(end-1)),'all')>5
                numCh = size(raw_unregistered,3);
                raw_unregistered = reshape(raw_unregistered,128,128,[]);
                raw_unregistered(:,:,2:end) = raw_unregistered(:,:,1:end-1);
                raw_unregistered = reshape(raw_unregistered,128,128,numCh,[]);
            end
            binnedRaw_cam2= raw_unregistered(:,:,sessionInfo.cam2Chan,:);
            clear raw_unregistered
            load(fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-tform.mat')),'mytform')
            raw  = registerCam2andCombineTwoCams(binnedRaw_cam1,binnedRaw_cam2,mytform,sessionInfo.cam1Chan,sessionInfo.cam2Chan);
            save(fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session,num2str(runInfo.run))),'raw','-v7.3')
        end
        
        raw = reshape(raw,size(raw,1),size(raw,2),runInfo.numCh,[]);
        
        % Remove dark  frames
        if ~isempty(runInfo.darkFramesInd)
            runInfo.darkFramesInd(runInfo.darkFramesInd == runInfo.invalidFramesInd) = [];
            darkFrame = squeeze(mean(raw(:,:,:,runInfo.darkFramesInd),4,'omitnan'));
            raw = double(raw) - double(repmat(darkFrame,1,1,1,size(raw,4)));
            raw(:,:,:,1:runInfo.darkFramesInd(end))=[];
        else
            raw(:,:,:,runInfo.invalidFramesInd) = [];%invalidFramesInd default is 1
        end
        clear darkFrame
        
        %QC raw data.
        disp('rawQC analysis...');
        if ~isempty(runInfo.fluorChInd)
            representativeCh = runInfo.fluorChInd;
        else
            %which channel are we using as a reference
            representativeCh = runInfo.hbChInd(1);
        end
        
        % Assume no dropped frames by default
        droppedFrames = [];
        rawTime = 1:size(raw,4);
        rawTime = rawTime/runInfo.samplingRate;
        load(runInfo.saveMaskFile)
        raw = single(raw);
        qcRawFig = qcRaw(runInfo,rawTime,raw,representativeCh,isbrain,droppedFrames,runInfo.system);
        savefig(qcRawFig,runInfo.saveRawQCFig);
        saveas(qcRawFig,[runInfo.saveRawQCFig '.png']);
        close(qcRawFig);
        
        if ~contains(runInfo.system,'EastOIS2') % ?did we figure out dropped frames for system 2?
            % Perform check for corrected dropped frames and QC fixed
            disp('checking for dropped frames...');
            fixedData = false; % flag to see if data has been fixed yet or not
            [dfIndRaw, dfIndFixed, fixedRaw, fixedRawTime] = fixDroppedFrames(runInfo,raw,rawTime,fixedData);
            
            if isempty(dfIndRaw) % if no dropped frames, continue as normal
                disp('No dropped frames detected.');
                clear dfIndRaw dfIndFixed fixedRaw fixedRawTime;
            else
                disp(['Detected and correcting for ' num2str(length(dfIndRaw)) ' dropped frames!']);
                raw = fixedRaw;
                rawTime = fixedRawTime;
                droppedFrames = dfIndFixed;
                
                % Re-run rawQC after correction
                representativeCh = runInfo.hbChInd(1);
                qcRawFig =qcRaw(runInfo,rawTime,raw,representativeCh,isbrain,droppedFrames,runInfo.system);
                newQCSaveLoc = [runInfo.saveRawQCFig '_dropCorr'];
                saveas(qcRawFig,[newQCSaveLoc '.png']);
                savefig(qcRawFig,newQCSaveLoc);
                close(qcRawFig);
                % Check to make sure all dropped frames have been removed
                fixedData = true;
                [dfIndRaw, ~, ~, ~] = fixDroppedFrames(runInfo,fixedRaw,fixedRawTime,fixedData);
                if ~isempty(dfIndRaw)
                    currMouseName = [runInfo.recDate '-' runInfo.mouseName '-' runInfo.session...
                        '' num2str(runInfo.run)];
                    disp(['<<< Potential uncorrected dropped frames in ' currMouseName ' >>>']);
                end
                clear fixedRaw fixedRawTimeWL;
                
            end
        end
        %Detrending
        for i=1: size(raw,3)
            raw(:,:,i,:)=temporalDetrend(raw(:,:,i,:),isbrain); %?!fix?!
        end
        %OPTICAL PROP AND SPECTROSCOPY
        %HB
        [op, E] = getHbOpticalProperties({runInfo.lightSourceFiles{runInfo.hbChInd}});
        datahb = procOIS(raw(:,:,runInfo.hbChInd,:),op.dpf,E,isbrain); %where is xform_isbrain made?
        datahb = smoothImage(datahb,runInfo.gbox,runInfo.gsigma);
        if istransform
            xform_datahb = affineTransform(datahb,I); %error here
            xform_isbrain_matrix = repmat(xform_isbrain,1,1,size(xform_datahb,3),size(xform_datahb,4));
            xform_datahb(logical(1-xform_isbrain_matrix)) = NaN;
            save(runInfo.saveHbFile,'xform_datahb','op','E','runInfo','-v7.3')
        else
            isbrain_matrix = repmat(isbrain,1,1,size(datahb,3),size(datahb,4));
            datahb(logical(1-isbrain_matrix)) = NaN;
            save(runInfo.saveHbFile,'datahb','op','E','runInfo','-v7.3')
        end
        
        
        %Fluoro (calcium)
        if ~isempty(runInfo.fluorChInd)
            [op_in, E_in] = getHbOpticalProperties({runInfo.lightSourceFiles{runInfo.fluorChInd}});
            [op_out, E_out] = getHbOpticalProperties(runInfo.fluorFiles);
            dpIn = op_in.dpf/2;
            dpOut = op_out.dpf/2;
            
            datafluor = procFluor(squeeze(raw(:,:,runInfo.fluorChInd,:)),mean(raw(:,:,runInfo.fluorChInd,:),4,'omitnan'));
            datafluor = smoothImage(datafluor,runInfo.gbox,runInfo.gsigma);
            if istransform
                xform_datafluor = affineTransform(datafluor,I);
            end
            %Fluoro-Correction
            if  strcmp('EastOIS1_Fluor',runInfo.system) %gcamp correction
                datafluorCorr = correctHb_differentBeta(datafluor,datahb,[E_in(1) E_out(1)],[E_in(2) E_out(2)],dpIn,dpOut,0.7,0.8);
            elseif contains(runInfo.system,'EastOIS2') %jRGECGO1a correction
                datafluorCorr = correctHb_differentBeta(datafluor,datahb,[E_in(1) E_out(1)],[E_in(2) E_out(2)],dpIn,dpOut,1,1);
            end
            datafluorCorr = smoothImage(datafluorCorr,runInfo.gbox,runInfo.gsigma);
            if istransform
                xform_datafluorCorr=affineTransform(datafluorCorr,I);
                save(runInfo.saveFluorFile,'xform_datafluor','xform_datafluorCorr','op_in', 'E_in', 'op_out', 'E_out','runInfo','-v7.3')
                clear xform_datafluor xform_datafluorCorr
            else
                save(runInfo.saveFluorFile,'datafluor','datafluorCorr','op_in', 'E_in', 'op_out', 'E_out','runInfo','-v7.3')
            end
            clear datafluor datafluorCorr
        end
        
        %FAD
        if ~isempty(runInfo.FADChInd)
            [op_inFAD, E_inFAD] = getHbOpticalProperties({runInfo.lightSourceFiles{runInfo.FADChInd}});
            [op_outFAD, E_outFAD] = getHbOpticalProperties(runInfo.FADFiles);
            dpIn = op_inFAD.dpf/2;
            dpOut = op_outFAD.dpf/2;
            
            dataFAD = procFluor(squeeze(raw(:,:,runInfo.FADChInd,:)),mean(raw(:,:,runInfo.FADChInd,:),4,'omitnan'));
            dataFAD = smoothImage(dataFAD,runInfo.gbox,runInfo.gsigma);
            if istransform
                xform_dataFAD = affineTransform(dataFAD,I);
            end
            clear raw
            %FAD-Correction
            dataFADCorr = correctHb_differentBeta(dataFAD,datahb,[E_inFAD(1) E_outFAD(1)],[E_inFAD(2) E_outFAD(2)],dpIn,dpOut,1,1);
            dataFADCorr = smoothImage(dataFADCorr,runInfo.gbox,runInfo.gsigma);
            
            if istransform
                xform_dataFADCorr=affineTransform(dataFADCorr,I);
                save(runInfo.saveFADFile,'xform_dataFAD','xform_dataFADCorr','op_inFAD', 'E_inFAD', 'op_outFAD', 'E_outFAD','runInfo','-v7.3')
                clear xform_dataFAD xform_dataFADCorr
            else
                save(runInfo.saveFADFile,'dataFAD','dataFADCorr','op_inFAD', 'E_inFAD', 'op_outFAD', 'E_outFAD','runInfo','-v7.3')
            end
            clear dataFAD dataFADCorr
        end
        
        clear datahb rawTime  isbrain
        
        %Laser Frame
        if ~isempty(sessionInfo.chLaser)
            datalaser =  squeeze(raw(:,:,sessionInfo.chLaser,:));
            if istransform
                xform_datalaser = affineTransform(datalaser,I);
                xform_isbrain_matrix = repmat(xform_isbrain,1,1,size(xform_datalaser,3));
                xform_datalaser(logical(1-xform_isbrain_matrix)) = NaN;
                save(strcat(runInfo.saveFilePrefix,'-dataLaser.mat'),'xform_datalaser','runInfo','-v7.3')
            else
                isbrain_matrix = repmat(isbrain,1,1,size(datalaser,3));
                datalaser(logical(1-isbrain_matrix)) = NaN;
                save(strcat(runInfo.saveFilePrefix,'-dataLaser.mat'),'datalaser','runInfo','-v7.3')
            end
        end
        
        %Pre-FC Process:
        disp(['Done Processing ', runInfo.saveHbFile ' on ', runInfo.system])
    end
end


%% Averaging
