poolobj = gcp('nocreate'); % If no pool, do not create new one.
numcores = feature('numcores');
if isempty(poolobj)
    parpool('local',numcores-1);
end
%% Get Run and System Info
% excelFile="C:\Users\xiaodanwang\Documents\GitHub\BauerLab\MATLAB\examples\Code Modification\exampleTiffOIS+Gcamp.xlsx";
% excelRows=[15];  % Rows from Excell Database
%clear all;close all;
excelFile="C:\Users\Nischal\Desktop\DOI_Database.xlsx";
excelRows=[2];

fRange_ISA = [0.01 0.08];
fRange_delta = [0.5 4];
runsInfo = parseRuns(excelFile,excelRows);
[row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!

%% Get Masks and WL
previousDate = []; %initialize the date
mouse_ind=start_ind_mouse';
for runInd=mouse_ind %for each mouse
    
    runInfo=runsInfo(runInd);
    currentDate = runInfo.recDate;
    
    if  ~contains(runInfo.system,'EastOIS2')
        
        if ~exist(runInfo.saveMaskFile,'file')
            raw = readtiff(runInfo.rawFile{1}); %partial load!
            raw = reshape(raw,size(raw,1),size(raw,2),runInfo.numCh,[]); %reshape everything
            [isbrain,xform_isbrain,I,WL,xform_WL] = getMask(raw,runInfo.darkFramesInd,runInfo.invalidFramesInd,runInfo.rgbInd);
            if ~exist(runInfo.saveFolder)
                mkdir(runInfo.saveFolder);
            end
            save(runInfo.saveMaskFile,'xform_WL','isbrain','I','xform_isbrain','WL','-v7.3');
        end
        
    else
        
        if ~exist(runInfo.saveMaskFile,'file')
            sessionInfo = sysInfo(runInfo.system);
            
            tmp=matfile(runInfo.rawFile{1});   %JPC partial load  -- has size size(test,'raw_unregistered')
            raw_unregistered=tmp.raw_unregistered(:,:,1:500*4); %JPC
            raw_unregistered=reshape(raw_unregistered,size(raw_unregistered,1),size(raw_unregistered,2),runInfo.numCh,[]);
            
            %%% NO longer checking if there is dark frame...Let's just
            %%% use the "second" frame!
            firstFrame_cam1  = squeeze(raw_unregistered(:,:,sessionInfo.rgb(2),runInfo.darkFramesInd(end)+2)); % JPC made this +2 instead of +1
            
            clear raw_unregistered tmp
            
            tmp=matfile(runInfo.rawFile{2});   %JPC partial load  -- has size size(test,'raw_unregistered')
            raw_unregistered=tmp.raw_unregistered(:,:,1:500*4); %JPC
            raw_unregistered=reshape(raw_unregistered,size(raw_unregistered,1),size(raw_unregistered,2),runInfo.numCh,[]);
            
            firstFrame_cam2  = squeeze(raw_unregistered(:,:,sessionInfo.rgb(1),runInfo.darkFramesInd(end)+2)); % JPC made this +2 instead of +1
            
            clear raw_unregistered tmp
            
            if ~strcmp(previousDate,currentDate)
                ind_recDate= strfind(runInfo.rawFile{1},runInfo.recDate);
                rawDataLoc = runInfo.rawFile{1}(1:ind_recDate(2)-1);
                if exist([rawDataLoc,runInfo.recDate,'-gridWL-cam1.mat'],'file') % use the paper grid to registration if there is one
                    load([rawDataLoc,runInfo.recDate,'-gridWL-cam1.mat'])
                    cam1 = mean(raw_unregistered(:,:,4,:),4);
                    cam1=normr(cam1);
                    clear raw_unregistered
                    load([rawDataLoc,runInfo.recDate,'-gridWL-cam2.mat'])
                    cam2 = mean(raw_unregistered(:,:,4,:),4);
                    cam2=normr(cam2); %JPC
                    clear raw_unregistered
                    [mytform,fixed_cam1,registered_cam2] = getTransformation(cam1,cam2);
                else
                    clear raw_unregistered
                    [mytform,fixed_cam1,registered_cam2] = getTransformation(firstFrame_cam1,firstFrame_cam2);%if no paper grid, use the mouse data.
                end
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
            [isbrain,xform_isbrain,I,seedcenter,WLcrop,xform_WLcrop,xform_WL] = getLandMarksandMask_Zyla(WL);
            save(runInfo.saveMaskFile,'xform_WL','isbrain','I','xform_isbrain','WL','-v7.3');
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
    clearvars -except runsInfo runNum runInd fRange_ISA fRange_delta
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
        load(runInfo.saveMaskFile);
        
        disp(['Processing ', runInfo.saveHbFile ' on ', runInfo.system])
        sessionInfo = sysInfo(runInfo.system);
        
        %READ DATA AND FIND DROPPED FRAMES, and reshape
        if ~contains(runInfo.system,'EastOIS2')
            raw = readtiff(runInfo.rawFile{1}); 
            raw=FindandFixDroppedFrames(raw,runInfo); %JPC added
        else
            load(fullfile(runInfo.rawFile{1}))
            fixed_by_channel_cam1=FindandFixDroppedFrames(raw_unregistered,runInfo);
            binnedRaw_cam1 = fixed_by_channel_cam1(:,:,sessionInfo.cam1Chan,:);
            clear raw_unregistered fixed_by_channel_cam1
            
            load(fullfile(runInfo.rawFile{2}))
            fixed_by_channel_cam2=FindandFixDroppedFrames(raw_unregistered,runInfo);
            binnedRaw_cam2= fixed_by_channel_cam2(:,:,sessionInfo.cam2Chan,:);
            clear raw_unregistered fixed_by_channel_cam2
            
            load(fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-tform.mat')),'mytform')
            raw  = registerCam2andCombineTwoCams(binnedRaw_cam1,binnedRaw_cam2,mytform,sessionInfo.cam1Chan,sessionInfo.cam2Chan);
            save(fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session,num2str(runInfo.run))),'raw','-v7.3')
        end
        
        %DARK FRAME REMOVAL        
        % Remove dark  frames (JPC changed to subtract average across all
        % darkframes
        if ~isempty(runInfo.darkFramesInd)
            runInfo.darkFramesInd(runInfo.darkFramesInd == runInfo.invalidFramesInd) = [];
            darkFrame = squeeze(mean(reshape(raw(:,:,:,runInfo.darkFramesInd),size(raw,1),size(raw,2),[]),3,'omitnan')); %JPC made it so that we took average across all 4 "channels"
            raw = double(raw) - double(repmat(darkFrame,1,1,4,size(raw,4)));
            raw(:,:,:,1:runInfo.darkFramesInd(end)+1)=[]; %JPC get rid of the dark frames, and the first-non Dark frame
            raw(:,:,:,end)=[]; %JPC get rid of last frame too
        else
            raw(:,:,:,runInfo.invalidFramesInd) = [];%invalidFramesInd default is 1
        end
        clear darkFrame
        
        % QC RAW DATA
        disp('rawQC analysis...');
        if ~isempty(runInfo.fluorChInd)
            representativeCh = runInfo.fluorChInd;
        else
            %which channel are we using as a reference
            representativeCh = runInfo.hbChInd(1);
        end
        
        rawTime = 1:size(raw,4);
        rawTime = rawTime/runInfo.samplingRate;
        load(runInfo.saveMaskFile)
        raw = single(raw);
        qcRawFig = qcRaw(runInfo,rawTime,raw,representativeCh,isbrain,runInfo.system);
        savefig(qcRawFig,runInfo.saveRawQCFig);
        saveas(qcRawFig,[runInfo.saveRawQCFig '.png']);
        close(qcRawFig);
        
        %Detrending
        for i=1: size(raw,3)
            raw(:,:,i,:)=temporalDetrend(raw(:,:,i,:),isbrain); %?!fix?!
        end
        %OPTICAL PROP AND SPECTROSCOPY
        %HB
        [op, E] = getHbOpticalProperties({runInfo.lightSourceFiles{runInfo.hbChInd}});
        datahb = procOIS(raw(:,:,runInfo.hbChInd,:),op.dpf,E,isbrain); %where is xform_isbrain made?
        datahb = smoothImage(datahb,runInfo.gbox,runInfo.gsigma);
        xform_datahb = affineTransform(datahb,I); %error here
        xform_isbrain_matrix = repmat(xform_isbrain,1,1,size(xform_datahb,3),size(xform_datahb,4));
        xform_datahb(logical(1-xform_isbrain_matrix)) = NaN;
        save(runInfo.saveHbFile,'xform_datahb','op','E','runInfo','-v7.3')
        
        
        %Fluoro (calcium)
        if ~isempty(runInfo.fluorChInd)
            [op_in, E_in] = getHbOpticalProperties({runInfo.lightSourceFiles{runInfo.fluorChInd}});
            [op_out, E_out] = getHbOpticalProperties(runInfo.fluorFiles);
            dpIn = op_in.dpf/2;
            dpOut = op_out.dpf/2;
            
            datafluor = procFluor(squeeze(raw(:,:,runInfo.fluorChInd,:)),mean(raw(:,:,runInfo.fluorChInd,:),4,'omitnan'));
            datafluor = smoothImage(datafluor,runInfo.gbox,runInfo.gsigma);
            xform_datafluor = affineTransform(datafluor,I);
            %Fluoro-Correction
            if  strcmp('EastOIS1_Fluor',runInfo.system) %gcamp correction
                datafluorCorr = correctHb_differentBeta(datafluor,datahb,[E_in(1) E_out(1)],[E_in(2) E_out(2)],dpIn,dpOut,0.7,0.8);
            elseif contains(runInfo.system,'EastOIS2') %jRGECGO1a correction
                datafluorCorr = correctHb_differentBeta(datafluor,datahb,[E_in(1) E_out(1)],[E_in(2) E_out(2)],dpIn,dpOut,1,1);
            end
            datafluorCorr = smoothImage(datafluorCorr,runInfo.gbox,runInfo.gsigma);
            xform_datafluorCorr=affineTransform(datafluorCorr,I);
            save(runInfo.saveFluorFile,'xform_datafluor','xform_datafluorCorr','op_in', 'E_in', 'op_out', 'E_out','runInfo','-v7.3')
            clear datafluor datafluorCorr xform_datafluor
        end
        
        %FAD
        if ~isempty(runInfo.FADChInd)
            [op_inFAD, E_inFAD] = getHbOpticalProperties({runInfo.lightSourceFiles{runInfo.FADChInd}});
            [op_outFAD, E_outFAD] = getHbOpticalProperties(runInfo.FADFiles);
            dpIn = op_inFAD.dpf/2;
            dpOut = op_outFAD.dpf/2;
            
            dataFAD = procFluor(squeeze(raw(:,:,runInfo.FADChInd,:)),mean(raw(:,:,runInfo.FADChInd,:),4,'omitnan'));
            dataFAD = smoothImage(dataFAD,runInfo.gbox,runInfo.gsigma);
            xform_dataFAD = affineTransform(dataFAD,I);
            clear raw
            %FAD-Correction
            dataFADCorr = correctHb_differentBeta(dataFAD,datahb,[E_inFAD(1) E_outFAD(1)],[E_inFAD(2) E_outFAD(2)],dpIn,dpOut,1,1);
            dataFADCorr = smoothImage(dataFADCorr,runInfo.gbox,runInfo.gsigma);
            xform_dataFADCorr=affineTransform(dataFADCorr,I);
            
            
            save(runInfo.saveFADFile,'xform_dataFAD','xform_dataFADCorr','op_inFAD', 'E_inFAD', 'op_outFAD', 'E_outFAD','runInfo','-v7.3')
            clear dataFAD xform_dataFAD dataFADCorr
        end
        
        clear datahb rawTime  isbrain
        
        %Laser Frame
        if ~isempty(sessionInfo.chLaser)
            datalaser =  squeeze(raw(:,:,sessionInfo.chLaser,:));
            xform_datalaser = affineTransform(datalaser,I);
            xform_isbrain_matrix = repmat(xform_isbrain,1,1,size(xform_datalaser,3));
            xform_datalaser(logical(1-xform_isbrain_matrix)) = NaN;
            save(strcat(runInfo.saveFilePrefix,'-dataLaser.mat'),'xform_datalaser','runInfo','-v7.3')
        end
    end
    
    
    
    
    %Pre-FC Process:
    if runInfo.session =="fc"
        %Initialize
        if ~isempty(runInfo.fluorChInd)
            data_full = nan(size(xform_datahb,1),size(xform_datahb,2),size(xform_datahb,4),4);
        elseif ~isempty(runInfo.FADChInd)
            data_full = nan(size(xform_datahb,1),size(xform_datahb,2),size(xform_datahb,4),5);
        else
            data_full = nan(size(xform_datahb,1),size(xform_datahb,2),size(xform_datahb,4),3);
        end
        %Shape Data
        data_full(:,:,:,1)=squeeze(xform_datahb(:,:,1,:));
        data_full(:,:,:,2)=squeeze(xform_datahb(:,:,2,:));
        data_full(:,:,:,3)=squeeze(xform_datahb(:,:,1,:))+squeeze(xform_datahb(:,:,2,:));
        if ~isempty(runInfo.fluorChInd)
            data_full(:,:,:,4) = xform_datafluorCorr;
        end
        if ~isempty(runInfo.FADChInd)
            data_full(:,:,:,5) = xform_dataFADCorr;
        end
        
        %Power Break- Can we initalize it?
        for i=1:size(data_full,4)
            [whole_spectra_map(:,:,:,i),powerMap(:,:,:,i),hz, global_sig_for(:,i),glob_sig_power(:,i)]= PowerAnalysis(squeeze(data_full(:,:,:,i)),runInfo.samplingRate,xform_isbrain); %anmol-semicolon
        end
        save(strcat(runInfo.saveFilePrefix,'-Power.mat'),'hz','whole_spectra_map','powerMap','global_sig_for','glob_sig_power','-v7.3')
        
        clear  xform_datafluorCorr xform_datahb xform_dataFADCorr ...
            whole_spectra_map powerMap global_sig_for glob_sig_power
        
        %Filter and downsample
        data_isa= nan(size(data_full));
        data_delta= nan(size(data_full));
        for ii=1:size(data_full,4)
            data_isa(:,:,:,ii)=filterData(data_full(:,:,:,ii),fRange_ISA(1),fRange_ISA(2),runInfo.samplingRate,xform_isbrain); %JPC changed filterData to have varargin for isbrain
            data_delta(:,:,:,ii)=filterData(data_full(:,:,:,ii),fRange_delta(1),fRange_delta(2),runInfo.samplingRate,xform_isbrain);
        end
        
        clear data_full
        
        %GSR
        for  ii=1:size(data_isa,4)
            data_isa(:,:,:,ii) = gsr(squeeze(data_isa(:,:,:,ii)),xform_isbrain);
            data_delta(:,:,:,ii) = gsr(squeeze(data_delta(:,:,:,ii)),xform_isbrain);
        end
        
        xform_isbrain_matrix = repmat(xform_isbrain,1,1,size(data_delta,3),size(data_delta,4));
        data_isa(logical(1-xform_isbrain_matrix)) = NaN;%
        data_delta(logical(1-xform_isbrain_matrix)) = NaN;
        %QC FC
        %ISA
        tic
        fStr = [num2str(fRange_ISA(1)) '-' num2str(fRange_ISA(2))];
        fStr(strfind(fStr,'.')) = 'p';
        ContFig2Save= {'HbO','HbT'};
        [fh, seedFC, seedFCMap,seedCenter,seedRadius,bilatFCMap]...
            = qcFC(data_isa,xform_isbrain,runInfo.Contrasts,runInfo,fRange_ISA,ContFig2Save);
        
        for contrastInd = 1:numel(ContFig2Save)
            saveFCQCFigName = [runInfo.saveFCQCFig '-' ContFig2Save{contrastInd} '-' fStr ];
            saveas(fh(contrastInd),[saveFCQCFigName '.png']);
            close(fh(contrastInd));
        end
        toc
        %Saving
        contrastName = runInfo.Contrasts;
        saveSeedFCName = [runInfo.saveFilePrefix '-seedFC-' fStr];
        save(saveSeedFCName,'contrastName','seedCenter','seedFC','seedFCMap','seedRadius','-v7.3');
        saveBilateralFCName = [runInfo.saveFilePrefix '-bilateralFC-' fStr];
        save(saveBilateralFCName,'contrastName','bilatFCMap','xform_isbrain','-v7.3');
        
        
        %Delta Band
        if ~isempty(runInfo.fluorChInd)
            
            fStr = [num2str(fRange_delta(1)) '-' num2str(fRange_delta(2))];
            fStr(strfind(fStr,'.')) = 'p';
            ContFig2Save= {'Calcium'};
            if ~isempty(runInfo.FADChInd)
                ContFig2Save= {'Calcium','FAD'};
            end
            clear fh
            if ~strcmp(runInfo.system,'EastOIS1')
                [fh, seedFC, seedFCMap,seedCenter,seedRadius,bilatFCMap]...
                    = qcFC(data_delta,xform_isbrain,runInfo.Contrasts,runInfo,fRange_delta,ContFig2Save);
                
                for contrastInd = 1:numel(ContFig2Save)
                    saveFCQCFigName = [runInfo.saveFCQCFig '-' ContFig2Save{contrastInd} '-' fStr ];
                    saveas(fh(contrastInd),[saveFCQCFigName '.png']);
                    close(fh(contrastInd));
                end
            end
            %Saving
            saveSeedFCName = [runInfo.saveFilePrefix '-seedFC-' fStr];
            save(saveSeedFCName,'contrastName','seedCenter','seedFC','seedFCMap','seedRadius','-v7.3');
            saveBilateralFCName = [runInfo.saveFilePrefix '-bilateralFC-' fStr];
            save(saveBilateralFCName,'contrastName','bilatFCMap','xform_isbrain','-v7.3');
        end
        
        
        
        
    elseif runInfo.session =="stim"
        load(runInfo.saveMaskFile,'xform_isbrain')
        %Building Data Full
        if ~isempty(runInfo.fluorChInd)
            data_full = nan(size(xform_datahb,1),size(xform_datahb,2),size(xform_datahb,4),4);
        elseif ~isempty(runInfo.FADChInd)
            data_full = nan(size(xform_datahb,1),size(xform_datahb,2),size(xform_datahb,4),5);
        else
            data_full = nan(size(xform_datahb,1),size(xform_datahb,2),size(xform_datahb,4),3);
        end
        %Shape Data
        data_full(:,:,:,1)=squeeze(xform_datahb(:,:,1,:))*10^6;
        data_full(:,:,:,2)=squeeze(xform_datahb(:,:,2,:))*10^6;
        data_full(:,:,:,3)=squeeze(xform_datahb(:,:,1,:))*10^6+squeeze(xform_datahb(:,:,2,:))*10^6;
        if ~isempty(runInfo.fluorChInd)
            data_full(:,:,:,4) = xform_datafluorCorr*100;
        end
        if ~isempty(runInfo.FADChInd)
            data_full(:,:,:,5) = xform_dataFADCorr*100;
        end
        
        %Power Break- Can we initalize it?
        for i=1:size(data_full,4)
            [whole_spectra_map(:,:,:,i),powerMap(:,:,:,i),hz, global_sig_for(:,i),glob_sig_power(:,i)]= PowerAnalysis(squeeze(data_full(:,:,:,i)),runInfo.samplingRate,xform_isbrain); %anmol-semicolon
        end
        save(strcat(runInfo.saveFilePrefix,'-Power.mat'),'hz','whole_spectra_map','powerMap','global_sig_for','glob_sig_power','-v7.3')
        
        clear whole_spectra_map powerMap global_sig_for glob_sig_power
        
        %Time Ranges
        if ~isempty(runInfo.FADChInd)
            peakTimeRange = {runInfo.stimEndTime-1:runInfo.stimEndTime+1,... %HbO
                runInfo.stimEndTime-1:runInfo.stimEndTime+1,... %HbR
                runInfo.stimEndTime-1:runInfo.stimEndTime+1,... $HbT
                runInfo.stimStartTime:runInfo.stimEndTime,... %Fluor
                runInfo.stimStartTime:runInfo.stimEndTime};% FAD
        elseif ~isempty(runInfo.fluorChInd)
            peakTimeRange = {runInfo.stimEndTime-1:runInfo.stimEndTime+1,...
                runInfo.stimEndTime-1:runInfo.stimEndTime+1,...
                runInfo.stimEndTime-1:runInfo.stimEndTime+1,...
                runInfo.stimStartTime:runInfo.stimEndTime};
        else
            peakTimeRange = {runInfo.stimEndTime-1:runInfo.stimEndTime+1,...
                runInfo.stimEndTime-1:runInfo.stimEndTime+1,...
                runInfo.stimEndTime-1:runInfo.stimEndTime+1};
        end
        
        %Setting up stim
        stimBlockSize=round(runInfo.blockLen*runInfo.samplingRate); %stim length in frames
        %why is this here?
        R=mod(size(data_full,3),stimBlockSize); %are there dropped frames?
        %if there is dropped frames
        if R~=0
            pad=stimBlockSize-R;
            disp(['** Non integer number of blocks presented. Padded with ' , num2str(pad), ' zeros **'])
            data_full(:,:,end:end+pad,:)=0;
            runInfo.appendedZeros=pad;
        end
        numBlocks = round(size(data_full,3)/(runInfo.blockLen*runInfo.samplingRate));%what if not integer
        
        
        % GSR function takes concatonated data
        data_full = reshape(data_full,size(data_full,1),size(data_full,2),[],length(runInfo.Contrasts));
        for ii = 1:size(data_full,4)
            data_full(:,:,:,ii) = gsr(squeeze(data_full(:,:,:,ii)),xform_isbrain);
        end
        
        data_full = reshape(data_full,size(data_full,1),size(data_full,2),[],numBlocks,size(data_full,4)); %reshape to pixel-pixel-blockSize-numblock-species
        %Baseline Subtract
        for ii = 1:numBlocks
            for jj = length(runInfo.Contrasts)
                meanFrame = squeeze(mean(data_full(:,:,1:runInfo.stimStartTime*runInfo.samplingRate,ii,jj),3));
                data_full_BaselineShift(:,:,:,ii,jj) = data_full(:,:,:,ii,jj) - repmat(meanFrame,1,1,size(data_full,3),1,1);
                data_full_BaselineShift(:,:,:,ii,jj) = data_full_BaselineShift(:,:,:,ii,jj).*xform_isbrain;
            end
        end
        clear xform_isbrain_matrix meanFrame
        %Reshape
        data_full_BaselineShift = reshape(data_full_BaselineShift,size(data_full_BaselineShift,1),size(data_full_BaselineShift,2),[],numBlocks,length(runInfo.Contrasts));
        %Create pres maps
        
        fh= generateBlockMap(data_full_BaselineShift,runInfo.Contrasts,runInfo,numBlocks,peakTimeRange,xform_isbrain);
        %save
        sgtitle([runInfo.saveFilePrefix(17:end),'GSR']) %anmol - suptitle to sgtitle
        saveName = strcat(runInfo.saveFilePrefix,'_GSR_BlockPeak');
        saveas(gcf,strcat(saveName,'.fig'))
        saveas(gcf,strcat(saveName,'.png'))
        close all
        
        
    end
    disp(['Done Processing ', runInfo.saveHbFile ' on ', runInfo.system])
    
end


