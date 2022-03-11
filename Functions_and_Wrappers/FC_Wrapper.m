close all;clear all;

excelFile = "L:\WhyDeOxyLowSNR\WhyDeOxyLowSNR.xlsx";
excelRows = 10;
fRange_ISA = [0.01 0.08];
fRange_delta = [0.5 4];
runsInfo = parseRuns(excelFile,excelRows);
[row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!

runNum = numel(runsInfo);


%% ISA
for runInd = 1:runNum
    clearvars -except runsInfo runNum runInd fRange_ISA fRange_delta
    runInfo=runsInfo(runInd);

%loading   
    disp('Loaded File')
    load(runInfo.saveHbFile,'xform_datahb','sessionInfo','op','E','runInfo','-v7.3')    
    
    data_full(:,:,1:size(xform_datahb,4),1)=squeeze(xform_datahb(:,:,1,:));
    data_full(:,:,1:size(xform_datahb,4),2)=squeeze(xform_datahb(:,:,2,:));
    data_full(:,:,1:size(xform_datahb,4),3)=squeeze(xform_datahb(:,:,1,:)+xform_datahb(:,:,2,:));    
    
    if ~isempty(runInfo.FADChInd)
    load(runInfo.saveFADFile,'xform_dataFAD','xform_dataFADCorr','op_inFAD', 'E_inFAD', 'op_outFAD', 'E_outFAD','runInfo','-v7.3')
    data_full(:,:,1:size(xform_FADCorr,3),4)=xform_FADCorr;
    data_full(:,:,1:size(xform_jrgeco1aCorr,3),5)=xform_jrgeco1aCorr;        
    end
    
    if ~isempty(runInfo.fluorChInd)
        load(runInfo.saveFluorFile,'xform_datafluor','xform_datafluorCorr','op_in', 'E_in', 'op_out', 'E_out','runInfo','-v7.3')
    end
    data_full(isnan(data_full)) = 0;
    data_full(isinf(data_full)) = 0;
    

%start FC Processing 
    for ii=1:size(data_full,4)
        data_delta(:,:,:,ii)=filterData(data_full(:,:,:,ii),0.4,4,20);
        data_delta(:,:,:,ii) = gsr(squeeze(data_delta(:,:,:,ii)),xform_isbrain);
        data_ISA(:,:,:,ii)=filterData(data_full(:,:,:,ii),0.009,0.08,20);
        data_ISA(:,:,:,ii) = gsr(squeeze(data_ISA(:,:,:,ii)),xform_isbrain);
    end
    data_full(data_full==0) = nan;
    clear data_full xform_datahb xform_FADCorr xform_jrgeco1aCorr
    
    %ISA FC
    data=data_ISA;
    clear data_ISA
    fRange =fRange_ISA;
    fStr = [num2str(fRange(1)) '-' num2str(fRange(2))];
    
    refseeds=GetReferenceSeeds_xw;
    seedCenter = flip(refseeds,2);
    
    seedNames = {'Fr-L'  'Cing-L' 'M-L' 'SS-L' 'RS-L' 'P-L' 'V-L' 'Aud-L'...
        'Fr-R' 'Cing-R' 'M-R' 'SS-R' 'RS-R' 'P-R' 'V-R' 'Aud-R'};
    seedNum = size(seedCenter,1);
    seedFC = nan(seedNum,seedNum,size(data,4));
    seedFCMap = nan(size(data,1), size(data,2),seedNum,size(data,4) );
    bilatFCMap =nan(size(data,1), size(data,2),size(data,4) );
    % run for each contrast
    
    for contrast = 1:size(data,4)
        cData = squeeze(data(:,:,:,contrast));
        
        %Get seeds
        numPixY = size(cData,1); % number of pixels in y axis
        numPixX = size(cData,2);
        sizeY = runInfo.window(2)-runInfo.window(1); % size in mm in y axis
        seedRadius = 0.25; % in mm
        seedRadius = round(seedRadius/sizeY*numPixY);
        
        %Get seed maps
        seedMap = false(numPixY,numPixX,seedNum);
        for seedInd = 1:seedNum
            seedCoor = circleCoor(seedCenter(seedInd,:),seedRadius);
            seedCoor = matCoor2Ind(seedCoor,[numPixY numPixX]);
            seedBool = false(numPixY,numPixX); seedBool(seedCoor) = true;
            seedMap(:,:,seedInd) = seedBool;
        end
        
        % get seed fc
        for seedInd = 1:seedNum
            seedFCMap(:,:,seedInd,contrast) = seedFCMap_fun(cData,seedMap(:,:,seedInd));
        end
        %fc matrix
        seedFC(:,:,contrast) = seedFC_fun(cData,seedMap);
        % get bilateral fc
        bilatFCMap(:,:,contrast) = bilateralFC_fun(cData);
        
    end
%Saving
        saveSeedFCName = [runInfo.saveFilePrefix '-seedFC-' fStr];
        save(saveSeedFCName,'contrastName','seedCenter','seedFC','seedFCMap','seedRadius','-v7.3');
        saveBilateralFCName = [runInfo.saveFilePrefix '-bilateralFC-' fStr];
        save(saveBilateralFCName,'contrastName','bilatFCMap','xform_isbrain','-v7.3');
            

    clear data
    
    
    
    
    
    %Delta FC
    data=data_delta;
    clear data_delta
    fRange =fRange_delta;
    fStr = [num2str(fRange(1)) '-' num2str(fRange(2))];
    
    refseeds=GetReferenceSeeds_xw;
    seedCenter = flip(refseeds,2);
    
    seedNames = {'Fr-L'  'Cing-L' 'M-L' 'SS-L' 'RS-L' 'P-L' 'V-L' 'Aud-L'...
        'Fr-R' 'Cing-R' 'M-R' 'SS-R' 'RS-R' 'P-R' 'V-R' 'Aud-R'};
    seedNum = size(seedCenter,1);
    seedFC = nan(seedNum,seedNum,size(data,4));
    seedFCMap = nan(size(data,1), size(data,2),seedNum,size(data,4) );
    bilatFCMap =nan(size(data,1), size(data,2),size(data,4) );
    % run for each contrast
    
    for contrast = 1:size(data,4)
        cData = squeeze(data(:,:,:,contrast));
        
        %Get seeds
        numPixY = size(cData,1); % number of pixels in y axis
        numPixX = size(cData,2);
        sizeY = runInfo.window(2)-runInfo.window(1); % size in mm in y axis
        seedRadius = 0.25; % in mm
        seedRadius = round(seedRadius/sizeY*numPixY);
        
        %Get seed maps
        seedMap = false(numPixY,numPixX,seedNum);
        for seedInd = 1:seedNum
            seedCoor = circleCoor(seedCenter(seedInd,:),seedRadius);
            seedCoor = matCoor2Ind(seedCoor,[numPixY numPixX]);
            seedBool = false(numPixY,numPixX); seedBool(seedCoor) = true;
            seedMap(:,:,seedInd) = seedBool;
        end
        
        % get seed fc
        for seedInd = 1:seedNum
            seedFCMap(:,:,seedInd,contrast) = seedFCMap_fun(cData,seedMap(:,:,seedInd));
        end
        %fc matrix
        seedFC(:,:,contrast) = seedFC_fun(cData,seedMap);
        % get bilateral fc
        bilatFCMap(:,:,contrast) = bilateralFC_fun(cData);
        
    end
    %Saving
        saveSeedFCName = [runInfo.saveFilePrefix '-seedFC-' fStr];
        save(saveSeedFCName,'contrastName','seedCenter','seedFC','seedFCMap','seedRadius','-v7.3');
        saveBilateralFCName = [runInfo.saveFilePrefix '-bilateralFC-' fStr];
        save(saveBilateralFCName,'contrastName','bilatFCMap','xform_isbrain','-v7.3');
            
    
end










