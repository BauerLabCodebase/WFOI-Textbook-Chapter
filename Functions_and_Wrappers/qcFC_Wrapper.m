clear all;close all;
excelFile="C:\Users\Nischal\Desktop\code mod\newpipeline\Code Modification\exampleTiffOIS+Gcamp.xlsx";
excelRows=2;
fRange_ISA = [0.01 0.08];
fRange_delta = [0.5 4];
runsInfo = parseRuns(excelFile,excelRows);
[row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!

runNum = numel(runsInfo);

for runInd = 1:runNum
    runInfo=runsInfo(runInd);
    saveBilateralFCName = [runInfo.saveFilePrefix '-bilateralFC-' fStr];
    saveSeedFCName = [runInfo.saveFilePrefix '-seedFC-' fStr];
    if isfile(saveBilateralFCName)
    disp(strcat(runInfo.mouseName, '-',runInfo.run, ' already processed.'))
       continue
    else
    
    
%Loading    
try
        load(runInfo.saveHbFile,'xform_datahb')
        load(runInfo.saveMaskFile,'xform_isbrain')
        if ~isempty(runInfo.fluorChInd)
            load(runInfo.saveFluorFile,'xform_datafluorCorr')
        end
        if ~isempty(runInfo.FADChInd)
            load(runInfo.saveFADFile,'xform_dataFADCorr')
        end    
catch
        warning(strcat(runInfo.mouseName, '-',runInfo.run, ' Not processed.'))
end

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

        clear  xform_datafluorCorr xform_datahb xform_dataFADCorr
            
        %Filter
        data_isa= nan(size(data_full));
        data_delta= nan(size(data_full));
        for ii=1:size(data_full,4) 
            data_isa(:,:,:,ii)=filterData(data_full(:,:,:,ii),fRange_ISA(1),fRange_ISA(2),runInfo.samplingRate,xform_isbrain);
            data_delta(:,:,:,ii)=filterData(data_full(:,:,:,ii),fRange_delta(1),fRange_delta(2),runInfo.samplingRate,xform_isbrain);
        end
                
        clear data_full
        
        %GSR
        for  ii=1:size(data_isa,4)
            data_isa(:,:,:,ii) = gsr(squeeze(data_isa(:,:,:,ii)),xform_isbrain);
            data_delta(:,:,:,ii) = gsr(squeeze(data_delta(:,:,:,ii)),xform_isbrain);
        end
        
            xform_isbrain_matrix = repmat(xform_isbrain,1,1,size(data_delta,3),size(data_delta,4));
            data_isa(logical(1-xform_isbrain_matrix)) = NaN;%do we
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
            save(saveSeedFCName,'contrastName','seedCenter','seedFC','seedFCMap','seedRadius','-v7.3');
            save(saveBilateralFCName,'contrastName','bilatFCMap','xform_isbrain','-v7.3');
        end
    end
        %%
end