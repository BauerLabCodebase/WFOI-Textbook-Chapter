function MouseAvgStim(excelFile,excelRows)

runsInfo = parseRuns(excelFile,excelRows);
[excel_row,first_ind_mouse,mouse_ids]=unique({runsInfo.excelRow_char});

%Mouse Averages
for mouse_indx=1:length(first_ind_mouse) %this is the number of mice
    inds=find(mouse_ids==mouse_indx);  %how many runs for each mice and which one is the associated runsInfo?  
    runInfo=runsInfo(first_ind_mouse(mouse_indx));      
    save_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgStim'); 

    
    load(runInfo.saveMaskFile,'xform_isbrain'); 

    for mouse_run=1:length(inds) %load for each run
            ind=inds(mouse_run);   %associated index in runsInfo
            load([runsInfo(ind).saveFilePrefix,'-Power','.mat']); %load

            tmp_whole_spectra_map(:,:,:,:,mouse_run)=whole_spectra_map.*noVascMask;
            tmp_powerMap(:,:,:,:,mouse_run)=powerMap.*noVascMask;
            tmp_global_sig_for(:,:,mouse_run)=global_sig_for;%anmol-deleted *.Novascmask
            tmp_glob_sig_power(:,:,mouse_run)=glob_sig_power;%.*noVascMask;
            
    end 

 load(runInfo.saveMaskFile,'xform_isbrain')
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