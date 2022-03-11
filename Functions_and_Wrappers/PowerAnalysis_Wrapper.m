function PowerAnalysis_Wrapper(excelFile,excelRows) % anmol- turned powerAnalysis_Wrapper into a function

runsInfo = parseRuns(excelFile,excelRows);
[row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!

runNum = numel(runsInfo);

for runInd = 1:runNum
    runInfo=runsInfo(runInd);
    if isfile([runInfo.saveFilePrefix,'-Power.mat'])
        disp(strcat('Power Analysis for ',runInfo.mouseName,' run ' ,num2str(runInfo.run),' already exists!'))
        continue
    end % anmol - don't need an else statement here
    disp(strcat('Running Power Analysis for ',runInfo.mouseName,' run ' ,num2str(runInfo.run))) %anmol- changed the order of display, so I wouldn't get redundant messages
    
    %Loading

            try
                load(runInfo.saveHbFile,'xform_datahb')
                load(runInfo.saveMaskFile,'xform_isbrain')
            catch
                warning(strcat(runInfo.mouseName, '-',runInfo.run, ' Not processed.'))
                continue
            end

            if ~isempty(runInfo.fluorChInd)
                load(runInfo.saveFluorFile,'xform_datafluorCorr');
            end
            if ~isempty(runInfo.FADChInd)
                load(runInfo.saveFADFile,'xform_dataFADCorr');
            end

    %Initializing    
            if ~isempty(runInfo.fluorChInd) && isempty(runInfo.FADChInd) % Xiaodan add && isempty(runInfo.FADChInd)
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
%Power Break
            for i=1:size(data_full,4)
        [whole_spectra_map(:,:,:,i),avg_cort_spec(:,i),powerMap(:,:,:,i),hz, global_sig_for(:,i),glob_sig_power(:,i)]= PowerAnalysis(squeeze(data_full(:,:,:,i)),runInfo.samplingRate,xform_isbrain); 
            end

            save(strcat(runInfo.saveFilePrefix,'-Power.mat'),'hz','whole_spectra_map','powerMap','global_sig_for','glob_sig_power','avg_cort_spec','-v7.3')
          end
    
    save(strcat(runInfo.saveFilePrefix,'-Power.mat'),'hz','whole_spectra_map','powerMap','global_sig_for','glob_sig_power','-v7.3')
end

end