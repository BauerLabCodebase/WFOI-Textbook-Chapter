%%Get list of mice and runs, and load into a cancat file
function MouseAvgFC(excelFile,excelRows,varargin)
if numel(varargin) > 0
    fRange = varargin{1};
else
    fRange = [0.01 0.08];
end

if (fRange(1)<=0.01 && fRange(2) >= 0.08)
range_name='ISA';
elseif (fRange(1)<=0.5 && fRange(2) >=4 )
range_name='Delta';
else
    disp('Please Double Check Frequency Range')
end
fStr = [num2str(fRange(1)) '-' num2str(fRange(2))];
fStr(strfind(fStr,'.')) = 'p';
runsInfo = parseRuns(excelFile,excelRows);

[excel_row,first_ind_mouse,mouse_ids]=unique({runsInfo.excelRow_char});

%loading fc params
paramPath = what('bauerParams');
seedsData = load(fullfile(paramPath.path,'seeds16.mat')); %make sure this is the right seeds!
seedNames = seedsData.seedNames;
seedNum = size(seedsData.seedCenter,1);


%Mouse Averages: FOR EACH MOUSE
    for mouse_indx=1:length(first_ind_mouse) %this is the number of mice
        inds=find(mouse_ids==mouse_indx);  %how many runs -- and which ones -- for each mice and which one is the associated runsInfo?
        runInfo=runsInfo(first_ind_mouse(mouse_indx));
        
        save_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgFC_', range_name);
        if isfile(strcat(save_name,'.mat')) 
            disp([strcat(save_name,'.mat') ' already exists'])
            continue
        end

        for mouse_run=1:length(inds) %load for each run
            ind=inds(mouse_run);   %associated index in runsInfo
            load([runsInfo(ind).saveFilePrefix,'-seedFC-',fStr,'.mat']); % loads seedcenter!
            load([runsInfo(ind).saveFilePrefix,'-bilateralFC-',fStr,'.mat']); %some sort of load
            load(runsInfo(ind).saveMaskFile,'xform_WL','xform_isbrain');

            disp(['Loaded ', runInfo.mouseName, ' runs: ' num2str(runsInfo(ind).run)])

            tmp_seedFC_mouse(:,:,:,mouse_run)=seedFC; 
            tmp_seedFCMAP_mouse(:,:,:,:,mouse_run)=seedFCMap.*xform_isbrain;
            tmp_bilatFCMap_mouse(:,:,:,mouse_run)=bilatFCMap.*xform_isbrain;
        end

        %FISHER SPACE
        seedFC_mouse_avg=squeeze(tanh(mean(atanh(tmp_seedFC_mouse),4,'omitnan')));
        seedFCMap_mouse_avg=squeeze(tanh(mean(atanh(tmp_seedFCMAP_mouse),5,'omitnan')));
        bilatFCMap_mouse_avg=squeeze(tanh(mean(atanh(tmp_bilatFCMap_mouse),4,'omitnan')));


        %Plotting
        for contrast = 1:size(seedFCMap_mouse_avg,4)
            % get the data
            speciesMat = squeeze(seedFC_mouse_avg(:,:,contrast));
            % get the mask
            cMask = xform_isbrain;
            seedMap=MakeSeedMap(128,128,seedCenter,seedRadius);
            %plot
            fh(contrast) = figure('Position',[100 100 1200 500]);
            %seedMap and Matrix
            seedFCPos = [0.02 0.05 0.6 0.85];
            fh(contrast) = plotSeedFC(fh(contrast),seedFCPos,...
                squeeze(seedFCMap_mouse_avg(:,:,:,contrast)),squeeze(seedFC_mouse_avg(:,:,contrast)),...
                seedCenter, seedMap, seedRadius,seedNames,cMask);
            %Bilat
            bilatFCPos = [0.62 0.05 0.36 0.85];
            MaskMirrored = cMask & fliplr(cMask);
            MaskMirrored(:,[63:66])=0; %anmol- masked the midline

            fh(contrast) = plotFCMap(fh(contrast),bilatFCPos,...
                squeeze(bilatFCMap_mouse_avg(:,:,contrast)),MaskMirrored);
            % add the title
            titleAxesHandle=axes('position',[0 0 1 0.95]);
            t = title(titleAxesHandle,[runInfo.mouseName ' ' runInfo.Contrasts{contrast} ' FC ' fStr 'Hz']);
            set(titleAxesHandle,'visible','off');
            set(t,'visible','on');
            % save fig
            saveFigName = strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-',runInfo.mouseName,'-', runInfo.Contrasts{contrast},'Average',range_name);
            saveas(gcf,strcat(saveFigName, '.png'));
            saveas(gcf,strcat(saveFigName, '.fig'));
            close(gcf);
        end

    %Saving each band
    save(save_name, ...
        'seedFC_mouse_avg','seedFCMap_mouse_avg','bilatFCMap_mouse_avg','seedNames','seedCenter','seedRadius','-v7.3')    

        clearvars -except paramPath seedsData seedNames seedNum excel_row first_ind_mouse mouse_ids range_name runsInfo fStr 

    end
    

end

