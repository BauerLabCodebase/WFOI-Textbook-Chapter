%% Get Run and System Info
function MouseAvgLag(excelFile,excelRows)

    runsInfo = parseRuns(excelFile,excelRows);
    [row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!
    runNum=numel(runsInfo);

    load(which('GoodWL.mat'),'xform_WL'); %anmol- changed hardcoding here  - more than 2 ?
    load(which('noVasculatureMask.mat')); %anmol- changed hardcoding here

%% Loading Individual Runs
    for mouse_indx=1:length(start_ind_mouse) 
        inds=find(numruns_per_mouse==mouse_indx);  %how many runs for each mice and which one is the associated runsInfo?  
        runInfo=runsInfo(start_ind_mouse(mouse_indx));    
        save_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgLag'); 
        
        if isfile(strcat(save_name,'.mat')) 
            disp([strcat(save_name,'.mat') ' already exists'])
            continue 
        end


        if isempty(runInfo.FADChInd)
            
            for mouse_run=1:length(inds) %load for each run
                ind=inds(mouse_run);   %associated index in runsInfo
                load([runsInfo(ind).saveFilePrefix '_lag.mat'],'lagTimeTrial_HbTCalcium', 'lagAmpTrial_HbTCalcium');
                tmp_lagTimeTrial_HbTCalcium(:,:,mouse_run)=lagTimeTrial_HbTCalcium;
                tmp_lagAmpTrial_HbTCalcium(:,:,mouse_run)=lagAmpTrial_HbTCalcium;
                
            end 
            
            mouse_lagTimeTrial_HbTCalcium=squeeze(mean(tmp_lagTimeTrial_HbTCalcium,3,'omitnan'));
            mouse_lagAmpTrial_HbTCalcium=squeeze(mean(tmp_lagAmpTrial_HbTCalcium,3,'omitnan'));
            
            
        else
            for mouse_run=1:length(inds) %load for each run
                ind=inds(mouse_run);   %associated index in runsInfo
                load([runsInfo(ind).saveFilePrefix '_lag.mat'],'lagTimeTrial_HbTCalcium', 'lagAmpTrial_HbTCalcium','lagTimeTrial_FADCalcium', 'lagAmpTrial_FADCalcium','lagTimeTrial_HbTFAD','lagAmpTrial_HbTFAD');
                tmp_lagTimeTrial_HbTCalcium(:,:,mouse_run)=lagTimeTrial_HbTCalcium;
                tmp_lagAmpTrial_HbTCalcium(:,:,mouse_run)=lagAmpTrial_HbTCalcium;
                tmp_lagTimeTrial_FADCalcium(:,:,mouse_run)=lagTimeTrial_FADCalcium;
                tmp_lagAmpTrial_FADCalcium(:,:,mouse_run)=lagAmpTrial_FADCalcium; 
                tmp_lagTimeTrial_HbTFAD(:,:,mouse_run)=lagTimeTrial_HbTFAD;
                tmp_lagAmpTrial_HbTFAD(:,:,mouse_run)=lagAmpTrial_HbTFAD;
            end
            
            mouse_lagTimeTrial_HbTCalcium=squeeze(mean(tmp_lagTimeTrial_HbTCalcium,3,'omitnan'));
            mouse_lagAmpTrial_HbTCalcium=squeeze(mean(tmp_lagAmpTrial_HbTCalcium,3,'omitnan'));
            
            mouse_lagTimeTrial_FADCalcium=squeeze(mean(tmp_lagTimeTrial_FADCalcium,3,'omitnan'));
            mouse_lagAmpTrial_FADCalcium=squeeze(mean(tmp_lagAmpTrial_FADCalcium,3,'omitnan'));
            
            mouse_lagTimeTrial_HbTFAD=squeeze(mean(tmp_lagTimeTrial_HbTFAD,3,'omitnan'));
            mouse_lagAmpTrial_HbTFAD=squeeze(mean(tmp_lagAmpTrial_HbTFAD,3,'omitnan'));
            
        end
   
    %WHY IS THIS HERE AND HOW CAN WE NOT HAVE USER INPUT HALF WAY INTO OUR CODE?
        edgeLen =1; % this means how many points to the left and right of the smapled maximum
        tZone = 4; %how many seconds to look at the lag (only minus) when use xcorr
        corrThr = 0; %if the maximum correlation is less than 0 at that pixel is less than this then eliminate (only look at positive correlation)
        validRange = - edgeLen: round(tZone*runInfo.samplingRate); %     How many data points to consider from zero lag. For example, -3:3 would mean looking at 7 points around zero lag.
        tLim = [-2 2]; %for visualizations for calcium and hb (calcium comes first)
        tLim_FAD = [0 0.3]; % %for visualizations
        rLim = [-1 1]; %for visualizations


        if isempty(runInfo.FADChInd)
            figure;
            colormap jet;
            subplot(2,1,1); imagesc(mouse_lagTimeTrial_HbTCalcium,tLim); axis image off;h = colorbar;ylabel(h,'t(s)');title('Calcium HbT');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,1,2); imagesc(mouse_lagAmpTrial_HbTCalcium,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            sgtitle(strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session))

            saveas(gcf,[save_name '.png']);% anmol - savename changes
            saveas(gcf,[save_name '.fig']);
            save([save_name,'.mat']...
                ,'mouse_lagTimeTrial_HbTCalcium', 'mouse_lagAmpTrial_HbTCalcium', '-v7.3');
            close all
        else %if we have FAD

            figure;
            colormap jet;
            subplot(2,3,1); imagesc(mouse_lagTimeTrial_HbTCalcium,tLim); axis image off;h = colorbar;ylabel(h,'t(s)');title('Calcium HbT');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,2); imagesc(mouse_lagTimeTrial_FADCalcium,tLim_FAD);axis image off;h = colorbar;ylabel(h,'t(s)');title('FAD Calcium');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,3); imagesc(mouse_lagTimeTrial_HbTFAD,[-1 1]);axis image off;h = colorbar;ylabel(h,'t(s)');title('FAD HbT');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);

            subplot(2,3,4); imagesc(mouse_lagAmpTrial_HbTCalcium,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,5); imagesc(mouse_lagAmpTrial_FADCalcium,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,6); imagesc(mouse_lagAmpTrial_HbTFAD,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            sgtitle(strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session))

            saveas(gcf,fullfile(save_name,'.png'));% anmol - savename  changes
            saveas(gcf,fullfile(save_name,'.fig'));

            save(fullfile(save_name,'.mat')...
                ,'mouse_lagTimeTrial_HbTCalcium', 'mouse_lagAmpTrial_HbTCalcium','mouse_lagTimeTrial_FADCalcium', ...
                'mouse_lagAmpTrial_FADCalcium','mouse_lagTimeTrial_HbTFAD','mouse_lagAmpTrial_HbTFAD','-v7.3');
            close all

        end                                
    end
end

