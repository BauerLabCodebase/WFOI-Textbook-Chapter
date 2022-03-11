%% Get Run and System Info
function GroupAvgLag(excelFile,excelRows,Group_directory)

    runsInfo = parseRuns(excelFile,excelRows);
    [row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); %Note that unique only takes characters! This makes it so that we only do landmark for one of the runs!
    runNum=numel(runsInfo);

    excel_row=str2double(row);
    excelData = readtable(excelFile);
    tableRows = excel_row - 1;


    load(which('GoodWL.mat'),'xform_WL'); 
    load(which('noVasculatureMask.mat'));
    
    
 %Get Group labels and ignore, non-labled mice
    groupLabel=cell(numel(excel_row),1);
    ct=1; 
 for row = tableRows

     if isnumeric(excelData{row,'Timepoint'}) & ~isnan(excelData{row,'Timepoint'}) %anmol- only numbr thing
        tp=num2str(excelData{row,'Timepoint'});
     elseif isnan(excelData{row,'Timepoint'})
         tp='-';
     else
        tp=excelData{row,'Timepoint'};
    end
    
    groupLabel(ct) = strcat(excelData{row,'Group'},'-',tp);
        if isempty(groupLabel{ct})
            groupLabel(ct)={'unassigned'};
            warning(['row ' num2str(row) ' has no group assignment. It will not be included.'])
        end
        ct=ct+1; 
    end

    [group_names,~,group_id]=unique(groupLabel);

    if ~exist(Group_directory)
       mkdir(Group_directory)
    end

%iterate over group    
    for group_indx=1:numel(group_names) %this is the number of group labels
        
        save_name=strcat(Group_directory,filesep,group_names{group_indx}, '-avgLag')

        tmp=[find(group_id==group_indx)]; %mice in the group
        runInfo=runsInfo( start_ind_mouse(1)); %load for the if statement -  assumes every type of mouse in a group is the same
        
        if isempty(runInfo.FADChInd)
            for i=1:numel(tmp) %load for each run
                runInfo=runsInfo(start_ind_mouse(tmp(i)));
                load_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgLag');

                load([load_name '.mat'],'mouse_lagTimeTrial_HbTCalcium', 'mouse_lagAmpTrial_HbTCalcium');
                tmp_lagTimeTrial_HbTCalcium(:,:,i)=mouse_lagTimeTrial_HbTCalcium;
                tmp_lagAmpTrial_HbTCalcium(:,:,i)=mouse_lagAmpTrial_HbTCalcium;

            end

            group_lagTimeTrial_HbTCalcium=squeeze(mean(tmp_lagTimeTrial_HbTCalcium,3,'omitnan'));
            group_lagAmpTrial_HbTCalcium=squeeze(mean(tmp_lagAmpTrial_HbTCalcium,3,'omitnan'));

        else
            for i=1:numel(tmp) %load for each run
                
                runInfo=runsInfo( start_ind_mouse(tmp(i)) );
                load_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgLag');
                
                load([load_name '.mat'],'mouse_lagTimeTrial_HbTCalcium', 'mouse_lagAmpTrial_HbTCalcium','mouse_lagTimeTrial_FADCalcium', 'mouse_lagAmpTrial_FADCalcium','mouse_lagTimeTrial_HbTFAD','mouse_lagAmpTrial_HbTFAD');
                tmp_lagTimeTrial_HbTCalcium(:,:,i)=mouse_lagTimeTrial_HbTCalcium;
                tmp_lagAmpTrial_HbTCalcium(:,:,i)=mouse_lagAmpTrial_HbTCalcium;
                tmp_lagTimeTrial_FADCalcium(:,:,i)=mouse_lagTimeTrial_FADCalcium;
                tmp_lagAmpTrial_FADCalcium(:,:,i)=mouse_lagAmpTrial_FADCalcium; 
                tmp_lagTimeTrial_HbTFAD(:,:,i)=mouse_lagTimeTrial_HbTFAD;
                tmp_lagAmpTrial_HbTFAD(:,:,i)=mouse_lagAmpTrial_HbTFAD;
            end

            group_lagTimeTrial_HbTCalcium=squeeze(mean(tmp_lagTimeTrial_HbTCalcium,3,'omitnan'));
            group_lagAmpTrial_HbTCalcium=squeeze(mean(tmp_lagAmpTrial_HbTCalcium,3,'omitnan'));

            group_lagTimeTrial_FADCalcium=squeeze(mean(tmp_lagTimeTrial_FADCalcium,3,'omitnan'));
            group_lagAmpTrial_FADCalcium=squeeze(mean(tmp_lagAmpTrial_FADCalcium,3,'omitnan'));

            group_lagTimeTrial_HbTFAD=squeeze(mean(tmp_lagTimeTrial_HbTFAD,3,'omitnan'));
            group_lagAmpTrial_HbTFAD=squeeze(mean(tmp_lagAmpTrial_HbTFAD,3,'omitnan'));
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
            subplot(2,1,1); imagesc(group_lagTimeTrial_HbTCalcium,tLim); axis image off;h = colorbar;ylabel(h,'t(s)');title('Lag: Calcium and HbT');hold on;imagesc(xform_WL,'AlphaData',1-mask_new); %changed title
            subplot(2,1,2); imagesc(group_lagAmpTrial_HbTCalcium,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            sgtitle(strcat(group_names{group_indx} , '-Lag'));


            saveas(gcf,strcat(save_name,'.png'));% anmol - savename changes
            saveas(gcf,strcat(save_name ,'.fig'));
            save(strcat(save_name,'.mat')...
                ,'group_lagTimeTrial_HbTCalcium', 'group_lagAmpTrial_HbTCalcium', '-v7.3');
            close all
        else %if we have FAD

            figure;
            colormap jet;
            subplot(2,3,1); imagesc(group_lagTimeTrial_HbTCalcium,tLim); axis image off;h = colorbar;ylabel(h,'t(s)');title('Lag: Calcium and HbT');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,2); imagesc(group_lagTimeTrial_FADCalcium,tLim_FAD);axis image off;h = colorbar;ylabel(h,'t(s)');title('FAD Calcium');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,3); imagesc(group_lagTimeTrial_HbTFAD,[-1 1]);axis image off;h = colorbar;ylabel(h,'t(s)');title('FAD HbT');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);

            subplot(2,3,4); imagesc(group_lagAmpTrial_HbTCalcium,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,5); imagesc(group_lagAmpTrial_FADCalcium,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            subplot(2,3,6); imagesc(group_lagAmpTrial_HbTFAD,rLim);axis image off;h = colorbar;ylabel(h,'r');hold on;imagesc(xform_WL,'AlphaData',1-mask_new);
            sgtitle(strcat(group_names{group_indx} , '-Lag'));

            saveas(gcf,strcat(save_name ,'.png'));% anmol - savename changes
            saveas(gcf,strcat(save_name ,'.fig'));

            save(strcat(save_name,'.mat')...
                ,'group_lagTimeTrial_HbTCalcium', 'group_lagAmpTrial_HbTCalcium','group_lagTimeTrial_FADCalcium', ...
                'group_lagAmpTrial_FADCalcium','group_lagTimeTrial_HbTFAD','group_lagAmpTrial_HbTFAD','-v7.3');
            close all

        end                                
    end
end

