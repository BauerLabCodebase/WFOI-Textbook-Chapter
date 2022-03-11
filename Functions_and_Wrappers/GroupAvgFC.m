function GroupAvgFC(excelFile,excelRows,Group_directory,varargin) %anmol - vararging deleted for some reason?
%note that group directory would be something like 'A:\JPC\Old
%Lis\Rad_project

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
    
%loading fc params
paramPath = what('bauerParams');
seedsData = load(fullfile(paramPath(1).path,'seeds16.mat')); 
seedNames = seedsData.seedNames;
seedNum = size(seedsData.seedCenter,1);

%getting excel/mouse info
runsInfo = parseRuns(excelFile,excelRows);
[excel_row,first_ind_mouse,~]=unique({runsInfo.excelRow_char});
excel_row=str2double(excel_row); 
excelData = readtable(excelFile);
tableRows = excel_row - 1;

%Get Group labels and ignore, non-labled mice
groupLabel=cell(numel(excel_row),1);
ct=1; 
for row = tableRows
    
     if isnumeric(excelData{row,'Timepoint'}) & ~isnan(excelData{row,'Timepoint'})
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

%Create group directory
if ~exist(Group_directory)
    mkdir(Group_directory)
end


%Group Averages: FOR EACH GROUP
    for group_indx=1:numel(group_names) %this is the number of group labels
        %loading all mice
        if ~contains(group_names{group_indx},'unassigned') %only load mice with labels
            xform_isbrain_intersect=ones(128,128); %get the intersection of all the isbrains
            tmp=[find(group_id==group_indx)]; %find which mice have are in the current group
            
            for i= 1:numel(tmp) %Load all of these mice
                runInfo=runsInfo( first_ind_mouse(tmp(i)) );
                disp(['Loaded ', runInfo.mouseName, ' runs: ' num2str(i) ', group:' num2str(group_indx)])
                
                load(runInfo.saveMaskFile)
                xform_isbrain_intersect=xform_isbrain_intersect.*xform_isbrain;
                
                load_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgFC_', range_name);
                load(load_name, ...
                    'seedFC_mouse_avg','seedFCMap_mouse_avg','bilatFCMap_mouse_avg','seedNames','seedCenter','seedRadius') %anmol-load in seedCenter
                
                tmp_seedFC_mouse(:,:,:,i)=seedFC_mouse_avg;
                tmp_seedFCMAP_mouse(:,:,:,:,i)=seedFCMap_mouse_avg;
                tmp_bilatFCMap_mouse(:,:,:,i)=bilatFCMap_mouse_avg;
                
            end
            
            %FISHER SPACE AVERAGING
            seedFC_group_avg=squeeze(tanh(mean(atanh(tmp_seedFC_mouse),4,'omitnan'))); %anmol- antanh to atanh
            seedFCMap_group_avg=squeeze(tanh(mean(atanh(tmp_seedFCMAP_mouse),5,'omitnan')));
            bilatFCMap_group_avg=squeeze(tanh(mean(atanh(tmp_bilatFCMap_mouse),4,'omitnan')));
            
            
            %PLOTTING
            for contrast = 1:size(seedFCMap_group_avg,4)
                % get the data
                speciesMat = squeeze(seedFC_group_avg(:,:,contrast));
                % get the mask
                cMask = xform_isbrain_intersect;
                seedMap=MakeSeedMap(128,128,seedCenter,seedRadius);
                %plot
                fh(contrast) = figure('Position',[100 100 1200 500]);
                %seedMap and Matrix
                seedFCPos = [0.02 0.05 0.6 0.85];
                fh(contrast) = plotSeedFC(fh(contrast),seedFCPos,...
                    squeeze(seedFCMap_group_avg(:,:,:,contrast)),squeeze(seedFC_group_avg(:,:,contrast)),...
                    seedCenter, seedMap, seedRadius,seedNames,cMask);
                %Bilat
                bilatFCPos = [0.62 0.05 0.36 0.85];
                MaskMirrored = cMask & fliplr(cMask);
                MaskMirrored(:,[63:66])=0; %anmol- masked the midline
                fh(contrast) = plotFCMap(fh(contrast),bilatFCPos,squeeze(bilatFCMap_group_avg(:,:,contrast)),MaskMirrored);
                % add the title
                titleAxesHandle=axes('position',[0 0 1 0.95]);
                t = title(titleAxesHandle,[group_names{group_indx} ' ' runInfo.Contrasts{contrast} ' FC ' fStr 'Hz']);
                set(titleAxesHandle,'visible','off');
                set(t,'visible','on');
                
                
                
                % save fig
                saveFigName=strcat(Group_directory,filesep, group_names{group_indx},'-', runInfo.Contrasts{contrast},'-Average',range_name)
                saveas(gcf,strcat(saveFigName, '.png'));
                saveas(gcf,strcat(saveFigName, '.fig'));
                close(gcf);
                
                
                %OnlyBilat
                fh1=figure;
                plotFCMap(fh1,[],squeeze(bilatFCMap_group_avg(:,:,contrast)),MaskMirrored);
                titleAxesHandle=axes('position',[0 0 1 0.95]);
                t = title([group_names{group_indx} ' ' runInfo.Contrasts{contrast} ' FC ' fStr 'Hz']);
                set(titleAxesHandle,'visible','off');
                set(t,'visible','on');
                
                saveFigName1=strcat(Group_directory,filesep, group_names{group_indx},'-', runInfo.Contrasts{contrast},'-Average-BilatOnly',range_name)
                saveas(gcf,strcat(saveFigName1, '.png'));
                saveas(gcf,strcat(saveFigName1, '.fig'));
                close(gcf);
            end
%Saving
        save_name=strcat(Group_directory,filesep,group_names{group_indx}, '-avgFC_',range_name)
        save(save_name, ...
                'seedFC_group_avg','seedFCMap_group_avg','bilatFCMap_group_avg','seedNames','seedCenter','seedRadius','xform_isbrain_intersect','-v7.3') 
        end
        
        clearvars -except paramPath seedsData seedNames seedNum excel_row first_ind_mouse mouse_ids range_name runsInfo fStr group_names group_id Group_directory 

    end

end