function GroupAvgPower(excelFile,excelRows,Group_directory)
runsInfo = parseRuns(excelFile,excelRows);
[excel_row,first_ind_mouse,~]=unique({runsInfo.excelRow_char});
excel_row=str2double(excel_row);
excelData = readtable(excelFile);
tableRows = excel_row - 1;

%no vasc mask
noVascMask=ones(128,128);
noVascMask(:,61:67)=0;

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

if ~exist(Group_directory)
   mkdir(Group_directory)
end

%group averages for each group!
for group_indx=1:numel(group_names) %this is the number of group labels
    %loading all mice
    if ~contains(group_names{group_indx},'unassigned')
        xform_isbrain_intersect=ones(128,128); 
        tmp=[find(group_id==group_indx)]; %mice in the group
        save_name=strcat(Group_directory,filesep,group_names{group_indx}, '-avgPower') 

        if isfile(strcat(save_name,'.mat')) 
            disp(strcat(save_name,'.mat already exists, loading'));
            load(strcat(save_name,'.mat'));
            runInfo=runsInfo( first_ind_mouse(1));
        else
        
            for i= 1:numel(tmp)
                runInfo=runsInfo( first_ind_mouse(tmp(i)) );
                load(runInfo.saveMaskFile)
                xform_isbrain_intersect=xform_isbrain_intersect.*xform_isbrain;

                load_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgPower'); 
                load(strcat(load_name,'.mat'), 'whole_spectra_map_mouse_avg' , 'power_map_mouse_avg', 'global_sig_for_mouse_avg',...
                    'glob_sig_power_mouse_avg','hz')
                tmp_whole_spectra_map(:,:,:,:,i)=whole_spectra_map_mouse_avg;%.*noVascMask;
                tmp_powerMap(:,:,:,:,i)=power_map_mouse_avg;%.*noVascMask;
                tmp_global_sig_for(:,:,i)=global_sig_for_mouse_avg;%anmol- deleted .*noVascMask;
                tmp_glob_sig_power(:,:,i)=glob_sig_power_mouse_avg;%
            end
            whole_spectra_map_group_avg=squeeze(mean(tmp_whole_spectra_map,5,'omitnan'));%this is now size, 128,128,hz,contrast
            power_map_group_avg=squeeze(mean(tmp_powerMap,5,'omitnan'));%this is now size, 128,128,f-band,contrast
            global_sig_for_group_avg=squeeze(mean(tmp_global_sig_for,3,'omitnan')); %this is now size hz, contrast
            glob_sig_power_group_avg=squeeze(mean(tmp_glob_sig_power,3,'omitnan')); %this is now size hz, contrast

            avg_power_across_cortex_group_avg=reshape(whole_spectra_map_group_avg,[128*128,numel(hz),numel(runInfo.Contrasts)]); 
            avg_power_across_cortex_group_avg=squeeze(mean(avg_power_across_cortex_group_avg,1,'omitnan')); 

            save(strcat(save_name,'.mat'), 'whole_spectra_map_group_avg' , 'power_map_group_avg', 'global_sig_for_group_avg',...
                'glob_sig_power_group_avg','hz','avg_power_across_cortex_group_avg','xform_isbrain_intersect','-v7.3')        
        end
%PLOTTING     
    %% Non-overlay contrast
        fbands={'Full','ISA','Delta'};
        lims={[-120 -105],[-120 -105],[-130 -123]; ...
              [-120 -105],[-120 -105],[-130 -123]; ...
              [-120 -105],[-120 -105],[-130 -123]; ...
              [-50 -40],  [-60 -40],  [-60 -50]}; %Manually determined ranges for contrast and band
        
        
        for i=1:numel(runInfo.Contrasts)
            fh=figure;
            tiledlayout(1,3); %Anmol Note ('Padding','none','TileSpacing', 'compact')

            for j=1:3 %through each freq band
                nexttile
                imagesc(10*log10(power_map_group_avg(:,:,j,i)).*xform_isbrain_intersect,'alphadata',xform_isbrain_intersect)
                title(fbands{j},'Position',[62 15]) 
                axis image
                axis off
                colormap; 
                cb=colorbar;
                caxis(lims{i,j});      
            end
            
            if i<=3 
                cb.Label.String='10*log10(Mol^2)/Hz'; 
            else
                cb.Label.String='10*log10(\DeltaF/F%)^2/Hz';
            end         
            sgtitle([group_names{group_indx}, ' Average Power Map:' runInfo.Contrasts{i}])
            savefig(fh,strcat(save_name,'-MAP-',runInfo.Contrasts{i}));
            saveas(fh,strcat(save_name,'-MAP-',runInfo.Contrasts{i} ,'.png')); 
            
            
            fh2=figure;
            loglog(hz,global_sig_for_group_avg(:,i)/global_sig_for_group_avg(1,i)) 
            grid on
            xlabel('Frequency (Hz)')
            if i<=3
                ylabel('Normalized (Mol^2/Hz)')
            else
                ylabel('Normalized (\DeltaF/F%)^2/Hz')
            end
            ylim([10^-5 10])
            
            xlim([10^(-2) runsInfo(1).samplingRate/2]) 
            
            sgtitle([group_names{group_indx}, ' Global Signal Spectra:' runInfo.Contrasts{i}])        
            savefig(fh2,strcat(save_name,'-Glob_Sig_Spec-',runInfo.Contrasts{i}));%anmol- parenthesis error    


            fh3=figure;
            loglog(hz,avg_power_across_cortex_group_avg(:,i)/avg_power_across_cortex_group_avg(1,i));
            grid on
            xlabel('Frequency (Hz)')
            if i<=3
                ylabel('Normalized (Mol^2/Hz)')
            else
                ylabel('Normalized (\DeltaF/F%)^2/Hz')
            end
            grid on
            ylim([10^-4 10])
            xlim([10^(-2) runsInfo(1).samplingRate/2]) 

            sgtitle([group_names{group_indx}, ' Average Cortical Power:' runInfo.Contrasts{i}])        
            savefig(fh3,strcat(save_name,'-Avg_Cort_Power-',runInfo.Contrasts{i}));
            saveas(fh3,strcat(save_name,'-Avg_Cort_Power-',runInfo.Contrasts{i} ,'.png'));      
        end
        
    %% Overlay Contrasts on One Figure
        
%Global Signal
        fh22=figure;
        r=sysInfo(runInfo.system);
                
        for i=1:numel(runInfo.Contrasts)
            loglog(hz,global_sig_for_group_avg(:,i)/global_sig_for_group_avg(1,i),r.colors{i});
            hold on
        end
        
        if i<=3
        ylabel('Normalized (Mol^2/Hz)')
        else
        ylabel('Normalized (\DeltaF/F%)^2/Hz')            
        end
        legend(runInfo.Contrasts)
        sgtitle([group_names{group_indx}, ' Global Signal Spectra'])   
        xlabel('Frequency (Hz)')
        xlim([10^(-2) runsInfo(1).samplingRate/2]) %anmol - changed runsinfo to index first samplingrate
        ylim([10^-5 10^2])
        grid on

        savefig(fh22,strcat(save_name,'-Glob_Sig_Spec_ALL'));
        saveas(fh22,strcat(save_name,'-Glob_Sig_Spec_ALL.png'));
        close(fh22);
        
%Average Cortical Poewr

        fh33=figure;
        r=sysInfo(runInfo.system);
        
        for i=1:numel(runInfo.Contrasts)
            loglog(hz,avg_power_across_cortex_group_avg(:,i)/avg_power_across_cortex_group_avg(1,i),r.colors{i});
            hold on
        end
        grid on
        xlabel('Frequency (Hz)')
        if i<=3
        ylabel('Normalized (Mol^2/Hz)')
        else
        ylabel('Normalized (\DeltaF/F%)^2/Hz')            
        end
        grid on
        ylim([10^-4 10]) 
        xlim([10^(-2) runsInfo(1).samplingRate/2]) 
        legend(runInfo.Contrasts)
        sgtitle([group_names{group_indx}, ' Average Cortical Power'])   

        savefig(fh33,strcat(save_name,'-Avg_Cort_Power_ALL'));
        saveas(fh33,strcat(save_name,'-Avg_Cort_Power_ALL.png'));  
        
        

        end
            close all;

    
    end
end