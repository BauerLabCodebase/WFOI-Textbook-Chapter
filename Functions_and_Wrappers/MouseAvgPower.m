%%Get list of mice and runs, and load into a cancat file
function MouseAvgPower(excelFile,excelRows)
%ALL of these should be ISA ( 0.01-0.08) and delta *(.5 to 4) because it
%uses JPC power_wrapper
noVascMask=ones(128,128);
noVascMask(:,61:67)=0;

runsInfo = parseRuns(excelFile,excelRows);
[excel_row,first_ind_mouse,mouse_ids]=unique({runsInfo.excelRow_char});

%Mouse Averages
for mouse_indx=1:length(first_ind_mouse) %this is the number of mice
    inds=find(mouse_ids==mouse_indx);  %how many runs for each mice and which one is the associated runsInfo?  
    runInfo=runsInfo(first_ind_mouse(mouse_indx));      
    save_name=strcat(runInfo.saveFolder ,filesep, runInfo.recDate, '-' ,runInfo.mouseName, '-avgPower'); 
%     if isfile(strcat(save_name,'.mat')) 
%         disp([strcat(save_name,'.mat') ' already exists'])
%         continue 
%     end
    load(runInfo.saveMaskFile,'xform_isbrain'); 

    for mouse_run=1:length(inds) %load for each run
            ind=inds(mouse_run);   %associated index in runsInfo
            load([runsInfo(ind).saveFilePrefix,'-Power','.mat']); %load

            tmp_whole_spectra_map(:,:,:,:,mouse_run)=whole_spectra_map.*noVascMask;
            tmp_powerMap(:,:,:,:,mouse_run)=powerMap.*noVascMask;
            tmp_global_sig_for(:,:,mouse_run)=global_sig_for;%anmol-deleted *.Novascmask
            tmp_glob_sig_power(:,:,mouse_run)=glob_sig_power;%.*noVascMask;
            
    end 

    disp(['Loaded ', runInfo.mouseName, ' runs: ' num2str(ind)])

    whole_spectra_map_mouse_avg=squeeze(mean(tmp_whole_spectra_map,5,'omitnan'));%this is now size, 128,128,hz,contrast
    power_map_mouse_avg=squeeze(mean(tmp_powerMap,5,'omitnan'));%this is now size, 128,128,f-band,contrast
    global_sig_for_mouse_avg=squeeze(mean(tmp_global_sig_for,3,'omitnan')); %this is now size hz, contrast
    glob_sig_power_mouse_avg=squeeze(mean(tmp_glob_sig_power,3,'omitnan')); %this is now size hz, contrast
        
    avg_power_across_cortex_Mouse_avg=reshape(whole_spectra_map_mouse_avg,[128*128,numel(hz),numel(runInfo.Contrasts)]);
    avg_power_across_cortex_Mouse_avg=squeeze(mean(avg_power_across_cortex_Mouse_avg,1,'omitnan')); %this is now size hz, contrast %removed a parenthesis
%Saving    
    save(strcat(save_name,'.mat'), 'whole_spectra_map_mouse_avg' , 'power_map_mouse_avg', 'global_sig_for_mouse_avg',...
    'glob_sig_power_mouse_avg','hz','avg_power_across_cortex_Mouse_avg','-v7.3'); 

%Plotting
fbands={'Full','ISA','Delta'};    
lims={[-120 -105],[-120 -105],[-130 -123]; ...
              [-120 -105],[-120 -105],[-130 -123]; ...
              [-120 -105],[-120 -105],[-130 -123]; ...
              [-50 -40],  [-60 -40],  [-60 -50]}; %manually determined ranges for contrast and band
        
        
    for i=1:numel(runInfo.Contrasts)
        fh1=figure; 
        tiledlayout(1,3)
        
        for j=1:3 %through each freq band
            nexttile
            imagesc(10*log10(power_map_mouse_avg(:,:,j,i)).*xform_isbrain,'alphadata',xform_isbrain)% Xiaodan delete _intersect
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
            cb.Label.String='10*log10(\DeltaF/F)^2/Hz';%Xiaodan deleted %
        end
        sgtitle([runInfo.mouseName, ' Average Power Map:' runInfo.Contrasts{i}])
        savefig(fh1,strcat(save_name,'-MAP-',runInfo.Contrasts{i})); 
        saveas(fh1,strcat(save_name,'-MAP-',runInfo.Contrasts{i} ,'.png'));  

        
        fh2=figure;
        loglog(hz,global_sig_for_mouse_avg(:,i)/global_sig_for_mouse_avg(1,i)); %Normalized by power at first datapoint
        grid on
        xlabel('Frequency (Hz)')
        if i<=3
        ylabel('Normalized (Mol^2/Hz)')
        else
        ylabel('Normalized (\DeltaF/F)^2/Hz')%Xiaodan deleted %            
        end
        ylim([10^-5 10])
        xlim([10^(-2) runInfo.samplingRate/2]) %Xiaodan change runsInfo to runInfo
        
        sgtitle([runInfo.mouseName, ' Global Signal Spectra:' runInfo.Contrasts{i}])        
        savefig(fh2,strcat(save_name,'-Glob_Sig_Spec-',runInfo.Contrasts{i}));
        saveas(fh2,strcat(save_name,'-Glob_Sig_Spec-',runInfo.Contrasts{i} ,'.png')); 

        fh3=figure();
        loglog(hz,avg_power_across_cortex_Mouse_avg(:,i)/avg_power_across_cortex_Mouse_avg(1,i));
        grid on
        xlabel('Frequency (Hz)')
        if i<=3
        ylabel('Normalized (Mol^2/Hz)')
        else
        ylabel('Normalized (\DeltaF/F)^2/Hz')  %Xiaodan deleted %          
        end
        grid on
        ylim([10^-4 10]) 
        xlim([10^(-2) runInfo.samplingRate/2]) %Xiaodan change runsInfo to runInfo
        
        sgtitle([runInfo.mouseName, ' Average Cortical Power:' runInfo.Contrasts{i}])
        savefig(fh3,strcat(save_name,'-Avg_Cort_Power-',runInfo.Contrasts{i}));
        saveas(fh3,strcat(save_name,'-Avg_Cort_Power-',runInfo.Contrasts{i} ,'.png'));            
    end
    close all; 

    
end