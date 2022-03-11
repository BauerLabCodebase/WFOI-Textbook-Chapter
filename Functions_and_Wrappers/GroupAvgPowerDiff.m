function GroupAvgPowerDiff(excelFile,excelRows,Group_directory)

    %loading wl params
    %no vasc mask
    noVascMask=ones(128,128);
    noVascMask(:,61:67)=0;

    %getting excel/mouse info
    runsInfo = parseRuns(excelFile,excelRows);
    [excel_row,first_ind_mouse,~]=unique({runsInfo.excelRow_char});
    excel_row=str2double(excel_row); %anmol- semicolon
    excelData = readtable(excelFile);
    tableRows = excel_row - 1;

    %Get Group labels and ignore, non-labled mice
    groupLabel=cell(numel(excel_row),1);
    groupLabelsOnly=cell(numel(excel_row),1);
    ct=1; %anmol- added iterator
    for row = tableRows
         
     if isnumeric(excelData{row,'Timepoint'}) & ~isnan(excelData{row,'Timepoint'})
        tp=num2str(excelData{row,'Timepoint'});
     elseif isnan(excelData{row,'Timepoint'})
         tp='-';
     else
        tp=excelData{row,'Timepoint'};
    end
    
    groupLabel(ct) = strcat(excelData{row,'Group'},'-',tp);
        groupLabelsOnly(ct)=excelData{row,'Group'};
        if isempty(groupLabel{ct})
            groupLabel(ct)={'unassigned'};
            warning(['row ' num2str(row) ' has no group assignment. It will not be included.'])
        end
        ct=ct+1;
    end

    [group_names,~,group_id]=unique(groupLabel);
    [gn,~,group_ido]=unique(groupLabelsOnly);
    
    
    qq=cellfun(@(s) strsplit(s,'-'),group_names,'UniformOutput',false);

    %Establish group directory
    if ~exist(Group_directory)
        mkdir(Group_directory)
    end
    
   

    %FOR EACH GROUP
    for rgroup=1:size(gn,1)
        runInfo=runsInfo( first_ind_mouse(1)); %load for the if statement -  assumes every type of mouse in a group is the same
      
        corre=cellfun(@(u) contains(u,gn{rgroup}),group_names);
        onegroup=group_names(corre);
        splitnames=cellfun(@(u) strsplit(u,'-'),onegroup,'UniformOutput',false);
        
        if numel(onegroup)~=1
            
            
            
            cMask=ones([128 128]); 
            for group_indx=1:numel(onegroup)
                load_name=strcat(Group_directory,filesep,onegroup{group_indx}, '-avgPower');
                tb=load(load_name,'xform_isbrain_intersect'); %anmol - union to intersect
                cMask=tb.xform_isbrain_intersect.*cMask;
            end

                
            

            for group_indx=1:numel(onegroup)-1 %for each timepoint

                load_name=strcat(Group_directory,filesep,onegroup{group_indx}, '-avgPower');
                load_name1=strcat(Group_directory,filesep,onegroup{group_indx+1}, '-avgPower');
                
                
                first=load(strcat(load_name,'.mat'),'whole_spectra_map_group_avg' , 'power_map_group_avg', 'global_sig_for_group_avg',...
            'glob_sig_power_group_avg','hz','avg_power_across_cortex_group_avg','xform_isbrain_union');
                second=load(strcat(load_name1,'.mat'),'whole_spectra_map_group_avg' , 'power_map_group_avg', 'global_sig_for_group_avg',...
            'glob_sig_power_group_avg','hz','avg_power_across_cortex_group_avg','xform_isbrain_union');
        
                power_map_group_sub=second.power_map_group_avg-first.power_map_group_avg;
                global_sig_for_group_sub=second.global_sig_for_group_avg-first.global_sig_for_group_avg;
                glob_sig_power_group_sub=second.glob_sig_power_group_avg-first.glob_sig_power_group_avg;
                avg_power_across_cortex_group_sub=second.avg_power_across_cortex_group_avg-first.avg_power_across_cortex_group_avg;
                


                % insert normal function stuff here 
                
                fbands={'Full','ISA','Delta'};
                lims={[-2 2],[-2 2],[-2 2]; ...
                  [-2 2],[-2 2],[-2 2]; ...
                  [-2 2],[-2 2],[-2 2]; ...
                  [-2.5 2.5],[-2.5 2.5],[-2.5 2.5]}; % anmol-
                
                
                for i=1:numel(runInfo.Contrasts)
                    fh=figure;
                    tiledlayout(1,3)

                    for j=1:3 %through each freq band
                        nexttile
                        noodle=10*log10(second.power_map_group_avg(:,:,j,i))-10*log10(first.power_map_group_avg(:,:,j,i)); %anmol- no normalization here!

                        imagesc(noodle,'alphadata',cMask)
                        title(fbands{j},'Position',[62 15]) 
                        axis image
                        axis off
                        colormap('redblue');
                        caxis(lims{i,j});
                    end
                    cb=colorbar;
                    if i<=3 %anmol - changed i =< 3 to i <= 3
                        %cb.Label.String='Normalized Log10(Mol^2/Hz)';
                        cb.Label.String='\Delta log10(Mol^2)'; %anmol - changed unitz 
                    else
                        %cb.Label.String='Normalized Log10(\DeltaF/F%)^2/Hz';
                        cb.Label.String='\Delta log10((\DeltaF/F%)^2)'; %anmol - no normalization here
                    end         
                    
                    letitle=strcat(splitnames{group_indx}{1},{' '},splitnames{group_indx+1}{2},'-',splitnames{group_indx}{2},' power-MAP-diff:',runInfo.Contrasts{i}); %anmol - minor title formatting
                    save_name=strcat(splitnames{group_indx}{1},{' '},splitnames{group_indx+1}{2},'-',splitnames{group_indx}{2},' power-MAP-diff-',runInfo.Contrasts{i}); %anmol - minor title formatting

                    
                    sgtitle(letitle)
                    savefig(fh,strcat(Group_directory,filesep,save_name));
                    saveas(fh,strcat(Group_directory,filesep,save_name,'.png')); 
                    close all
                end
            
                % end normal function stuff
            end
            
            
            
%last point vs first point
            load_name=strcat(Group_directory,filesep,onegroup{1}, '-avgPower');
            load_name1=strcat(Group_directory,filesep,onegroup{numel(onegroup)}, '-avgPower');
            save_name=strcat(Group_directory,filesep, splitnames{1}{1},' ',splitnames{numel(onegroup)}{2},'-',splitnames{1}{2},' lag-diff');
   
            first=load(strcat(load_name,'.mat'),'whole_spectra_map_group_avg' , 'power_map_group_avg', 'global_sig_for_group_avg',...
            'glob_sig_power_group_avg','hz','avg_power_across_cortex_group_avg','xform_isbrain_union');
            second=load(strcat(load_name1,'.mat'),'whole_spectra_map_group_avg' , 'power_map_group_avg', 'global_sig_for_group_avg',...
            'glob_sig_power_group_avg','hz','avg_power_across_cortex_group_avg','xform_isbrain_union');
        
            power_map_group_sub=second.power_map_group_avg-first.power_map_group_avg;
            global_sig_for_group_sub=second.global_sig_for_group_avg-first.global_sig_for_group_avg;
            glob_sig_power_group_sub=second.glob_sig_power_group_avg-first.glob_sig_power_group_avg;
            avg_power_across_cortex_group_sub=second.avg_power_across_cortex_group_avg-first.avg_power_across_cortex_group_avg;
  

           fbands={'Full','ISA','Delta'};    
                for i=1:numel(runInfo.Contrasts)
                    fh=figure;
                    tiledlayout(1,3)

                    for j=1:3 %through each freq band
                        nexttile
                        
%                         noodle=10*log10(second.power_map_group_avg(:,:,j,i).*1/max(max(second.power_map_group_avg(:,:,j,i).*1)))- ...
%                             10*log10(first.power_map_group_avg(:,:,j,i).*1/max(max(first.power_map_group_avg(:,:,j,i).*1)));
                        noodle=10*log10(second.power_map_group_avg(:,:,j,i))-10*log10(first.power_map_group_avg(:,:,j,i)); %anmol- no normalization here!
                        
                        
                        imagesc(noodle,'alphadata',cMask)
                        title(fbands{j},'Position',[62 15]) 
                        axis image
                        axis off
                        colormap('redblue');
                        caxis(lims{i,j});
                    end
                    cb=colorbar;
                    if i<=3 %anmol - changed i =< 3 to i <= 3
                        %cb.Label.String='Normalized Log10(Mol^2/Hz)';
                        cb.Label.String='Mol^2'; 
                    else
                        %cb.Label.String='Normalized Log10(\DeltaF/F%)^2/Hz';
                        cb.Label.String='\Delta log10((\DeltaF/F%)^2)'; %anmol - changed units
                    end         
                    
                    letitle=strcat(splitnames{1}{1},{' '},splitnames{numel(onegroup)}{2},'-',splitnames{1}{2},' power-MAP-diff:',runInfo.Contrasts{i});
                    savename=strcat(splitnames{1}{1},{' '},splitnames{numel(onegroup)}{2},'-',splitnames{1}{2},' power-MAP-diff-',runInfo.Contrasts{i}); %anmol-savename=/=title
                    
                    sgtitle(letitle)
                    savefig(fh,strcat(Group_directory,filesep,savename));
                    saveas(fh,strcat(Group_directory,filesep,savename,'.png')); 
                    close all;
                end
    

        elseif numel(onegroup)==1
            disp('only one group found')
            continue
        end

    end
    
end
    
    