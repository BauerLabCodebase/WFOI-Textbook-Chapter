%Cross Mag Spectrun
excelFile="C:\Users\Nischal\Desktop\code mod\newpipeline\Code Modification\exampleTiffOIS+Gcamp.xlsx";
excelRows=8:12;

load("A:\ToJonah\190627-R5M2285-fc1_processed.mat",'xform_datahb','xform_jrgeco1aCorr','runInfo','xform_isbrain')
runInfo.samplingRate=25;
xform_datafluorCorr=xform_jrgeco1aCorr;

%Loading Info
runsInfo = parseRuns(excelFile,excelRows);
runNum=numel(runsInfo);
load(which('AtlasandIsbrain.mat'),'AtlasSeeds_Big')
Seed_names={'olf','fr','cing','m','ss','rs','p','vis','aud','ass'};
AtlasSeeds_Big(AtlasSeeds_Big==0)=nan;
 
% Loading Variables
    for runInd = 1:runNum
        runInfo=runsInfo(runInd);
        load(runInfo.saveMaskFile,'xform_isbrain');
        %make mask nan
        xform_isbrain=single(xform_isbrain);
        tmp=find(xform_isbrain==0);
        xform_isbrain(tmp)=nan;

        
        %load Hb and Calc. 
        load(runInfo.saveHbFile,'xform_datahb');
                xform_total = squeeze(xform_datahb(:,:,1,:)+ xform_datahb(:,:,2,:));
                xform_total(isinf(xform_total)) = 0;
                xform_total(isnan(xform_total)) = 0;
       load(runInfo.saveFluorFile,'xform_datafluorCorr') 
                xform_datafluorCorr(isinf(xform_datafluorCorr)) = 0;
                xform_datafluorCorr(isnan(xform_datafluorCorr)) = 0;                

                
for dx=1:size(xform_datafluorCorr,1)
    for dy=1:size(xform_datafluorCorr,2)
        a=squeeze(xform_total(dx,dy,:));
        b=squeeze(xform_datafluorCorr(dx,dy,:));
[MagSqCoher_Calc_HbT(dx,dy,:),freq(dx,dy,:)]=mscohere(a,b,[],[],[],runInfo.samplingRate  );
[CrossSpec_Calc_HbT(dx,dy,:),freq(dx,dy,:)] = cpsd(a,b,[],[],[],runInfo.samplingRate  ); 
tmp=find(squeeze(MagSqCoher_Calc_HbT(dx,dy,:)) < 0.2);

CrossSpec_Calc_HbT(dx,dy,tmp) = 0; %Neglect cross spectrum when coherence is small.

CrossSpecPhase(dx,dy,:)=angle(CrossSpec_Calc_HbT(dx,dy,:)/pi);

    end
end


                
%Rotating Things 
xform_isbrain_rot=flipdim(imrotate(xform_isbrain,-90),2);

CrossSpecPhase_rotate1=flipdim(imrotate(CrossSpecPhase,-90),2);
CrossSpecPhase_rotate=reshape(CrossSpecPhase_rotate1,[size(CrossSpecPhase_rotate1,1)*size(CrossSpecPhase_rotate1,2)],size(CrossSpecPhase_rotate1,3));


MagSqCoher_Calc_HbT_rotate1=flipdim(imrotate(MagSqCoher_Calc_HbT,-90),2);
MagSqCoher_Calc_HbT_rotate=reshape(MagSqCoher_Calc_HbT_rotate1,[size(MagSqCoher_Calc_HbT_rotate1,1)*size(MagSqCoher_Calc_HbT_rotate1,2)],size(MagSqCoher_Calc_HbT_rotate1,3));


    keep_rot=[];
    indx=cell(size(Seed_names,2),1);
    indx_rot=cell(size(Seed_names,2),1);
    sort_ind_rot=[];    
    a_rot=0;      
    
    %only consider things that are on the brain map
    for i=1:size(Seed_names,2)       
        indx_rot{i}=find((flipdim(imrotate(AtlasSeeds_Big.*xform_isbrain,-90),2)==i)==1);  
        [x{i},y{i}]=ind2sub([160,160],indx_rot{i});
        keep_rot=[keep_rot;indx_rot{i}]; %if we want all L and then all right  
        a_rot=a_rot+length(indx_rot{i});
        sort_ind_rot=[sort_ind_rot; a_rot]; 
        
    end
           
    
    
%Plotting
tiledlayout(2,2)
nexttile

%Mag-Sq Cohernce over Cortex
imagesc(MagSqCoher_Calc_HbT_rotate(keep_rot,:))
colormap default
h=colorbar;
set(get(h,'title'),'string','Coherence');
hold on
for i=1:length(sort_ind_rot)
line([1 6000], [sort_ind_rot(i),sort_ind_rot(i)], 'Color', 'r');
end
ylim([sort_ind_rot(1)+1 sort_ind_rot(10)])
yticks(sort_ind_rot(1:end))
yticklabels(Seed_names(1:end))
ytickangle(35)
caxis([0 1])
title('Magnitude-Square Coherence')
xlabel('Frequnecy (Hz)')
xticks([0:size(MagSqCoher_Calc_HbT_rotate,2)/5:size(MagSqCoher_Calc_HbT_rotate,2)])
xtickz=xticks;
xticklabels([xtickz.*runInfo.samplingRate./xtickz(end)]')
               
    
%Average Mag-Sq Cohernce over Cortex
nexttile
plot(mean(MagSqCoher_Calc_HbT_rotate(keep_rot,:)))
xticks([0:size(MagSqCoher_Calc_HbT_rotate,2)/5:size(MagSqCoher_Calc_HbT_rotate,2)])
xtickz=xticks;
xticklabels([xtickz.*runInfo.samplingRate./xtickz(end)]')
title('Average Over Cortex:Magnitude-Square Coherence')
xlabel('Frequnecy (Hz)')
grid on



%Cross Spectrum Phase Over Cortex
nexttile
imagesc(CrossSpecPhase_rotate(keep_rot,:))
colormap default
h=colorbar;
set(get(h,'title'),'string','Lag (\times\pi rad)');
hold on
for i=1:length(sort_ind_rot)
line([1 6000], [sort_ind_rot(i),sort_ind_rot(i)], 'Color', 'r');
end
ylim([sort_ind_rot(1)+1 sort_ind_rot(10)])
yticks(sort_ind_rot(1:end))
yticklabels(Seed_names(1:end))
ytickangle(35)
caxis([-pi pi])
title('Cross Spectrum Phase Over Cortex')
xlabel('Frequnecy (Hz)')
xticks([0:size(MagSqCoher_Calc_HbT_rotate,2)/5:size(MagSqCoher_Calc_HbT_rotate,2)])
xtickz=xticks;
xticklabels([xtickz.*runInfo.samplingRate./xtickz(end)]')
               
    
%Average Mag-Sq Cohernce over Cortex
nexttile
plot(mean(CrossSpecPhase_rotate(keep_rot,:)))
xticks([0:size(MagSqCoher_Calc_HbT_rotate,2)/5:size(MagSqCoher_Calc_HbT_rotate,2)])
xtickz=xticks;
xticklabels([xtickz.*runInfo.samplingRate./xtickz(end)]')
title('Average Over Cortex:Cross Spectrum Phase Over Cortex')
xlabel('Frequnecy (Hz)')
grid on


    end






