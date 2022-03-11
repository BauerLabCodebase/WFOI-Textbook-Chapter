function GammaFit_wrapper(excelFile,excelRows)
%%Parallel Pool
poolobj = gcp('nocreate'); % If no pool, do not create new one.
numcores = feature('numcores');
if isempty(poolobj)
    parpool('local',numcores-1);
end

%Loading Info
    runsInfo = parseRuns(excelFile,excelRows);
    [row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); 
    runNum=numel(runsInfo);

    load(which('GoodWL.mat'),'xform_WL'); 
    load(which('noVasculatureMask.mat')); 

% Processing 
    for runInd = 1:runNum
        runInfo=runsInfo(runInd);
        savename = fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session,num2str(runInfo.run)','_gamma')); %anmol- added a savename here to check if it has already run
        if isfile(strcat(savename,'.mat')) 
            disp([strcat(savename,'.mat') ' already exists'])
            continue 
        end
        
    disp('loading processed data')
        load(runInfo.saveMaskFile,'xform_isbrain');
        %load Hb
        load(runInfo.saveHbFile,'xform_datahb');
                xform_total = squeeze(xform_datahb(:,:,1,:)+ xform_datahb(:,:,2,:));
                xform_total(isinf(xform_total)) = 0;
                xform_total(isnan(xform_total)) = 0;    
       %Load Calcium
       load(runInfo.saveFluorFile,'xform_datafluorCorr') 
                xform_datafluorCorr(isinf(xform_datafluorCorr)) = 0;
                xform_datafluorCorr(isnan(xform_datafluorCorr)) = 0;    
%        %Load FAD
%         if ~isempty(runInfo.FADChInd)
%             load(runInfo.saveFluorFile,'xform_dataFADCorr')
%                 xform_dataFADCorr(isnan(xform_dataFADCorr)) = 0;
%                 xform_dataFADCorr(isinf(xform_dataFADCorr)) = 0;            
%         end
        OGsize=size(xform_datafluorCorr);
    disp('filtering')
%Filtering    
        onlyBrain=find(xform_isbrain==1);

        xform_datafluorCorr = filterData(double(xform_datafluorCorr),0.02,2,runInfo.samplingRate);
        xform_datafluorCorr = reshape(xform_datafluorCorr,OGsize(1)*OGsize(2),[]);
        xform_datafluorCorr=resample(xform_datafluorCorr',5,runInfo.samplingRate)';
        
        xform_total = filterData(double(xform_total),0.02,2,runInfo.samplingRate);% a 0.02-Hz high-pass filter (HPF) to remove slow drifts, as well as a 2-Hz low-pass filter (LPF) to reduce physiological noise
        xform_total = reshape(xform_total,OGsize(1)*OGsize(2),[]);
        xform_total=resample(xform_total',5,runInfo.samplingRate)';
        
        xform_datafluorCorr = normRow(xform_datafluorCorr);
        xform_total = normRow(xform_total);
        
        xform_datafluorCorr = reshape(xform_datafluorCorr,OGsize(1),OGsize(2),[]);
        xform_total = reshape(xform_total,OGsize(1),OGsize(2),[]);

%setting things up
        time_epoch=30;
        t= [0:1/runInfo.samplingRate:time_epoch]; %30 seconds7
        t=linspace(0,time_epoch,time_epoch*runInfo.samplingRate *(5/runInfo.samplingRate));%% force it to be 5 hz
       
     tic   
    [T,W,A,r,r2,hemoPred,objective_vals] = GammaFit(xform_datafluorCorr,xform_total,t,xform_isbrain);
    time_taken=toc
    save(fullfile(savename,'.mat'),'T','W','A','r','r2','hemoPred','-append')

end




end