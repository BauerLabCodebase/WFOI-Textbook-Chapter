%% Get Run and System Info
function Lag_analysis_wrapper(excelFile,excelRows)
    % excelFile="C:\Users\xiaodanwang\Documents\GitHub\BauerLab\MATLAB\examples\Code Modification\exampleTiffOIS+Gcamp.xlsx";
    % excelRows=[15];  % Rows from Excell Database
    
%Loading Info
    runsInfo = parseRuns(excelFile,excelRows);
    [row,start_ind_mouse,numruns_per_mouse]=unique({runsInfo.excelRow_char}); 
    runNum=numel(runsInfo);

    load(which('GoodWL.mat'),'xform_WL'); 
    load(which('noVasculatureMask.mat')); 

% Processing 
    for runInd = 1:runNum
        runInfo=runsInfo(runInd);
        savename = fullfile(runInfo.saveFolder,strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session,num2str(runInfo.run)','_lag')); %anmol- added a savename here to check if it has already run
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
       %Load FAD
        if ~isempty(runInfo.FADChInd)
            load(runInfo.saveFluorFile,'xform_dataFADCorr')
                xform_dataFADCorr(isnan(xform_dataFADCorr)) = 0;
                xform_dataFADCorr(isinf(xform_dataFADCorr)) = 0;            
        end% % 
    disp('filtering')
%Filtering    
        xform_total = filterData(double(xform_total),0.02,2,runInfo.samplingRate);% a 0.02-Hz high-pass filter (HPF) to remove slow drifts, as well as a 2-Hz low-pass filter (LPF) to reduce physiological noise
        xform_datafluorCorr = filterData(double(xform_datafluorCorr),0.02,2,runInfo.samplingRate);
        if ~isempty(runInfo.FADChInd)
                xform_dataFADCorr = filterData(double(xform_dataFADCorr),0.02,2,runInfo.samplingRate);
        end
        
%Setting up for Lag Algorithm
    %WHY IS THIS HERE AND HOW CAN WE NOT HAVE USER INPUT HALF WAY INTO OUR CODE?
        edgeLen =1; % this means how many points to the left and right of the smapled maximum
        tZone = 4; %how many seconds to look at the lag (only minus) when use xcorr
        corrThr = 0; %if the maximum correlation is less than 0 at that pixel is less than this then eliminate (only look at positive correlation)
        validRange = - edgeLen: round(tZone*runInfo.samplingRate); %     How many data points to consider from zero lag. For example, -3:3 would mean looking at 7 points around zero lag.
        tLim = [-2 2]; %for visualizations for calcium and hb (calcium comes first)
        tLim_FAD = [0 0.3]; % %for visualizations
        rLim = [-1 1]; %for visualizations


    disp(strcat('Lag analysis on ', runInfo.recDate, '', runInfo.mouseName, ' run#', num2str(runInfo.run)))
    %
        [lagTimeTrial_HbTCalcium, lagAmpTrial_HbTCalcium,covResult_HbTCalcium] = dotLag(...
            xform_total,xform_datafluorCorr,edgeLen,validRange,corrThr, true,true);
        lagTimeTrial_HbTCalcium = lagTimeTrial_HbTCalcium./runInfo.samplingRate; %convert to seconds
%Plotting
        if isempty(runInfo.FADChInd) %If we dont have FAD
            clear xform_datahb xform_datafluorCorr            
            %Plotting                 
            figure;
            
            subplot(2,1,1) 
            imagesc(lagTimeTrial_HbTCalcium,tLim) 
            axis image off
            h = colorbar;ylabel(h,'t(s)')
            title('Calcium HbT')
            hold on
            colormap jet            
            imagesc(xform_WL,'AlphaData',1-mask_new)
            
            subplot(2,1,2)
            imagesc(lagAmpTrial_HbTCalcium,rLim)
            axis image off
            h = colorbar
            ylabel(h,'r')
            hold on
            imagesc(xform_WL,'AlphaData',1-mask_new)
            sgtitle(strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session,num2str(runInfo.run)))
            colormap jet            
            %saving
            saveas(gcf,[savename '.png']);
            saveas(gcf,[savename '.fig']);
            save([savename,'.mat']...
                ,'lagTimeTrial_HbTCalcium', 'lagAmpTrial_HbTCalcium', '-v7.3');
            close all
        else %if we have FAD

            [lagTimeTrial_FADCalcium, lagAmpTrial_FADCalcium,covResult_FADCalcium] = dotLag(...
                xform_dataFADCorr,xform_datafluorCorr,edgeLen,validRange,corrThr, true,true);
            lagTimeTrial_FADCalcium= lagTimeTrial_FADCalcium./runInfo.samplingRate;

            [lagTimeTrial_HbTFAD, lagAmpTrial_HbTFAD,covResult_HbTFAD] = dotLag(...
                xform_total,xform_dataFADCorr,edgeLen,validRange,corrThr, true,true);
            lagTimeTrial_HbTFAD= lagTimeTrial_HbTFAD./runInfo.samplingRate;


            clear xform_datahb xform_datafluorCorr xform_dataFADCorr

            figure;
            subplot(2,3,1)
            imagesc(lagTimeTrial_HbTCalcium,tLim)
            axis image off
            h = colorbar;ylabel(h,'t(s)')
            title('Calcium HbT')
            hold on
            colormap jet;            
            imagesc(xform_WL,'AlphaData',1-mask_new)
            
            subplot(2,3,2)
            imagesc(lagTimeTrial_FADCalcium,tLim_FAD)
            axis image off
            h = colorbar
            ylabel(h,'t(s)')
            title('FAD Calcium')
            hold on
            imagesc(xform_WL,'AlphaData',1-mask_new)
            
            subplot(2,3,3)
            imagesc(lagTimeTrial_HbTFAD,[-1 1])
            axis image off
            h = colorbar;ylabel(h,'t(s)')
            title('FAD HbT')
            hold on
            imagesc(xform_WL,'AlphaData',1-mask_new)

            subplot(2,3,4)
            imagesc(lagAmpTrial_HbTCalcium,rLim)
            axis image off
            h = colorbar
            ylabel(h,'r')
            hold on
            imagesc(xform_WL,'AlphaData',1-mask_new)
            
            subplot(2,3,5)
            imagesc(lagAmpTrial_FADCalcium,rLim)
            axis image off
            h = colorbar
            ylabel(h,'r')
            hold on
            imagesc(xform_WL,'AlphaData',1-mask_new)
            
            subplot(2,3,6)
            imagesc(lagAmpTrial_HbTFAD,rLim)
            axis image off
            h = colorbar
            ylabel(h,'r')
            hold on
            imagesc(xform_WL,'AlphaData',1-mask_new)
            
            sgtitle(strcat(runInfo.recDate,'-',runInfo.mouseName,'-',runInfo.session,num2str(runInfo.run))) 

            saveas(gcf,fullfile(savename,'.png'));
            saveas(gcf,fullfile(savename,'.fig'));

            save(fullfile(savename,'.mat')...
                ,'lagTimeTrial_HbTCalcium', 'lagAmpTrial_HbTCalcium','lagTimeTrial_FADCalcium', 'lagAmpTrial_FADCalcium','lagTimeTrial_HbTFAD','lagAmpTrial_HbTFAD','-v7.3');
            close all

        end                                
    end
end

