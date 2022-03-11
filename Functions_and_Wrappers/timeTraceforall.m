function fh = timeTraceforall(data,runInfo,ROI,colors,contrastName,goodBlocks)
    data = squeeze(mean(data(:,:,:,goodBlocks,:),4));
    ibi = reshape(ROI,1,[]);
    data = reshape(data,size(data,1)*size(data,2),size(data,3),size(data,4));
    data_ROI = squeeze(mean(data(ibi,:,:),1));
    time = (1:size(data_ROI,1))/runInfo.samplingRate;
    cMax = nan(1,3);
    fh = figure;
    for ii = 1:3
        cMax(ii) = max(abs(data_ROI(:,ii)),[],'all');
    end
    cMax = max(cMax,[],'all');
    for ii = 1:3
        plot(time,data_ROI(:,ii),colors{ii})
        hold on 
    end
    xlabel('Time(s)','FontSize',12)
    ylabel('Hb(\Delta\muM)','FontSize',12)
    ylim([-1.1*cMax 1.1*cMax])
    
    
    if ~isempty(runInfo.fluorChInd)
        fMax = max(data_ROI(:,4),[],'all');
        yyaxis right
        plot(time,data_ROI(:,4),colors{4})
        ylabel('Fluorescence(\DeltaF/F)','FontSize',12);
        ylim([-1.1*fMax 1.1*fMax])
        
    end
    legend(contrastName)
    set(findall(gca, 'Type', 'Line'),'LineWidth',1);