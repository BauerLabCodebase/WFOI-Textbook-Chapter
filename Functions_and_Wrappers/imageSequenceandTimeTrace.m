function fh = imageSequenceandTimeTrace(data,runInfo,xform_isbrain,goodBlocks,contrastName,ROI,colors,ylabelName)
fh = figure('units','normalized','outerposition',[0 0 1 1]);
subplot('position',[0.05,0.08,0.55,0.35])

data = squeeze(mean(data(:,:,:,goodBlocks,:),4));
ibi = reshape(ROI,1,[]);
data = reshape(data,size(data,1)*size(data,2),size(data,3),size(data,4));
data_ROI = squeeze(nanmean(data(ibi,:,:),1));
time = (1:size(data_ROI,1))/runInfo.samplingRate;
m = nan(1,length(contrastName));
for ii = 1:length(contrastName)
    cMax_temp(ii) = max(abs(data_ROI(:,ii)),[],'all');
end
cMax = max(cMax_temp,[],'all');
for ii = 1:length(contrastName)
    plot(time,data_ROI(:,ii),colors{ii})
    hold on
end
xlabel('Time(s)','FontSize',12)
ylabel(ylabelName,'FontSize',12)
ylim([-1.1*cMax 1.1*cMax])

data = reshape(data,size(xform_isbrain,1),size(xform_isbrain,2),size(data,2),size(data,3));
downSampled = nan(size(xform_isbrain,1),size(xform_isbrain,2),size(data,3)/runInfo.samplingRate*1,size(data,4));
for ii = 1:size(data,4)
    downSampled(:,:,:,ii) = resampledata(squeeze(data(:,:,:,ii)),runInfo.samplingRate,1);
end

load('C:\Users\xiaodanwang\Documents\GitHub\BauerLabXiaodanScripts\GoodWL','WL')
for ii = 1:size(data,4)
    for b=runInfo.stimStartTime-2: runInfo.stimEndTime+5
        p = subplot('position', [0.16+(b-runInfo.stimStartTime)*0.07 0.80-(ii-1)*0.16 0.07 0.12]);
        imagesc(downSampled(:,:,b,ii), [-1.1*cMax_temp(ii) 1.1*cMax_temp(ii)]);
        if b == runInfo.stimEndTime+5
            colorbar
            set(p,'Position',[0.16+(b-runInfo.stimStartTime)*0.07 0.80-(ii-1)*0.16 0.07 0.12]);
        end
        axis image
        set(gca, 'XTick', []);
        set(gca, 'YTick', []);
        if b == runInfo.stimStartTime-2
            ylabel(contrastName{ii},'FontWeight','bold')
        end
        hold on
        imagesc(WL,'AlphaData',1-xform_isbrain)
        
        if ii ==1
            title([num2str(b),'s']);
            ax = gca;
            ax.FontSize = 10;
            
        end
        colormap jet
    end
end

subplot('position',[0.64,0.27,0.13,0.18])
imagesc(mean(downSampled(:,:,runInfo.stimStartTime+2:runInfo.stimEndTime,1),3),...
    [-cMax_temp(1)*0.8,cMax_temp(1)*0.8])
colorbar
hold on
imagesc(WL,'AlphaData',1-xform_isbrain)
hold on
contour(ROI,'k')
axis image off
title(strcat(contrastName{1},'(',ylabelName,')'));

subplot('position',[0.64,0.05,0.13,0.18])
imagesc(mean(downSampled(:,:,runInfo.stimStartTime+2:runInfo.stimEndTime,2),3),...
    [-cMax_temp(2)*1.2,cMax_temp(2)*1.2])
hold on
imagesc(WL,'AlphaData',1-xform_isbrain)
colorbar
hold on
contour(ROI,'k')
axis image off
title(strcat(contrastName{2},'(',ylabelName,')'));

if length(contrastName)>2
    subplot('position',[0.79,0.05,0.13,0.18])
    
    imagesc(mean(downSampled(:,:,runInfo.stimStartTime+2:runInfo.stimEndTime,3),3),...
        [-cMax_temp(3)*1.2,cMax_temp(3)*1.2
        ])
    
    colorbar
    hold on
    imagesc(WL,'AlphaData',1-xform_isbrain)
    hold on
    contour(ROI,'k')
    axis image off
    title(strcat(contrastName{3},'(',ylabelName,')'))
    colormap jet
end
