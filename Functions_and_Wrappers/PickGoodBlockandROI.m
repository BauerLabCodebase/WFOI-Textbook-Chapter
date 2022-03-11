function [fh,goodBlocks,ROI] = PickGoodBlockandROI(data,contrastName,range,runInfo,numBlocks,xform_isbrain)
% pres of each block during stim period
numRows = length(contrastName);
fh = figure('units','normalized','outerposition',[0 0 1 1]);
colormap jet
for ii = 1:numBlocks
    for jj = 1:length(contrastName)
        p = subplot(numRows,numBlocks,numBlocks*(jj-1)+ii);
        pre = squeeze(mean(data(:,:,runInfo.stimStartTime*runInfo.samplingRate+1:runInfo.stimEndTime*runInfo.samplingRate,ii,jj),3));
        imagesc(pre,[-range(jj),range(jj)]);
                    axis image
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
        if jj==1
        title(strcat('Pres',{' '},num2str(ii)))
        end
        if ii == 1
            ylabel(contrastName{jj})
        end
        if ii == numBlocks
            Pos = get(p,'Position');
            colorbar
            set(p,'Position',Pos)
        end
    end
end

% Pick good block
prompt = {'Which block is good?'};
title1 = 'Pick block';
dims = [1 35];
definput = {'[]'};
answer = inputdlg(prompt,title1,dims,definput);
goodBlocks_pre = str2num(answer{1});



% generate ROI
temp = squeeze(mean(data(:,:,runInfo.stimStartTime*runInfo.samplingRate+1:runInfo.stimEndTime*runInfo.samplingRate,goodBlocks_pre,end),3));
peakMap = squeeze(mean(temp,3));
f=figure;
imagesc(peakMap.*xform_isbrain)
colormap jet
colorbar
axis image off
title(contrastName(end));
[x1,y1] = ginput(1);
[x2,y2] = ginput(1);
[X,Y] = meshgrid(1:128,1:128);
radius = sqrt((x1-x2)^2+(y1-y2)^2);
ROI = sqrt((X-x1).^2+(Y-y1).^2)<radius;
max_ROI = prctile(peakMap(ROI),99);
temp = peakMap.*ROI;
ROI = temp>0.75*max_ROI;
hold on
contour(ROI,'k')
pause(0.1)
close(f)
% other way to find good blocks: mean>2*standard deviation
data = reshape(data,size(data,1)*size(data,2),[],numBlocks,size(data,5));
iri = reshape(ROI,1,[]);
goodBlocks_temp = [];
for ii = 1:numBlocks     
    timetrace = mean(data(iri,:,ii,end),1);
    meanValue = mean(timetrace(runInfo.stimStartTime*runInfo.samplingRate+1:runInfo.stimEndTime*runInfo.samplingRate));
    stdValue= std(timetrace(1:runInfo.stimStartTime*runInfo.samplingRate));
    if meanValue>2*stdValue
        goodBlocks_temp = [goodBlocks_temp,ii];
    end
end
goodBlocks = intersect(goodBlocks_pre,goodBlocks_temp);