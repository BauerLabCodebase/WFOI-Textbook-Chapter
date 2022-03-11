function [fh] = generateBlockMap(data,contrastName,runInfo,numBlocks,peakRange,xform_isbrain)
% pres of each block during stim period
parampath=what('bauerparams');
load(strcat(parampath.path,'\noVasculatureMask.mat')) %anmol - point to local novascmask
mask_new = logical(mask_new);
numRows = length(contrastName);
fh = figure('units','normalized','outerposition',[0 0 1 1]);
colormap jet
for jj = 1:length(contrastName)
     peakMap = squeeze(mean(data(:,:,:,:,jj),4));
     peakMap = squeeze(mean(peakMap(:,:,(round(peakRange{jj}(1)*runInfo.samplingRate)):(round(peakRange{jj}(end)*runInfo.samplingRate))),3));
     maxVal = prctile(abs(peakMap(mask_new)),90,'all')*2.5;
     maxVal = round(maxVal,2,'significant');
    for ii = 1:numBlocks+1
        p = subplot(numRows,numBlocks+1,(numBlocks+1)*(jj-1)+ii);
        if ii == numBlocks+1
            imagesc(peakMap,'AlphaData',xform_isbrain,[-maxVal,maxVal])
            axis image
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            title('Averaged')
            Pos = get(p,'Position');
            cb = colorbar;
            caxis([-maxVal maxVal])
            set(p,'Position',Pos)
            set(cb,'YTick',[-maxVal,0,maxVal]);

            if sum(contains({'HbO','HbR','HbT'},contrastName{jj}))>0
            set(get(cb,'label'),'string','Hb(\DeltaM)');
            elseif sum(contains({'Calcium','FAD'},contrastName{jj}))>0
                set(get(cb,'label'),'string','Fluorescence(\DeltaF/F)');
            end
        else
            pre = squeeze(mean(data(:,:,(round(peakRange{jj}(1)*runInfo.samplingRate)):(round(peakRange{jj}(end)*runInfo.samplingRate)),ii,jj),3));
            imagesc(pre,'AlphaData',xform_isbrain,[-maxVal,maxVal]);
            axis image
            set(gca, 'XTick', []);
            set(gca, 'YTick', []);
            if jj==1
                title(strcat('Pres',{' '},num2str(ii)))
            end
            if ii == 1
                ylabel(contrastName{jj})
            end
        end
    end
end
