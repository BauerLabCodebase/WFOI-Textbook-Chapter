function [lagTime,lagAmp,covResult] = dotLag(data1,data2,edgeLen,validRange,varargin)
%dotLag Finds lag within each pixel of two data sets
%   Inputs:
%       data1 = pix dim1 x pix dim2 x time
%       data2 = pix dim1 x pix dim2 x time

%       edgeLen = data points from peak to consider for quadratic fit (if 3, 7
%       points around are considered)

%       validRange = How many data points to consider from zero
%       lag. For example, -3:3 would mean looking at 7 points around zero lag.

%       corrThr (optional) = the threshold at which correlation is considered
%       significant. (default = 0.3)
%       positiveSignOnly (optional) = whether to ignore negative xcorr values (default = true)
%       quadFitUse (optional) = whether to use quadratic fitting (default =
%       true)
%   Outputs:
%       lagTime = pix dim1 x pix dim2

%       lagAmp = pix dim1 x pix dim2

%       covResult = (dim1*dim2) x time frames

if numel(varargin) > 0
    corrThr = varargin{1};
else
    corrThr = 0.3;
end
if numel(varargin) > 1
    positiveSignOnly = varargin{2};
else
    positiveSignOnly = true;
end
if numel(varargin) > 2
    quadFitUse = varargin{3};
else
    quadFitUse = true;
end

% validInd = [];
covResult = nan(size(data1,1)*size(data1,2),numel(validRange));

% initialize outputs
lagTime = nan(size(data1,1),size(data1,2)); %size of the image
lagAmp = nan(size(data1,1),size(data1,2)); %size of the image



for pixRow = 1:size(data1,1) % for each pixel
    for pixCol = 1:size(data1,2)
        lagData1 = squeeze(data1(pixRow,pixCol,:));
        lagData2 = squeeze(data2(pixRow,pixCol,:));
        [lagTime(pixRow,pixCol),lagAmp(pixRow,pixCol),covResultPix]...
            = findLag(lagData1,lagData2,quadFitUse,positiveSignOnly,validRange,...
            edgeLen,corrThr); % finding lag; data 1 lags data 2 by how much
        
        covResult(pixRow+(pixCol-1)*size(data1,1),:) = covResultPix;
%         if ~isnan(lagTime(pixRow,pixCol))
%             covResult(pixRow+(pixCol-1)*size(data1,1),:) = covResultPix;
%         end
    end
end
end

