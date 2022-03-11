function [isbrain,WL] = getMask_noLandmark(raw,DarkFrameInd,InvalidInd,rgbOrder)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
import mouse.*

badDataInd = unique([DarkFrameInd InvalidInd]);
realDataStart = max(badDataInd) + 1;
raw = raw(:,:,:,realDataStart:end);
raw = raw(:,:,:,size(raw,4));
raw = single(raw);
WL = getWL(raw,rgbOrder);
isbrain = getIsBrain(WL);

end


function wl = getWL(data,rgbOrder)
%getWL Gets white light image from data
%   data = 3D raw data
%   rgbOrder = 1x3 vector that holds indices of red LED, green LED, and
%   blue LED

frameData = double(data);
frameData = frameData./repmat(max(max(frameData)),size(frameData,1),size(frameData,2));

ledNum = size(data,numel(size(data)));

availableLEDs = 1:ledNum;
nonNanInd = ~isnan(rgbOrder);
availableLEDs(nonNanInd) = [];

% replace nan values in rgbOrder
nanInd = 0;
for rgbInd = 1:numel(rgbOrder)
    if isnan(rgbOrder(rgbInd))
        nanInd = nanInd + 1;
        if nanInd > max(availableLEDs)
            nanInd = 1;
        end
        rgbOrder(rgbInd) = availableLEDs(nanInd); % replace nan with first non nan index
    end
end

wl = squeeze(frameData(:,:,rgbOrder,:));

end



function isbrain = getIsBrain(WL)
% loads mask file in and outputs white light image while getting landmarks
% for affine transform
%   Input:
%       WL = white light image (n x n x 3)
%   Output:
%       isbrain = logical array (n x n)

polyFig = figure;
disp('Click the boundary between brain and non-brain to create a polygon describing area containing the brain.');
mask = roipoly(WL);
close(polyFig);

if ~any(any(mask))
    error('No mask given from the gui.');
else
    isbrain = single(uint8(mask));
end

isbrain = isbrain > 0;

end