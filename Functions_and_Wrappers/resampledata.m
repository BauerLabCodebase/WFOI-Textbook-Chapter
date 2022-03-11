function [dataOut, tOut] = resampledata(dataIn,tIn,tOut)
%resampledata Resamples data. This resampling algorithm first low-pass
%filters data, then interpolates the data over output time.
% Inputs:
%   dataIn = nD matrix. Last dimension is time.
%   tIn = time points of input. This should be the same length as the
%   length in the time dimension of dataIn
%   tOut = time points of output. If this is 1 value, then it is considered
%   fOut (output frequency).

dataIn(isinf(dataIn)) = 0;
dataIn(isnan(dataIn)) = 0;

sizeData = size(dataIn);
dataIn = reshape(dataIn,[],sizeData(end));

if numel(tOut) == 1
    fOut = tOut;
    tOut = min(tIn):1/fOut:max(tIn);
else
    fOut = 1/median(diff(tOut));
end
fIn = 1/median(diff(tIn));

%% band pass the data

typeFcn = @single;
if strcmp(string(class(dataIn)),"double")
    typeFcn = @double;
end

if size(dataIn,numel(size(dataIn))) > 15
    if (fOut < fIn/2) % if butterworth filtering is actually possible
        omegaNy = fOut*(2/fIn);
        [b,a] = butter(5,omegaNy);
        for i = 1:size(dataIn,1)
            dataIn(i,:) = typeFcn(filtfilt(b,a,double(dataIn(i,:))));
        end
    else
        disp('during resampling process filtering was not done since the output frequency is more than half of input frequency.');
    end
else
    disp('during resampling process filtering was not done since the number of time frames given is too small.');
end

dataOut = typeFcn(nan(size(dataIn,1),numel(tOut)));
%% interpolate over sinc kernel

for pix = 1:size(dataIn,1)
    dataOut(pix,:) = interp1(tIn,dataIn(pix,:),tOut,'linear');
end

dataOut = reshape(dataOut,[sizeData(1:end-1) numel(tOut)]);
end

