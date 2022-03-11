function [lagTime,lagAmp,covResult,sign] = findLag(data1,data2,varargin)
%findLag Calculates cross covariance and finds the peak to obtain lag.
%   There is also the option to use quadratic fit to get the lag time.
%       data1 = time series, 1D, already zero-meaned
%       data2 = time series, 1D, already zero-meaned
%       quadFitUse (optional) = whether to use quadratic fitting (default =
%       true)
%       positiveSignOnly (optional) = whether to ignore negative xcorr values (default = true)
%       validRange (optional) = How many data points to consider from zero
%       lag. For example, -3:3 would mean looking at 7 points around zero lag.
%       Empty implies default. (default = all values considered)
%       edgeLen (optional) = data points to actually consider in quadratic fit.
%       If zero, then quadratic fit isn't used.
%       (default = 3)
%       corrThr (optional) = the threshold at which correlation is considered
%       significant. (default = 0.3)
%   Output:
%       lagTime = if positive, data1 lags behind data2.
%       lagAmp = how strong the dependence is.
%       covResult = cross covariance result
%       sign = whether the lag is from positive correlation or negative
%       correlation. (1 = correlation, 0 = anti-correlation)
%
% (c) 2018 Washington University in St. Louis
% All Right Reserved
%
% Licensed under the Apache License, Version 2.0 (the "License");
% You may not use this file except in compliance with the License.
% A copy of the License is available is with the distribution and can also
% be found at:
%
% http://www.apache.org/licenses/LICENSE-2.0
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE, OR THAT THE USE OF THE SOFTWARD WILL NOT INFRINGE ANY PATENT
% COPYRIGHT, TRADEMARK, OR OTHER PROPRIETARY RIGHTS ARE DISCLAIMED. IN NO
% EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
% OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
% ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

if size(data1,1) == 1
    data1 = data1'; % make it time x 1
end
if size(data2,1) == 1
    data2 = data2'; % make it time x 1
end

if numel(varargin) > 0
    quadFit = varargin{1};
else
    quadFit = true;
end
if numel(varargin) > 1
    positiveSignOnly = varargin{2};
else
    positiveSignOnly = true;
end
if numel(varargin) > 2
    validRange = varargin{3};
else
    validRange = -size(data1,1):size(data1,1);
end
if numel(varargin) > 3
    edgeLen = varargin{4};
else
    edgeLen = 3;
end
if numel(varargin) > 4
    corrThr = varargin{5};
else
    corrThr = 0.3;
end

if edgeLen == 0 % if the number of data in edge is zero
    quadFit = false; % just don't quadratic fit
end

% unifying input data format to column vector
if ~isempty(data2) && size(data2,1) == 1
    data2 = data2';
end

maxValidRange = max(abs(validRange));
[covResult,lags] = xcorr(data1,data2,maxValidRange,'coeff'); % how much data1 lags behind data2
validInd = (min(validRange)+maxValidRange+1):(maxValidRange+1+max(validRange));

% vectorize cross correlation
covResult = covResult(validInd); covResult = covResult(:);
lags = lags(validInd); lags = lags(:);

[peakInd,sign] = findPeak(covResult,positiveSignOnly);
[lagTime, lagAmp] = getLagTimeAndAmp(covResult,lags,quadFit,peakInd,edgeLen,corrThr);
lagTime = real(lagTime); lagAmp = real(lagAmp);
end

function [peakInd,sign] = findPeak(covResult,positiveSignOnly)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% find peaks
posCov = covResult;
posCov(posCov < 0) = 0;
y = diff(posCov);
posExtremaLoc = find(y(1:end-1) > 0 & y(2:end) < 0) + 1;
posExtrema = posCov(posExtremaLoc);

if positiveSignOnly
    negExtrema = [];
    negExtremaLoc = [];
else
    negCov = -covResult;
    negCov(negCov < 0) = 0;
    y = diff(negCov);
    negExtremaLoc = find(y(1:end-1) > 0 & y(2:end) < 0) + 1;
    negExtrema = negCov(negExtremaLoc);
%     [negExtrema, negExtremaLoc] = findpeaks(negCov);
end

% find the highest peaks in both positive and negative
posExtremaInd = find(posExtrema == max(posExtrema),1);
negExtremaInd = find(negExtrema == max(negExtrema),1);
posExtrema = posExtrema(posExtremaInd);
posExtremaLoc = posExtremaLoc(posExtremaInd);
negExtrema = negExtrema(negExtremaInd);
negExtremaLoc = negExtremaLoc(negExtremaInd);

% choose positive or negative peak
if numel(negExtrema)>0 && numel(posExtrema)>0
    if posExtrema >= negExtrema
        peakInd = posExtremaLoc;
        sign = 1;
    else
        peakInd = negExtremaLoc;
        sign = 0;
    end
elseif numel(posExtrema)>0
    peakInd = posExtremaLoc;
    sign = 1;
elseif numel(negExtrema)>0
    peakInd = negExtremaLoc;
    sign = 0;
else
    peakInd = NaN;
    sign = NaN;
end
end

function [lagTime, lagAmp] = getLagTimeAndAmp(covResult,lags,quadFitUse,peakInd,edgeLen,corrThr)
if isnan(peakInd)
    lagTime = NaN;
    lagAmp = NaN;
else
    if quadFitUse
        % find the length of edges to look at for polyfit (edge cases
        % considered)
        leftEdge = edgeLen;
        if peakInd < edgeLen+1
            leftEdge = peakInd - 1;
        end
        rightEdge = edgeLen;
        if numel(covResult) - peakInd < edgeLen+1
            rightEdge = numel(covResult) - peakInd;
        end
        
        % quadratic fit
        lagInd = lags(peakInd-leftEdge:peakInd+rightEdge);
        lagInd = lagInd - lagInd(1);
        [~,~,peakX,peakY] = quadFit(lagInd,covResult(peakInd-leftEdge:peakInd+rightEdge));
        peakX = peakX + lags(peakInd-leftEdge);
    else
        peakX = lags(peakInd);
        peakY = covResult(peakInd);
    end
    
    lagTime = peakX;
    lagAmp = peakY;
end

% % if lag amplitude isn't significant enough, the lag does not exist.
% if abs(lagAmp) < corrThr
%     lagTime = NaN;
% end
end