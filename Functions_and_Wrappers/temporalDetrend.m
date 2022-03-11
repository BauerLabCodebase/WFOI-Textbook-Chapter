function result = temporalDetrend(inputData,isbrain,varargin)
%temporalDetrend Temporally detrends by fiting to a polynomial.
% input:
%   inputData = nD matrix. Last dim is time.
%   usePoly (optional) = Boolean determining whether to use polyfit option.
%   Default is false.

% result = inputData;

if numel(varargin) > 0
    usePoly = varargin{1};
else
    usePoly = false;
end


raw_isbrain_tmp=zeros(size(inputData,1),size(inputData,2),size(inputData,3),size(inputData,4));
for i=1:size(inputData,3)
raw_isbrain_tmp(:,:,i,:)=squeeze(inputData(:,:,i,:)).*repmat(isbrain,1,1,size(inputData,4));
end

dataSize = size(inputData);
data = reshape(raw_isbrain_tmp,[],dataSize(end));
clear raw_isbrain
isbrain_tmp=repmat(isbrain,1,1,size(inputData,3),size(inputData,4));
isbrain_tmp=reshape(isbrain_tmp,[],dataSize(end));

for c=1:size(data,1) %each pixel
    if isbrain_tmp(c,1)
    data(c,:)=data(c,:)-polyval(polyfit(1:dataSize(end), data(c,:), 5), 1:dataSize(end))+mean(data(c,:));
    warning('off')
    end
end

clear raw_isbrain raw_isbrain_tmp

result = reshape(data,dataSize);

end