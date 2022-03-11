function normData = dataNormalizer(data)
% dataNormalizer(data) normalizes each row vectors of the input data
% matrix. 

%zscore each row

normData = zeros(size(data));
ceteredData = data - mean(data,2);
for i = 1: size(data,1)
    normData(i,:) = ceteredData(i,:)./norm(ceteredData(i,:),2);
end