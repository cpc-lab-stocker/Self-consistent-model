function [dataNew, indexOutlier, valueOutlier] = removeOutlierCarling(dataIn)
%%%%% Implement Carling (2000) method shown in Wilcox (2005) book p.100 %%%%%
%%% Input: dataIn - a nxm matrix with the data to be processed along the
%%%                 column 
%%%        dataNew - processed data with outliers replaced by NaN
%%%        indexOutlier - the indices of outliers
%%%        valueOutlier - the values of outliers
%
dimDataIn = size(dataIn);

% Check for empty data
if sum(dimDataIn) == 0
    dataNew = dataIn;
    indexOutlier = [];
    valueOutlier = [];
    return
end

% Start the removal algorithm
if dimDataIn(1) == 1
    dataIn = dataIn';
else
   nLoops = dimDataIn(2);
end
dataNew = NaN(dimDataIn);
indexOutlier = cell(1, dimDataIn(2));
valueOutlier = cell(1, dimDataIn(2));
for ii = 1 : nLoops
    tempData = dataIn(:,ii);
    nDataSingular = sum(isnan(tempData) | isinf(tempData));
    tempData(isnan(tempData)) = [];
    nSample = length(tempData);    
    tempData = sort(tempData, 'ascend');

    % Calculate the 2 quantiles
    j = floor(nSample/4+5/12);
    if j < 1 || j > nSample
        dataNew(:, ii) = dataIn(:, ii);
        indexOutlier = [];
        valueOutlier = [];
    else        
        h = nSample/4 + 5/12 - j;
        q1 = (1-h) * tempData(j) + h * tempData(j+1);
        k = nSample - j + 1;
        q2 = (1-h) * tempData(k) + h * tempData(k-1);

        % Calculate the removal bounds
        k = (17.63*nSample - 23.64) / (7.74*nSample - 3.71);
        lowerBound = median(tempData) - k * (q2-q1);
        upperBound = median(tempData) + k * (q2-q1);

        % Detect and remove the outliers
        indexOutlier{ii} = (tempData > upperBound) | (tempData < lowerBound);
        valueOutlier{ii} = tempData(indexOutlier{ii});
        tempData(indexOutlier{ii}) = NaN;
        dataNew(:,ii) = [tempData; NaN(nDataSingular,1)];
    end
end
end