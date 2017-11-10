% OPTBINS computes the optimal number of bins for a given one-dimensional 
% data set. This optimization is based on the posterior probability for 
% the number of bins
%
% Usage:
%           optM = OPTBINS(data,maxM);
%           [optM, optlogp] = OPTBINS(data,maxM);
%
% Where:
%           data is a (1,N) vector of data points
%           maxM is the maximum number of bins to consider
%
%
% This algorithm uses a brute-force search trying every possible bin number  
% in the given range.  This can of course be improved.
% Generalizations to multidimensional data sets is straightforward.
%
% Created by Kevin H. Knuth on 17 April 2003
% Modified by Kevin H. Knuth on 17 March 2006: added optlogp
% Modified by Kevin H. Knuth on 16 October 2006: added isROUNDED


function [optM, optlogp] = OPTBINS(data, maxM)

if numel(size(data))>2 || (size(data,1)>1)
    error('data dimensions must be (1,N)');
end

N = size(data,2);

% Loop through the different numbers of bins
% and compute the posterior probability for each.

logp = zeros(1,maxM);
data(isnan(data) | isinf(data)) = [];
for M = 1:maxM
    
        % Bin the data (equal width bins here)
        % to use mex file, comment out the line below and uncomment the
        % following line
        counts = histogram(data, M);        %   compute likelihood   
%       [counts, binwidth, centers] = histomex(data, M);
        
    logp(M) = N*log(M) + gammaln(M/2) - gammaln(N+M/2) - ...
                M*gammaln(1/2) + sum(gammaln(counts+0.5));
    
end

[maximum, optM] = max(logp);
optlogp = logp(optM);

% Test to see if data are excessively rounded or truncated
rounded = isROUNDED(data, optM);
if rounded
    disp('The data are excessively rounded.');
    warning('Add a small number to each datum vector to obtain a reasonable result.');
end

return