% mOPTBINS computes the optimal number of bins for a given 
% multi-dimensional data set. This optimization is based on the posterior 
% probability for the number of bins
%
% Usage:
%           moptM = mOPTBINS(data, maxM);
%
% Where:
%           data is a (D,N) vector of data points
%           maxM is an (1,D) array of the maximum number of bins to 
%                   consider for each dimension
%
%
% This algorithm uses a brute-force search trying every possible bin number  
% in the given range.  This can of course be improved.
% Generalizations to multidimensional data sets is straightforward.
%
% Created by Kevin H. Knuth on 17 April 2003
% Modified for use on multi-dimensional data on 18 Aug 2003
% Modified by Knuth on 12 Aug 2004 - commented on storing logP
% Modified by Knuth on 18 Oct 2006 - trimmed down and cleaned up


function moptM = mOPTBINS(data, maxM)

[D, N] = size(data);
if length(maxM) < D
    error('Dimensions of the bin limit does not agree with the data');
end


if (D == 1)     % if it is a 1D problem, use optBINS()
    moptM = optBINS(data, maxM);
else
    
    % BRUTE FORCE: Loop through the different numbers of bins
    % and compute the posterior probability for each.
    
    % setup indexing factors
    % Programming and Data Types: M-File Programming: Advanced Indexing
    factors = ones(1,D);
    for d = 1:D-1
        factors(d+1:D) = factors(d+1:D)*maxM(d);
    end
    
    
    % preliminaries
    M = ones(1,D);
    maxLOGP = -Inf;
    moptM = M;
    
    keepgoing = 1;
    while (keepgoing)   % loop until the bin assignments are all maxed out
        
        % Bin the data (equal width bins here)
        % to use mex file, comment out the line below and uncomment the
        % following line
        counts = histogram(data, M);        %   compute likelihood   
%       [counts, binwidth, centers] = histomex(data, M);
        
        MM = prod(M);
        % compute log prob
        logprob = N*log(MM) + gammaln(MM/2) - gammaln(N+MM/2) + ...
            - MM*gammaln(1/2) + sum(reshape(gammaln(counts+0.5),1,MM));
        
        % update optimal value
        if logprob > maxLOGP
            moptM = M;
            maxLOGP = logprob;
        end
        
        if (sum(M < maxM))  % keep going, we are not done
            
            % increment the bin sizes
            d = 1;  % start at the first dimension of the space
            carryflag = 1;  % signal need to increment the next higher dim
            while (d <= D) && (carryflag == 1)
                M(d) = M(d)+1;
                if (M(d) > maxM(d))     % if dim d has exceeded the limit
                    M(d) = 1;           % start it over
                    d = d+1;            % inc bins in the next higher dim
                else
                    carryflag = 0;      % no carry necessary -> stop
                end
            end
            
        else    % we are all done
            keepgoing = 0;
        end 
    end
end


% Test to see if data are excessively rounded or truncated
rounded = isROUNDED(data, moptM);
if rounded
    disp('The data are excessively rounded.');
    warning('Add a small number to each datum vector to obtain a reasonable result.');
end

return