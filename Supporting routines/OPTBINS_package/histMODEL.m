%histMODEL  Creates a multi-dimensional histmodel.
%   MODEL = histMODEL(Y,M) bins the elements of the array Y into a set of containers
%                 in a multi-dimensional space, and returns the number of elements 
%                 in each container.
%
%   Y is a D x N array, where Y(i,:) is the ith datum point.
%   M is a 1 x D array, where M(j) is the number of bins in the jth dimension
%
%   model is a cell array representing a histogram model.
%       It contains:
%           model{1} = D is the dimensionality of the space
%           model{2} = N is the total number of counts
%           model{3} = M is a 1 x D array of the number of bins along each dimension d
%           model{4} = ps is a multidimensional array of mean bin probabilities with dimensions dictated by M
%           model{5} = vs is a multidimensional array of bin probability variances with dimensions dictated by M
%           model{6} = binwidth in the dth dimension is binwidth(d)
%           model{7} = centers;
%
%   Created by Kevin H. Knuth on 15 August 2003
%   Notes modified by KHK 4 April 2005

function model = histMODEL(Y,M)

if (nargin ~= 2)
    error('Requires two input arguments.')
end

[D, N] = size(Y);

if (length(M) ~= D)
    error('The dimensions of the data and the bin number array disagree.');
end

% find the extrema
ymin = zeros(1,D);
ymax = zeros(1,D);
for d = 1:D
    ymax(d) = max(Y(d,:));
    ymin(d) = min(Y(d,:));
end

% define the bin edges
% the extreme bins go off to infinity
edges = zeros(D,max(M));
binwidth = (ymax-ymin)./M;
for d = 1:D
    m = M(d)-1;
    edges(d,1:M(d)-1) = binwidth(d)*(1:m)+ymin(d);
    edges(d,M(d)) = Inf;  % the rightmost bin edge is at negative Infinty
end




% obtain the counts %%%%%%%%%%%%%%%%%%%%%%%%%
if D == 1
    counts = zeros(1,M);
else
    counts = zeros(M);
end

% setup indexing factors
% See: Programming and Data Types: M-File Programming: Advanced Indexing
factors = ones(1,D);
for d = 1:D-1
    factors(d+1:D) = factors(d+1:D)*M(d);
end

for n = 1:N

    % get the multidimensional bin coordinates of each data point
    bin = zeros(1,D);
    for d = 1:D
        % determine which bin the count belongs in for each dimension (these are its coords in the count array)
        % I go through these from right to left to ensure that this algorithm is consistent with hist.m
        % I can put in a binary search algorithm here which will be more efficient
        for m = 1:M(d)  %loop through all the bins in this dimension
            if (Y(d,n) <= edges(d,m))
                bin(d) = m;
                break;
            end
        end
    end
    
    % now add a count to that bin
    index = factors * (bin-1)' + 1;
    counts(index) = counts(index)+1;
end

% define the bin centers
centers = cell(D);
for d = 1:D
    centers(d) = {(binwidth(d)*((1:M(d))-1)+ymin(d))+binwidth(d)/2};
end

% make ps
vol = prod(binwidth);
B = prod(M);  % get the total number of bins
ps = ((counts+0.5)/(N+B/2))/vol; 
% must scale by volume of a bin to assure that probs sum to zero
% remember that this is a model of a pdf---not a histogram

% make vs
vs = ps .* (N-counts+(B-1)/2)/((N+B/2)*(N+B/2+1));

model{1} = D;
model{2} = N;
model{3} = M;
model{4} = ps;
model{5} = vs;
model{6} = binwidth;
model{7} = centers;

return;
