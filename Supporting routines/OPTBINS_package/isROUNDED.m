% isROUNDED
% GNU General Public License Software Copyright: Kevin H. Knuth 2006
%
% isROUNDED determines whether the data are excessively rounded or
% truncated
%
% Usage:
%           flag = isROUNDED(data, M)
%           
% Where:
%           M  - is the optimal number of bins
%           flag - 1 for yes, 0 for no
%
% Created by:   Kevin Knuth
%           	16 Oct 2006
% Modified by: KHK 27 Feb 2007
%               Bug Fix regarding instance count
%

function result = isROUNDED(data, m)

[D, N] = size(data);

% find out how many unique data vectors there are
S = unique(data','rows')';
P = length(S);
% count the numbers of instances
n = zeros(1,P);
for i = 1:P
    summ = zeros(1,N);
    for j = 1:D
        summ = summ + ismember(data(j,:),S(j,i));
    end
    n(i)=numel(find(summ==D));
end

% compute asymptotic logP value
logs = 0;
for i = 1:P
    for s = n(i):2*n(i)-1
        logs = logs + log(s);
    end
end
asymptote = (P-N)*log(2) + logs;


% find optimal logP
counts = histogram(data, m);    %   compute likelihood

% compute log prob
MM = prod(m);
optlogP = N*log(MM) + gammaln(MM/2) - gammaln(N+MM/2) + ...
    - MM*gammaln(1/2) + sum(reshape(gammaln(counts+0.5),1,MM));


result = 1; %assume the data is bad then determine if it is good
if optlogP >= asymptote %log p is a value
    result = 0;
end

return