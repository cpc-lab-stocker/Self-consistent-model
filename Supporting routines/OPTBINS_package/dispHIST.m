% dispHIST plots a bar graph with mean bin heights and error bars
%
%   h = dispHIST(data, M, style)
%
%   INPUT:
%   data is a 1 x N array of values
%   M is the number of bins
%   style is either 'errorbars' or 'thick'
%   
%   OUTPUT:
%   h is the handle to the figure
%
%   Created by Kevin H. Knuth on 17 February 2005

function h = dispHIST(data, M, varargin)

if numel(size(data))>2 || (size(data,1)>1)
    error('data dimensions must be (1,N)');
end

% handle extra conditions
THICK = 1;
ERRORBARS = 2;
style = ERRORBARS;
k = 1;
while k <= length(varargin)
    switch lower(varargin{k})
        case 'thick'
            style = THICK;
        case 'errorbars'
            style = ERRORBARS;
        otherwise
            error('dispHIST:invalidargument','Argument not recognized');
    end
    k = k+1;
end


N = size(data,2);

[n binwidth centers] = histomex(data,M);            % determine the bin counts
c = centers{1};
p = ((n+0.5)/(N+M/2)/binwidth);                                % determine the mean bin probabilities
s = (sqrt(p .* (N-n+(M-1)/2)/((N+M/2+1)*(N+M/2)))/binwidth);   % determine the std devs

h = figure;
if M==1
    % Fake a bar graph
    x = [c-binwidth/2, c+binwidth/2];
    y = [p,p];
    area(x,y,0);
    axis([c-2*binwidth/3 c+2*binwidth/3 0 7*p/6]);
    if (style == THICK)
        colormap copper;
    end
else
    % Make a bar graph
    if (style == THICK)
        b = zeros(M,3);
        b(:,1) = (p-s)';
        b(:,2) = s';
        b(:,3) = s';
        bar(c, b, 1, 'stack');
        colormap copper;
    else
        bar(c, p, 1);

        for m = 1:M
            line([c(m); c(m)], [p(m); p(m)-s(m)], 'Color', 'w');
            line([c(m)- binwidth/6; c(m)+ binwidth/6], [p(m)-s(m); p(m)-s(m)], 'Color', 'w');
            line([c(m); c(m)], [p(m); p(m)+s(m)], 'Color', 'k');
            line([c(m)- binwidth/6; c(m)+ binwidth/6], [p(m)+s(m); p(m)+s(m)], 'Color', 'k');
        end
        colormap('default')
    end
end

return;