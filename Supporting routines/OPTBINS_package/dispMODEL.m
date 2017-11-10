% dispMODEL plots a bar graph with mean bin heights and error bars from a
% 1-dimensional histogram model
%
%   h = dispMODEL(histmodel, style)
%
%   INPUT:
%   histmodel is a 1-dimensional histogram model
%   style is either 'errorbars' or 'thick'
%   
%   OUTPUT:
%   h is the handle to the figure
%
%   Created by Kevin H. Knuth on 17 February 2005
%   Modified to dispMODEL on 3 May 2005

function h = dispMODEL(histmodel, varargin)

D = histmodel{1};

if D>2
    error('this will only plot 1-D histogram models');
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
            error('dispMODEL:invalidargument','Argument not recognized');
    end
    k = k+1;
end


% N = histmodel{2};
M = histmodel{3};
p = histmodel{4};
s = sqrt(histmodel{5});
binwidth = histmodel{6};
c = histmodel{7}{1};

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