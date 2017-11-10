function y = TukeyWindow(range, breakLeft, fractionTaper, theta, symmetricWindow)
% Determine full or half window
if ~exist('symmetricWindow', 'var')
    symmetricWindow = 0;
end

% Find the end points of Tukey window in true coordinate
r = fractionTaper;
A = range(1);
B = range(2);
a = r/4;
b = 1 - a;
endPointNorm = [0 1];
endPointTrue = (B-A) * (endPointNorm - a) / (b-a) + A;

% Create the Tukey window
y = zeros(size(theta));
nPoints = sum(theta >= endPointTrue(1) & theta <= endPointTrue(2));
y(theta >= endPointTrue(1) & theta <= endPointTrue(2)) = tukeywin(nPoints, r);

% Create the sharp edge on one side of the window
if symmetricWindow == 0
    switch breakLeft
        case 1
            midPoint = (A - B)/2;
            y(theta > B) = 0;
            y(theta >= midPoint & theta <= B) = max(y);
        case 0
            midPoint = (B - A)/2;        
            y(theta < A) = 0;
            y(theta >= A & theta <= midPoint) = max(y);
        case 'NA'
    end
end

% Normalize to make a probability distribution
scalingFactor = trapz(theta, y);
y = y / scalingFactor;
end