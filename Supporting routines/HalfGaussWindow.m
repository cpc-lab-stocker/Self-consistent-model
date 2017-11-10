function y = HalfGaussWindow(range, breakLeft, stdGauss, theta)
% Initialize the window 
y = zeros(size(theta));

% Create the window
switch breakLeft
    case 1
        indexGauss = theta <= range(1); 
        y(indexGauss) = normpdf(theta(indexGauss), range(1), stdGauss);
        y(theta >=range(1) & theta <= range(2)) = max(y(indexGauss));
    case 0
        indexGauss = theta >= range(2); 
        y(indexGauss) = normpdf(theta(indexGauss), range(2), stdGauss);
        y(theta >= range(1) & theta <= range(2)) = max(y(indexGauss));        
    case 'NA'
end
    
% Normalize to make a probability distribution
scalingFactor = trapz(theta, y);
y = y / scalingFactor;
end