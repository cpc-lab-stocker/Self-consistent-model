function alpha = alphaCompute(priorRange, fractionTaper)
priorRange = [-priorRange zeros(length(priorRange), 1)];
theta = linspace(-50, 50, 100000);
alpha = NaN(length(priorRange), 1);
for ii = 1 : length(priorRange)
    y = TukeyWindow(priorRange(ii,:), 1, fractionTaper(ii), theta);
    y(theta>=0) = [];
    indexPeak = find(y==max(y),1);
    indexZero = find(y==0,1,'last');
    alpha(ii) = theta(indexPeak) -  theta(indexZero);
end
