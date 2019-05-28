% conditional observer model with memory recall and motor noise
% for orientation stimulus (infinite space approximation) 
% astocker - lluu
% 06.2017

function selfConsistentModel
flagDecisionGiven = 0; % 1: stimulus category (cw/ccw) are given (Exp. 3)
                       % 0: observer makes categorical decision (Exp. 1 and 2)

dstep = 0.1;
paramsAll = [3.66      5.50      9.24           0.0000       39.18       5.35     0.32      3.60];
lapseRate = paramsAll(4);

% stimulus orientation
thetaStim = -21:0.1:21; % 
thetaStim = round(thetaStim, -log10(dstep));

% sensory noise
stdSensory = paramsAll(1:3);

% memory recall noise
stdMemory = paramsAll(6);

% motor noise;
stdMotor = paramsAll(8);

% priors
pC = [0.5, 0.5]'; % [cw ccw]
pthcw = paramsAll(5);
pthccw = -paramsAll(5);
smoothFactor = paramsAll(7);

%% LOOP - noise levels
rangeth = [-42 42];
th = rangeth(1):dstep:rangeth(2);
th = round(th, -log10(dstep));
nth = length(th);

pthGC = zeros(2,nth);

pthGC(1,:) = TukeyWindow([0 pthcw], 0, smoothFactor, th);
pthGC(2,:) = TukeyWindow([pthccw 0], 1, smoothFactor, th);

figure;
fontSize = 15;
for kk=1:length(stdSensory)  
    rangeM = [min(thetaStim)-5*stdSensory(kk) max(thetaStim)+5*stdSensory(kk)];
    if rangeM(2) < rangeth(2)
        rangeM = rangeth;
    end
    nm = 1000;
    m = linspace(rangeM(1), rangeM(2), nm);

    nmm = 1200;
    rangeMM = [min(rangeM)-6*stdMemory max(rangeM)+6*stdMemory];
    if rangeMM(2) < rangeth(2)
        rangeMM = rangeth;
    end        
    mm = linspace(rangeMM(1), rangeMM(2), nmm);

    M = repmat(m',1,nth);
    MM_m = repmat(mm',1,nm);
    MM_th = repmat(mm',1,nth); 
    THm = repmat(th, nm, 1); 
    THmm = repmat(th, nmm, 1); 
    
    %% Generative (forward)
    % orientation noise p(m|th)
    pmGth = exp(-((M-THm).^2)./(2*stdSensory(kk)^2));
    pmGth = pmGth./(repmat(sum(pmGth,1),nm,1)); 

    %% Inference
    % 1: categorical judgment
    if ~flagDecisionGiven
        PCGm = (pthGC * pmGth') .* repmat(pC,1,nm);
        % fix the issue when sensory noise is too low
        indFirstNonZero = find(PCGm(2,:), 1);
        PCGm(2, 1: indFirstNonZero-1) = PCGm(2, indFirstNonZero);
        indLastNonZero = find(PCGm(1,:), 1, 'last');
        PCGm(1, indLastNonZero+1:end) = PCGm(1, indLastNonZero);
        PCGm = PCGm./(repmat(sum(PCGm,1),2,1));
        % max posterior decision
        PChGm = round(PCGm);
        % marginalization
        PChGtheta = PChGm * pmGth(:, ismember(th, thetaStim));
        PChGtheta_lapse = lapseRate + (1 - 2*lapseRate) * PChGtheta;
        PChGtheta_lapse = PChGtheta_lapse ./ repmat(sum(PChGtheta_lapse, 1), 2, 1);
    else
        PChGtheta_lapse = NaN(2, length(thetaStim));
        PChGtheta_lapse(1, thetaStim > 0) = 1;
        PChGtheta_lapse(1, thetaStim < 0) = 0;
        PChGtheta_lapse(1, thetaStim == 0) = 0.5;
        PChGtheta_lapse(2, :) = 1 - PChGtheta_lapse(1, :);
    end        
    % 2: estimation
    pmmGth = exp(-((MM_th-THmm).^2)./(2*(stdSensory(kk)^2 + stdMemory^2))); % p(mm|th) = N(th, sm^2 + smm^2)
    pmmGth = pmmGth./(repmat(sum(pmmGth,1),nmm,1)); 
    pthGmmChcw = (pmmGth.*repmat(pthGC(1,:),nmm,1))';
    pthGmmChcw = pthGmmChcw./repmat(sum(pthGmmChcw,1),nth,1);
    pthGmmChcw(isnan(pthGmmChcw)) = 0;

    pthGmmChccw = (pmmGth.*repmat(pthGC(2,:),nmm,1))';
    pthGmmChccw = pthGmmChccw./repmat(sum(pthGmmChccw,1),nth,1);
    pthGmmChccw(isnan(pthGmmChccw)) = 0;

    EthChcw = th * pthGmmChcw;
    EthChccw = th * pthGmmChccw;
    
    a = 1./gradient(EthChcw,dstep);
    if ~flagDecisionGiven
        % memory noise
        pmmGm = exp(-((MM_m-repmat(m, nmm, 1)).^2)./(2*stdMemory^2)); 
        pmmGm = pmmGm./(repmat(sum(pmmGm,1),nmm,1));   

        % attention marginalization: compute distribution only over those ms that lead to cw decision!
        pmmGthChcw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(1,:)',1,length(thetaStim)));
        pmmGthChcw = pmmGthChcw ./ repmat(sum(pmmGthChcw,1),nmm,1);
        b = repmat(a',1,length(thetaStim)) .* pmmGthChcw;        
    else
        b = repmat(a',1,length(thetaStim)) .* pmmGth(:, ismember(th, thetaStim));   
    end
    pthhGthChcw = interp1(EthChcw,b,th,'linear','extrap');
    pthhGthChcw(pthhGthChcw < 0) = 0; 
    % add motor noise
    pthhGthChcw = conv2(pthhGthChcw,pdf('norm',th,0,stdMotor)','same');
    pthhGthChcw(pthhGthChcw < 0) = 0; 
    
    a = 1./gradient(EthChccw,dstep);
    if ~flagDecisionGiven
        % attention marginalization: compute distribution only over those ms that lead to ccw decision!
        pmmGthChccw = pmmGm * (pmGth(:, ismember(th, thetaStim)).*repmat(PChGm(2,:)',1,length(thetaStim))); 
        pmmGthChccw = pmmGthChccw ./ repmat(sum(pmmGthChccw,1),nmm,1);
        b = repmat(a',1,length(thetaStim)) .* pmmGthChccw;        
    else
        b = repmat(a',1,length(thetaStim)) .* pmmGth(:, ismember(th, thetaStim));
    end  
    pthhGthChccw = interp1(EthChccw,b,th,'linear','extrap');
    pthhGthChccw(pthhGthChccw < 0) = 0; 
    % add motor noise
    pthhGthChccw = conv2(pthhGthChccw,pdf('norm',th,0,stdMotor)','same');
    pthhGthChccw(pthhGthChccw < 0) = 0; 
    
    pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
    pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);            

    % modify psychometric curve p(Chat|theta, Congruent) ~ p(Congruent| Chat, theta) * p(Chat|Theta)
    pCongruentGcwTh = sum(pthhGthChcw(th' >= 0, :));
    pCongruentGccwTh = sum(pthhGthChccw(th' <= 0, :));
    PChGtheta_lapse_new = PChGtheta_lapse .* [pCongruentGcwTh; pCongruentGccwTh];
    PChGtheta_lapse_new = PChGtheta_lapse_new ./ repmat(sum(PChGtheta_lapse_new, 1), 2, 1);

    % modify the estimate distribution p(thetaHat|theta, Chat, Congrudent)
    pthhGthChccw(th'>= 0, :) = 0;
    pthhGthChcw(th'< 0, :) = 0;
    pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
    pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);   

    pthhGthChcw = pthhGthChcw./repmat(sum(pthhGthChcw,1),nth,1); % normalize - conv2 is not    
    pthhGthChccw = pthhGthChccw./repmat(sum(pthhGthChccw,1),nth,1);        
    mthhGthChcw = th * pthhGthChcw;
    mthhGthChccw = th * pthhGthChccw;

    pthhGth = pthhGthChcw.*repmat(PChGtheta_lapse_new(1,:),nth,1) + pthhGthChccw.*repmat(PChGtheta_lapse_new(2,:),nth,1);
    
    %% plot
    showrange = [-21 21];

    subplot(3,3,1);
    plot(thetaStim,PChGtheta_lapse(1,:),'r-');
    hold on;
    if kk==1
        plot([rangeth(1) rangeth(2)],[0.5 0.5],'k:');
        plot([0 0],[0 1],'k:');
    end
    axis([showrange(1) showrange(2) 0 1]);
    xlabel('Stimulus orientation (deg)')
    ylabel('Fraction "cw"')
    set(gca, 'FontSize', fontSize)
    
    subplot(3,3,2+(kk-1)*3);
    pthhGth = max(pthhGth(:)) - pthhGth;
    xRange = [-22 22];
    indX = find(thetaStim >= xRange(1) & thetaStim <= xRange(2));
    xMax = length(indX);
    yRange = [-40 40];
    indY = find(th >= yRange(1) & th <= yRange(2));
    yMax = length(indY);
    thNew = th(indY);
    indYStart = find(thNew == xRange(1));
    indYEnd = find(thNew == xRange(2));
    imagesc(pthhGth(indY, indX));
    hold on;
    axis xy;
    colormap('gray');
    plot([1 xMax],[round(yMax/2) round(yMax/2)],'k:', 'LineWidth', 1);
    plot([round(xMax/2) round(xMax/2)],[1 yMax],'k:', 'LineWidth', 1);
    plot([1 xMax],[indYStart indYEnd],'w--', 'LineWidth', 1.5);
    set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
        'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
        'FontSize', fontSize)
    xlabel('Stimulus orientation (deg)')
    ylabel('Estimated orientation (deg)')
    
    
    ax=subplot(3,3,kk*3);
    contour(pthhGth(indY, indX),5);
    colormap(ax,'parula');
    hold on;
    axis xy;
    plot([1 xMax],[round(yMax/2) round(yMax/2)],'k:', 'LineWidth', 1.5);
    plot([round(xMax/2) round(xMax/2)],[1 yMax],'k:', 'LineWidth', 1.5);
    plot([1 xMax],[indYStart indYEnd],'b:', 'LineWidth', 2);
    set(gca, 'ylim', [1 yMax], 'xlim', [1 xMax], ...
        'XTick', round(linspace(1,xMax,5)), 'XTickLabel', num2cell([-22 -11 0 11 22]),...
        'YTick', round(linspace(1,yMax,5)), 'YTickLabel', num2cell(round(linspace(yRange(1),yRange(2),5))), ...
        'FontSize', fontSize)
    xlabel('Stimulus orientation (deg)')
    ylabel('Estimated orientation (deg)')


    subplot(3,3,4);
    pthres = 0.075;
    ind = find(PChGtheta_lapse(1,:)>pthres);
    plot(thetaStim(ind),mthhGthChcw(ind),'c-','linewidth',2);
    hold on;
    ind = find(PChGtheta_lapse(2,:)>pthres);
    plot(thetaStim(ind),mthhGthChccw(ind),'g-','linewidth',2);
    axis([showrange(1) showrange(2) rangeth(1) rangeth(2)]);
    if kk==1
        plot(thetaStim,zeros(1,length(thetaStim)),'k:');
        plot([rangeth(1) rangeth(2)],[rangeth(1) rangeth(2)],'k--');
        plot([0 0],[rangeth(1) rangeth(2)],'k:');
    end
    xlabel('Stimulus orientation (deg)')
    ylabel('Estimated orientation (deg)')
    set(gca, 'FontSize', fontSize)

    subplot(3,3,7);
    pthres = 0.075;
    ind = find(PChGtheta_lapse(1,:)>pthres);
    plot(thetaStim(ind),mthhGthChcw(ind)-thetaStim(ind),'c-','linewidth',2);
    hold on;
    %ind = find(PChGth(2,:)>pthres);
    %plot(th(ind),mthhGthChccw(ind)-th(ind),'g-','linewidth',2);
    if kk==1
        plot(thetaStim,zeros(1,length(thetaStim)),'k:');
    end
    axis([0 showrange(2) -10 15]);
    xlabel('Stimulus orientation (deg)')
    ylabel('Bias (deg)')
    set(gca, 'FontSize', fontSize)
end
end

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
