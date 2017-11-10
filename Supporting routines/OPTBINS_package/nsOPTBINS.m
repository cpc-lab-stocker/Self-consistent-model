% nsOPTBINS
%
% Usage:    M = nsOPTBINS(data, 'options');
% Where:    M is the number of bins
%           data is a (D,N) vector of data points
%           options (described below) may be appended in any order
%
% Options:  'max'  returns the most probable number of bins (default)
%           'mode' same as 'max'
%           'mean' returns the mean number of bins
%           'verbose' produces verbose comments during the analysis
%           'diagnostics' produces diagnostic comments and plots
%
%
%	nsOPTBINS(data) returns the optimal (mode) number of bins
%   nsOPTBINS(data, 'mean') returns the mean number of bins
%   nsOPTBINS(data, 'verbose') displays progress statements
%   nsOPTBINS(data, 'diagnostics') displays diagnostic information
%   nsOPTBINS(data, 'mean', 'diagnostics') returns the mean number of bins 
%       and displays diagnostic information
%
%
% Created by:   Kevin Knuth
%           	2 August 2006 
%
% Last Modified: Kevin Knuth
%                16 Oct 2006
%
% This function uses nested sampling exploration with a copy operation
% to explore the posterior probability of the number of bins.  The data set
% can be multi-dimensional.  Memory availability is the only constraint.
% Average circumstances result in a limit of about 4-dimensional data.  The
% default operation is to return the mode (most probable solution) as the 
% optimal number of bins. The optimal number of bins is found when all 
% samples have converged at the same number of bins.  Nested sampling is 
% used in this case so that all objects may converge to the same solution 
% after adequate exploration of the logP space.  In the event that the mean
% is chosen, the samples drawn from the nested sampling algorithm will
% result in the mean number of bins along with the standard deviation.
%
% To speed up the algorithm, we store the log likelihood for evaluated bins
% making repeated binning and evaluation unnecessary.  However, the program
% may search many many bins before converging.  For this reason, we will
% only store a limit of L = 10000 bins, which are recycled based on when 
% last updated.
%
% The nested sampling algorithm was modified from the code in Sivia and
% Skilling, "Data Analysis: A Bayesian Tutorial", 2nd Ed., Oxford Univ.
% Press, 2006, pp. 188-189.  The evidence is not computed due to the fact
% that an artificial cutoff for the maximum number of bins is applied (see
% maxM below) and the fact that excessively rounded data produce
% stereotypic logP behavior that is difficult to manage (see Knuth, Castle
% and Wheeler 2006 below).
%
% It is respectfully that requested reference #1 below is cited in publications
% involving the use of this algorithm.
%
% References:
% 1. Knuth K.H. Optimal data-based binning for histograms, 
% http://arxiv.org/physics/0605197 
% Please cite this paper to reference optBINS in publications.
% Please examine arxiv.org for updated reference.
%
% 2. Knuth K.H., Castle J.P., Wheeler K.R. 2006. Identifying excessively
% rounded or truncated data. (Invited paper), Proceedings of the 
% 17th meeting of The International Association for Statistical Computing—
% European Regional Section: Computational Statistics (COMPSTAT 2006), 
% Rome, Italy, Aug 2006.



function optM = nsOPTBINS(data, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%% set up variables

% preliminaries
MEAN = 1;
MAX = 2;
EMPTY = 0;
verbose = 0;                % control verbose reporting
diagnostics = 0;            % control diagnostic evaluation
mode = EMPTY;
[D, N] = size(data);        % Dimensions and Number of data points
computedLOGP = 0;           % number of times logP was computed and stored


k = 1;
while k <= length(varargin)
    switch lower(varargin{k})
        case 'diagnostics'
            diagnostics = 1;
            verbose = 1;
        case 'verbose'
            verbose = 1;
        case 'mean'
            if mode == EMPTY
                mode = MEAN;
            else
                error('OPTBINS:multiplespecs','Input arguments set variable twice');
            end
        case 'max'
            if mode == EMPTY
                mode = MAX;
            else
                error('OPTBINS:multiplespecs','Input arguments set variable twice');
            end
        case 'mode'
            if mode == EMPTY
                mode = MAX;
            else
                error('OPTBINS:multiplespecs','Input arguments set variable twice');
            end
        otherwise
            error('OPTBINS:invalidargument','Argument not recognized');
    end
    k = k+1;
end
if mode == EMPTY    % set the mode to its default value
    mode = MAX;
end



% range of bins to consider...be somewhat conservative
maxM = ceil(5*(N^(1/3)))*ones(1,D);   % Maximum number if bins to consider

% mcmc variables
% the number of objects will be between 10 and 100, but we will aim 1/10
% along each dimension
n = min(100,max(10,ceil(prod(maxM)/(10^D))));       % Number of Objects
MCMC = 10;                % MCMC counter (pre-judged number of steps)
if diagnostics
    disp([blanks(5) 'There are ' num2str(n) ' objects']);
    disp([blanks(5) 'There are ' num2str(MCMC) ' intermediate MCMC steps']);
end

% number of bin values we will store
% chosen to be large, but not too large to swamp someone's memory
L = 10000;

% maximum number of iterations in the event that the mean is desired
MAXITERS = D*10000;
if diagnostics
    disp([blanks(5) 'Storage for ' num2str(L) ' logP values']);
    if mode == MAX
        disp([blanks(5) 'Maximum of ' num2str(MAXITERS) ' iterations']);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%% set up the Objects
fieldnames = {'m', 'logP', 'logWt'};
f = size(fieldnames,2);

cObj = cell([n, f]);
Obj = cell2struct(cObj, fieldnames, 2);
clear cObj;

cSamples = cell([n, f]);
Samples = cell2struct(cSamples, fieldnames, 2);
clear cSamples;

%%%%%%%%%%%%%%%%%%%%%%%%%%% set initial objects
% in the multi-dimensional case, we keep a list of all the bin numbers that
% were tried, and search it each time we want to take a step.
% We will store at most L bin values
% and will limit ourselves to the appropriate size integer for speed.

% binarray is an array containing the coordinates of the bins that have
% been tried.  Its index corresponds to logP from logParray.
% binarray(:,1) = trial last visited (used to get rid of old values)
% binarray(:,2:D+1) = bin dims

binarray = zeros(L,D+1,'uint16');
logParray = zeros(1,L);

H = 0.0;            % Information, initially 0
logZ = -realmax;    % Evidence Z, initially 0 (use exp(-realmax))
logwidth = log(1.0 - exp(-1.0/n)); % Outermost interval of prior mass

% initialize objects
% since we are not computing the evidence, I will set a few key values
% manually.  This will speed things up and help ensure that they are
% investigated.
iter = 1;

if verbose
    disp('Initializing Objects');
end

for i = 1:n

    if ((i > 5) || strcmp(upper(mode),'MEAN'))  % if we are to report the mean
        m = ceil(maxM.*rand(1,D));  %   choose bin numbers uniformly
    else                            % if we are to report the max
        m = ceil(maxM.^((i-1)/i));  %   choose a few values along the diag
    end
    
    % work with logP
    Obj(i).m = m;
    bindex = isEVALUATED(Obj(i).m, binarray);
    if bindex                               % if it has been saved
        Obj(i).logP = logParray(bindex);    %   look up logP
        binarray(index,1) = iter;           %   save the iteration number
    else                                    % otherwise
        
        computedLOGP = computedLOGP + 1;
        
        % Bin the data (equal width bins here)
        % to use mex file, comment out the line below and uncomment the
        % following line
        counts = histogram(data, m);        %   compute likelihood   
%       [counts, binwidth, centers] = histomex(data, m);  

        mm = prod(m);
        Obj(i).logP = N*log(mm) + gammaln(mm/2) - gammaln(N+mm/2) - ...
            mm*gammaln(1/2) + sum(reshape(gammaln(counts+0.5),1,mm));
        
                                            % save likelihood in table
        [trial index] = min(binarray(:,1)); % get next open index
        binarray(index,2:D+1) = Obj(i).m;   %   save the bin dimensions
        binarray(index,1) = iter;           %   save the iteration number
        logParray(index) = Obj(i).logP;     %   save logP
    end
end
if diagnostics
    disp([blanks(5) 'There are ' num2str(index-1) ' logP values stored']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% MCMC Loop
if verbose
    disp('MCMC Loop');
end

converged = 0;
notdone = 1;
while notdone

    iter = iter+1;
    
    % find worst object in collection
    worst = 1;
    for i = 2:n
        if (Obj(i).logP < Obj(worst).logP)
            worst = i;
        end
    end

    if diagnostics
        disp('---------------');
        disp(['Worst = ' mat2str(worst)])
        for i = 1:n
            disp(Obj(i).m);
        end
        disp('---------------');

        disp(['logwidth = ' num2str(logwidth)]);
        disp(['worst logP = ' num2str(Obj(worst).logP)]);
    end
    
    
    % if we want the mean, we have to keep track of all this
    if mode == MEAN
        % Weight = width * Likelihood
        Obj(worst).logWt = logwidth + Obj(worst).logP;
        
        if diagnostics
            disp(['iter = ' num2str(iter)]);
            disp(['logwidth = ' num2str(logwidth)]);
            disp(['worst logP = ' num2str(Obj(worst).logP)]);
            disp(['logWt = ' num2str(Obj(worst).logWt)]);
        end

        % Update evidence Z and information H
        if (logZ > Obj(worst).logWt)
            logZnew = logZ + log(1 + exp(Obj(worst).logWt-logZ));
        else
            logZnew = Obj(worst).logWt + log(1 + exp(logZ - Obj(worst).logWt));
        end
        H = exp(Obj(worst).logWt - logZnew) * Obj(worst).logP + ...
            exp(logZ - logZnew) * (H + logZ) - logZnew;
        logZ = logZnew;
    end
    
    %%%%%% store posterior sample
    Samples(iter-1) = Obj(worst);       % we start the iters at 2
    
    if verbose
        if rem(iter,10) == 0
            disp(['Iteration = ' num2str(iter)]);
            disp(['     Bins = ' mat2str(Samples(iter-1).m) ' and LogP = ' num2str(Samples(iter-1).logP)]);
        end
    end
    
    %%%%%% kill worst object in favor of a copy of a different survivor
    copy = ceil(n * rand());        % choose an object between 1 and n
    while ((copy == worst) && n>1)  % if it is the same as it was before
        copy = ceil(n*rand);        % then choose another object
    end
    logPstar = Obj(worst).logP;     % new logP constraint
    Obj(worst) = Obj(copy);         % overwrite worst object

    if diagnostics
        disp(['ccd:     Bins = ' mat2str(Obj(worst).m) ' and LogP = ' num2str(Obj(worst).logP)]);
    end
    
    
    %%%%%% Evolve copied object
    accept = 0;         % number of MCMC acceptances
    reject = 0;         % number of MCMC rejections
    step = 10*ones(1,D);          % initial step-size
    
    for i = 1:MCMC
        
        if diagnostics
            disp(['step = ' num2str(step)]);
        end
        
        %%%%% try a move
        notvalid = 1;
        while notvalid
            m = ceil(Obj(worst).m + step.*(2.0*randn(1,D)-1.0));  % |move| < step
            m = rem(m, maxM);       % wraparound to stay within range

            % make sure that no dimensions are negative
            for d = 1:D
                if m(d) <= 0
                    break;
                end
                if d == D
                    notvalid = 0;   % it is within the range
                end
            end
        end
        
        if diagnostics
            disp(['try m = ' num2str(m)]);
        end
        
        % work with logP
        bindex = isEVALUATED(m, binarray);
        if bindex                           % if it has been saved
            logP = logParray(bindex);       %   look up logP
            binarray(index,1) = iter;       %   save the iteration number
        else                                % otherwise
            computedLOGP = computedLOGP + 1;
            
            counts = histogram(data, m);    %   compute likelihood
%           [counts, binwidth, centers] = histomex(data, m);
            mm = prod(m);
            logP = N*log(mm) + gammaln(mm/2) - gammaln(N+mm/2) - ...
                mm*gammaln(1/2) + sum(reshape(gammaln(counts+0.5),1,mm));

            % save likelihood in table
            [trial index] = min(binarray(:,1));  % get next open index
            binarray(index,2:D+1) = m;        % save the bin dimensions
            binarray(index,1) = iter;         % save iteration number
            logParray(index) = logP;          % save logP
        end
        
        % accept only within hard likelihood constraint
        if (logP > logPstar)
            Obj(worst).m = m;
            Obj(worst).logP = logP;
            accept = accept + 1;
        else
            reject = reject + 1;
        end

        % refine step size to let acceptance ratio converge around 50%
        if (accept > reject)
            step = step * exp(1.0/accept);
        end
        if (accept < reject)
            step = step / exp(1.0/reject);
        end
        
    end
    
    % Test for convergence
    if rem(iter,10) == 0     % if we desire the maximum 
        m = Obj(1).m;
        for i = 2:n
            test = sum(m == Obj(i).m);
            if test ~= D
                break;
            end
            if i == n           % the samples have converged
                converged = converged + 1;    % we have converged
            end
        end
    end
    
    % See if we should end early
    if mode == MAX               % if we want the maximum and we have converged
        if converged
            notdone = 0;
        end
    else     % if we desire the mean: test if we are done
        if (iter == MAXITERS) || (converged == 10)
            notdone = 0;        % Here we go a bit longer to ensure exploration
        end
        logwidth = logwidth - 1.0/n;    % shrink the interval
    end

end


% Find the solution
% For the MAX, this is easy
% For MEAN we need to compute a weighted sum
if mode == MAX     % if we desire the maximum
    optM = Obj(1).m;
else                       % if we desire the mean
    optM  = zeros(1,D);
    mu2   = zeros(1,D);
    w = zeros(1,iter-1);
    for i = 1:iter-1  % in C indices go from 0:nest-1, here from 1:nest
        w(i) = exp(Samples(i).logWt - logZ);   % proportional weight
        if diagnostics
            disp(['w = ' num2str(w(i))])
        end
        optM  = optM  + w(i) * Samples(i).m;           % weighted sum
        mu2 = mu2 + w(i) * (Samples(i).m .* Samples(i).m);
    end
    stddev = round(100*sqrt(mu2-optM.^2))/100;
    optM = round(optM);

    if diagnostics
        figure; 
        plot(1:iter-1, w);
    end
end

% Verbose Display
if verbose
    if mode == MAX     % if we desire the maximum
        disp(['optM = ' mat2str(optM)]);
    else
        disp(['optM = ' mat2str(optM) ' +- ' mat2str(stddev)]);
    end

    if isempty(find(optM,1))   % only display if all bins are > 1
        % display if 1d
        if D == 1
            figure;
            hist(data,optM);
        end

        % display if 2d
        if D == 2
            counts = histogram(data, optM);
            figure;
            bar3(counts,1,'c');
        end
    end
    
    if mode == MEAN
        disp(['Number of iterations = ' num2str(iter)]);
        disp(['Evidence: ln(Z) = ' num2str(logZ) ' +- ' num2str(sqrt(H/n))]);
        disp(['Information: H = ' num2str(H) 'nats = ' num2str(H/log(2)) 'bits']);
    end
end

% Test to see if data are excessively rounded or truncated
rounded = isROUNDED(data, m);
if rounded
    disp('The data are excessively rounded.');
    warning('Add a small number to each datum vector to obtain a reasonable result.');
end

disp(['Computed ' num2str(computedLOGP) ' logP values.']);

return





% isEVALUATED
% isEVALUATED determines whether this binning strategy has been evaluated
% recently
%
% Usage:
%           flag = isEVALUATED(M, binarray)
%           
% Where:
%           M  - is the array of the number of bins
%           binarray - is the array of bin values evaluated recently
%           flag - 1 for yes, 0 for no
%
% Created by:   Kevin Knuth
%           	10 Sept 2006

function result = isEVALUATED(m, binarray)

D = length(m);

d = 1;
result = find(binarray(:,d+1) == m(d));
for d = 2:D
    if isempty(result) 
        break;
    else
        result = intersect(result, find(binarray(:,d+1) == m(d)));
    end
end

if isempty(result)
    result = 0;
else
    if length(result) > 1
        warning('Multiple entries in binarray');
    end
end

return


