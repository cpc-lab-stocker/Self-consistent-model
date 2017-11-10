
OPTBINS binning package v1.0

Developed by: Kevin H Knuth
Made public:  19 October 2006

It is respectfully requested that the reference below is cited in 
publications involving the use of this algorithm.

References:
Knuth K.H. Optimal data-based binning for histograms, 
http://arxiv.org/physics/0605197 
Please cite this paper to reference optBINS in publications.
Please examine arxiv.org for updated reference.

The test for excessive rounding was first published in:
Knuth K.H., Castle J.P., Wheeler K.R. 2006. Identifying excessively
rounded or truncated data. (Invited paper), Proceedings of the 
17th meeting of The International Association for Statistical Computing—
European Regional Section: Computational Statistics (COMPSTAT 2006), 
Rome, Italy, Aug 2006.

The code can be sped up by commenting out the references to histogram() and uncommenting histomex()


PACKAGE CONTENTS:

dispHIST	- displays a 1D histogram of a density model
dispMODEL	- displays a 1D density model
histMODEL	- develop a cell array model of a density model
histogram	- bins the data in multiple dimensions
histomex	- mex version of histogram
isROUNDED	- determines whether the data are excessively rounded
moptBINS	- brute force algorithm to determine the number of bins in multiple dimensions
nsoptBINS	- nested sampling algorithm to determine the number of bins in multiple dimensions
optBINS		- brute force algorithm to determine the number of bins for one dimension


RECOMMENDED STARTER CODE:

% playing with a Gaussian density

data = randn(1,1000);
m = OPTBINS(data, 20)           % try OPTBINS
m = mOPTBINS(data, 20)          % try mOPTBINS
m = nsOPTBINS(data)             % try nsOPTBINS

h1 = dispHIST(data, m, 'errorbars');
h2 = dispHIST(data, m, 'thick');

m = nsOPTBINS(data,'mean')      % try nsOPTBINS

m = nsOPTBINS(data,'verbose')	% try nsOPTBINS

model = histMODEL(data,nsOPTBINS(data));
dispMODEL(model,'errorbars');



% playing with a 2D Gaussian density

data = randn(2,5000);
m = mOPTBINS(data, [20, 20])    % try mOPTBINS
m = nsOPTBINS(data)             % try nsOPTBINS

% Make a 2D bar graph
[counts binwidth centers] = histogram(data, m);
figure; bar3(counts,1,'c');




% test to see if data are excessively rounded or truncated

data = rand(1,100);
m = nsOPTBINS(data)




% playing with a Uniform density

data = rand(1,1000);
m = OPTBINS(data, 20)


data = rand(2,1000);
m = nsOPTBINS(data)




