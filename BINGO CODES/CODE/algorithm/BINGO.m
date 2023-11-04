%% BINGO UNMIXING
%inptut:
% V    image matrix
% r    number of endmembers
% spH  sparseness value
%
%output:
% W  abundance matrix
% H0 feature matrix (estimated spectra)

function [W,H] = BINGO(V,r,spH)
options.lambda = size(V,1);
options.verbose = 2;
options.max_epoch = 100;
options.sH = spH;
% options.stepsizeH = 1e-1;
options.x_init.H = [];

fprintf(' Sparseness Setting of H: spH=%.5f\n',options.sH);
[X,~] = nmf_sc(V,r,options); 
W = X.W./(max(max(X.W)));
H = X.H;
end