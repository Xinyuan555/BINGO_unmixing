%% nNMF
%inptut:
% V  image matrix
% r  number of endmembers
%
%output:
% W  abundance matrix
% H0 feature matrix (estimated spectra)




function [W,H] = nnmf_unmixing(V,r)

opt = statset('Maxiter',200,'TolFun',1e-7,'TolX ',1e-7,'Display','final');
[W,H] = nnmf(V,r,'replicates',1,'options',opt, 'algorithm','als');
end