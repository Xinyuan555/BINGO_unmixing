%%Linear unmixing
%
%LU(V,H0)
%
%inptut:
% V  image matrix
% H0 feature matrix (reference spectra)
%
%output:
% W  abundance matrix


function [W,H] = LU(V,H0)


W_lsq = zeros(size(V,1),size(H0,1));
for i = 1:size(V,1)
    mid = lsqnonneg(H0',V(i,:)');
    W_lsq(i,:) = mid';
end
W = W_lsq;
H = H0;
W = imnoise(W,"gaussian",0,1e-3);
W = imnoise(W,"poisson");
end