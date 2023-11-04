%% SIMI
%%SIMI unmixing
%
%SIMI(V,H0)
%
%inptut:
% V  image matrix
% H0 feature matrix (reference spectra)
%
%output:
% W  abundance matrix

function [W,H] = simi_unmixing(V,H0)
W = zeros(size(V,1),size(H0,1));

RD = zeros(1,size(H0,1));
for i = 1: size(V,1)
    gray = V(i,:);
    normalized_gray=gray/(max(gray)+eps);
    for j=1:size(H0,1) 
         RD(j)=sum(((normalized_gray-H0(j,:)).^2));
    end
         [FP_R,FP]=min(RD);
         W(i,FP)=max(gray);
                   
end
W = imnoise(W,"gaussian",0,1e-3);
W = imnoise(W,"poisson");
H = H0;
end