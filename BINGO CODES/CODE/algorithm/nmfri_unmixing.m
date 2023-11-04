%%NMF-RI UNMIXING

function [W, H] = nmfri_unmixing(V,H0)
Y = V';
A0 = H0';
% % Data Preprocessing.
% [Y_sps, H_sps] = dataPreprocessing(A0,Y);
   
% Unmix data using NMF-RI
[AnmfRI, HmfRI] = NMF_RI(Y,A0);
W =  HmfRI';
H = AnmfRI';
end