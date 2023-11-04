function [W,H] = SVD(A,k)
%
% This function implements the NNDSVD algorithm described in [1] for
% initialization of Nonnegative Matrix Factorization Algorithms.
%
% [W,H] = nndsvd(A,k,flag);
%
% INPUT
% ------------
%
% A    : the input nonnegative m x n matrix A
% k    : the rank of the computed factors W,H
% flag : indicates the variant of the NNDSVD Algorithm
%        flag = 0 --> NNDSVD
%        flag = 1 --> NNDSVDa
%        flag = 2 --> NNDSVDar
%
% OUTPUT
% -------------
%   
% W   : nonnegative m x k matrix
% H   : nonnegative k x n matrix
%


%----------------------check the input matrix------------------------------
if numel(find(A<0)) > 0
    error('The input matrix contains negative elements !')
end
%--------------------------------------------------------------------------

% SVD --> partial SVD rank-k to the input matrix A. 
[U,S,V] = svds(A,k);

%the matrices of the factorization
W = pos(U)*sqrt(S);
H = sqrt(S)*pos(V);
end
%end of the svd initialization



%This function sets to zero the negative elements of a matrix
%--------------------------------------------------------------------------
function [Ap] = pos(A)
Ap = (A>=0).*A;
end
%--------------------------------------------------------------------------
