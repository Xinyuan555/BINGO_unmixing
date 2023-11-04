function [sparseness,L0,L1,L2] = sparseness(X)
% X可以是向量或是矩阵
% th = max(max(X))*0.1;
% X(X>th) = 1;
[m,n] = size(X);
num = m*n;
%L0范数
nz=(X~=0);
s0=sum(nz(:));
% 分别计算L1和L2范数
s1 = 0;
s2 = 0;
for i=1:m
    for j=1:n
        s1 = s1+X(i,j);
        s2 = s2+X(i,j)^2;
    end
end

% 计算稀疏度
s2 = sqrt(s2);
c = s1/s2;
a = sqrt(num)-c;
b = sqrt(num)-1;
sparseness = a/b;
L0=s0;
L1=s1;
L2=s2;
end