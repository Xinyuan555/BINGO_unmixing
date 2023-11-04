function [ss, cor2, HistDist, dice, MIoU] = similarity(I_in,I_out)
ss = zeros(1,size(I_in,2));
cor2 = zeros(1,size(I_in,2));
HistDist = zeros(1,size(I_in,2));
dice = zeros(1,size(I_in,2));
MIoU = zeros(1,size(I_in,2));
confu = zeros(2,2,size(I_in,2));
for i = 1:size(I_in,2)
    %SSIM
    ss(i) = ssim(I_in{i},I_out{i});
    %corr2
    cor2(i) = corr2(I_in{i},I_out{i});
    %histdist
   [c1,~]=imhist(I_in{i});
   c1=c1/size(I_in{i},1)/size(I_in{i},2);
   [c2,~]=imhist(I_out{i});
   c2=c2/size(I_out{i},1)/size(I_out{i},2);
SumI=sum(c1);
 SumJ=sum(c2);
 %Sumup = sqrt(Count1.*Count2);
 %SumDown = sqrt(Sum1*Sum2);
 Sumup = sqrt(c1.*c2);
 SumDown = sqrt(SumI*SumJ);
 Sumup = sum(Sumup);
 HistDist(i)=sqrt(1-Sumup/SumDown);
 %dice
 bw1 = imbinarize(I_in{i});
 bw2 = imbinarize(I_out{i});
 inters = nnz(bw1&bw2);
 union = nnz(bw1)+nnz(bw2);
 dice(i) = 2*inters/union;
 %MIoU
TP = inters;
FP = nnz(~bw1&bw2);
FN = nnz(bw1&~bw2);
TN = nnz(~bw1&~bw2);
confu(:,:,i) = [TP FN;FP TN];
MIoU(i) = TP/(TP+FN+FP);
end
end