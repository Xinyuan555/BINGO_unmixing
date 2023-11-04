function [I_unmix_order,name,tag,SAD,SID,RMSE,corr,HN] = evluation(V,W,H,H0,inimg,mark) 
name = cell(1,size(W,2));
 img_num = size(H,2);
 r = size(H,1);
 pixel = sqrt(size(W,1));
 W_order = zeros(size(W));
 H_order = zeros(size(H));
 x = 1:img_num; 
L = zeros(1,r);
 Lr =zeros(r,r);
 SAD = zeros(1,r);
 SID=zeros(1,r);
tag = zeros(1,r);
 %% RMSE
  RE_V = W*H;
  RMSE = sum(sum((V-RE_V).^2))/((sqrt(img_num))*pixel*pixel);
%% reconstruction
I_unmix=cell(1,r);I_unmix_order=cell(1,r);
I_RE = cell(1,img_num);
PT=zeros(pixel,pixel);
corr_af = zeros(size(W,2),size(W,2));
        for i = 1:size(W,2)
            for j = i:size(W,2)
                a = corrcoef(W(:,i),W(:,j));
                corr_af(i,j) = a(1,2);
            end
        end
for n=1:r
    for j=1:pixel
        PT(:,j)=W((j-1)*pixel+1:pixel*j,n);
    end
    PT = PT/max(max(PT));
    I_unmix{n}=uint8(round(PT*255));%将data转化至0-255之间
end
PT=zeros(pixel,pixel);
for n=1:img_num
    for j=1:pixel
        PT(:,j)=RE_V((j-1)*pixel+1:pixel*j,n);
    end
    I_RE{n}=uint8(round(PT*255));%将data转化至0-255之间
end
 %% precision evaluation
 

if mark == 1
%% discriminate the fluo
for i = 1:r
    [L(i),tag(i)] = max(H(i,:));
end

 

  a = unique(tag);
  for n = 1:size(a,2)
      y = find(tag==a(n));
      b = 2;
      while size(y,2) > 1
          sub = size(y,2);
          L_or = zeros(sub,size(H,2));
          for m = 1:sub
          L_or(m,:) = sort(H(y(m),:),'descend');
          end
          [~,inn] = min(L_or(:,b));
          y(inn) = [];
          L_or(inn,:) = []; 
          for m = 1:sub-1
          [s,~] = max(L_or(m,b:end));
          tag(y(m)) = find(H(y(m),:)==s);
          end
          b = b+1;
       end
  end 
  for i = 1: r
      I_unmix_order{tag(i)} = I_unmix{i};
      H_order(tag(i),:) = H(i,:);
  end
%   W = W_order;
  H = H_order;
% %% DISCRIMINATE MANUALY
%     figure
% subplot(2,6,1);imshow(I_unmix{1});
% subplot(2,6,2);imshow(I_unmix{2});    
% subplot(2,6,3);imshow(I_unmix{3});
% subplot(2,6,4);imshow(I_unmix{4});
% subplot(2,6,5);imshow(I_unmix{5});
% subplot(2,6,6);imshow(I_unmix{6});
% subplot(2,6,7);imshow(I_unmix{7});
% subplot(2,6,8);imshow(I_unmix{8});
% subplot(2,6,9);imshow(I_unmix{9});
% subplot(2,6,10);imshow(I_unmix{10});
% 
%     for i = 1:r
% 
%      prompt = "What is the tag? ";
%       t(i) = input(prompt);
%     end
%   for i = 1: r
%       I_unmix_order{t(i)} = I_unmix{i};
%       H_order(t(i),:) = H(i,:);
%   end
%   H = H_order;
%  
% 
else
    I_unmix_order = I_unmix;
end
  tag = 1:r;
  
  for i = 1: r
      Dab = 0;
      Dba = 0;
     H(i,:) = H(i,:)-min(H(i,:));
     H(i,:) = H(i,:)/max(H(i,:))+eps;
  SAD(i)=acos(H0(i,:)*H(i,:)'/(norm(H0(i,:),2)*norm(H(i,:),2)));


%   计算光谱信息散度

  for j = 1:r
      for l = 1:r
      p = H0(tag(i),l)/sum(H0(tag(i),:))+eps;
      q = H(i,l)/sum(H(i,:))+eps;
      Dab = Dab + p*log(p/q);
      Dba = Dba + q*log(q/p);
      end
      SID(i) = Dab + Dba;
  end
%   [SID(i),tag(i)] = min(SIDr(i,:));
 %记录精度 
 switch tag(i)
    case 1
        name{i} = 'blue';
     case 2
        name{i} = 'bluegreen';
    case 3
        name{i} = 'green';
    case 4
        name{i} = 'yellowgreen';
    case 5
        name{i} = 'yellow';
    case 6
        name{i} = 'orange';
    case 7
        name{i} = 'redorange';
    case 8
        name{i} = 'red';
    case 9
        name{i} = 'crimson';
    case 10
        name{i} = 'scarlet';
end
 end
SAD_AVE = mean(SAD);SID_AVE = mean(SID);
HN = H;
 %% RMSE
  RE_V = W*H;
  RMSE = sum(sum((V-RE_V).^2))/((sqrt(img_num))*pixel*pixel);

%% caculate corrcoef
corr = zeros(size(W,2),size(W,2));
        for i = 1:size(W,2)
            for j = i:size(W,2)
              corr(i,j) = corr2(I_unmix_order{i},I_unmix_order{j});

            end
        end
%        corr = mean(mean(corr-diag(diag(corr))));
I_unmix_order = subb(I_unmix_order);
for  i = 1:r
    I_unmix_order{i} = imadjust(I_unmix_order{i},[0.01,1],[]);
end
end