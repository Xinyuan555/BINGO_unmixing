function [I_stack,MIP] = merge(I_unmix,tag,fluo,pixelx,pixely)
%% image display
I_color = I_unmix;
I_enhanced =I_unmix;
name = cell(1,fluo);
 for n = 1:fluo %逐一显示图像
%      if tag(n) == 1
%       I_unmix{n} = imadjust( I_unmix{n},[0 0.8],[]);
%      end
     gamma = 0.8;
     sigma = 0.75;
     [I_color{n},name{n}] = color_ge(I_unmix{n},tag(n));
 end
 I_stack = zeros(pixelx,pixely);
MIP = zeros(pixelx,pixely);
 for n = 1:fluo %逐一显示图像
     I_stack = I_stack + gamma.*I_color{n};
     MIP = max(MIP,I_color{n});
%      figure(12);
%      imshow(I_stack);
%      figure(13);    
%      imshow(MIP);     
 end
%  ima = im2double(I_stack);
%  imm = im2double(MIP);
% im_gr = rgb2gray(ima);
% im_low = imgaussfilt(im_gr,pixel(1));
% I_stack = ima./(im_low*5+eps);
% MIP = imm./(im_low*5+eps);
%  for i = 1:3
% %  MIP(:,:,i) = imgaussfilt3(MIP(:,:,i), sigma);
% %  I_stack(:,:,i) = imadjust(I_stack(:,:,i),[0.01 1]);
% %  I_stack(:,:,i) = imgaussfilt3(I_stack(:,:,i), sigma);
%  end
end