close all
clear;
clc;

addpath(pwd);
ori_path = pwd;
addpath functions/;
addpath algorithms/;
addpath(genpath(pwd));

%% input data
file_path0 = '../DEMO/BEADS/';
num_fluo = 10;
rect = cell(1,num_fluo);
gtfolder = 'ground truth';rawfolder = 'raw';
load('Fluosphere.mat');
folderinfo0 = dir(fullfile(file_path0,rawfolder,'*jpg*'));
gtinfo = dir(fullfile(file_path0,gtfolder,'*tif*'));
img_num = length(folderinfo0);

left = 400;
right = 700;
max_ch = img_num;
min_bindwith = 4;
H0 = zeros(num_fluo,img_num);
Hm = zeros(num_fluo,img_num);
for gap = 1%:4
    if img_num/gap == fix(img_num/gap)
        channel_num = fix(img_num/gap);
        bandwidth = fix((right-left)/channel_num/gap);
        if bandwidth<min_bindwith
            break;
        end
I=cell(1,channel_num); 
I_color=cell(1,channel_num);  
H0_NOR = zeros(num_fluo,channel_num);
Hm = zeros(num_fluo,channel_num);
spH = 0.1;%sparseness setting（adjustable）
TH = 0.00;%threshold setting of display（adjustable）
spectra = cell(1,num_fluo);
filter = cell(1,channel_num);
signal = zeros(num_fluo,channel_num);
spectra{1} = blue;
spectra{2} = bluegreen;
spectra{3} = green;
spectra{4} = yellowgreen;
spectra{5} = yellow;
spectra{6} = orange;
spectra{7} = redorange;
spectra{8} = red;
% spectra{8} = carmine;
spectra{9} = crimson;
spectra{10} = scarlet;
snr = zeros(1,channel_num);
%% calculate spectra overlap
overlap = zeros(1,num_fluo-1);
self = zeros(1,num_fluo-1);
for i =  1: num_fluo-1
    for x = 1:450
        self(1,i) = self(1,i)+spectra{i}(x,2);
      if spectra{i}(x,2)>=spectra{i+1}(x,2)
          overlap(1,i) = overlap(1,i)+spectra{i+1}(x,2);
      else
          overlap(1,i) = overlap(1,i)+spectra{i}(x,2);
      end
    end
    overlap(1,i) =overlap(1,i)/self(1,i); 
end
%% Pseudo-color
if channel_num == 3
    map_index=[0 0 1 ; 0 1 0 ; 1 0 0];
    for i =1:channel_num
    map{i} = customcolormap([0 0.5 1], [1 1 1; map_index(i,:); 0 0 0]);
    end
else
map = cell(1,channel_num);
map_index = zeros(channel_num,3);
Nc = ceil(channel_num/6);
    map_index(1,:) = [0.5 0 1];
    map_index(2,:) = [0 0 1];
    map_index(3,:) = [0 1 1];
    map_index(4,:) = [0 1 0.5];
    map_index(5,:) = [0 1 0];
    map_index(6,:) = [0.5 1 0];
    map_index(7,:) = [1 1 0];
    map_index(8,:) = [1 0.5 0];
    map_index(9,:) = [1 0 0];
    map_index(10,:) = [1 0 1];
    for i = 1:channel_num
    map{i} = customcolormap([0 0.5 1], [1 1 1; map_index(i,:); 0 0 0]);
    end
end

%% Filter setting
    channel = left:right;
    
    for i = 1:channel_num
        mid = left+fix(((right-left)/channel_num)/2)*(2*(i-1)+1);
        filter{i} = mid-fix(bandwidth/2):mid+fix(bandwidth/2);
    end

 for j = 1:num_fluo
     fluo = spectra{j};
     signal(j,:) = fluorescence_collection(fluo,filter);
 end

 for m = 1:num_fluo
     signal(m,:) = signal(m,:)./max(max(signal(m,:)));
 end
 H0 = signal;

 figure(1);
for i = 1: num_fluo
    [~,p(i)] = max(H0(i,:));
    peak(i) = wavelength(p(i));peak = peak';
    plot(spectra{i}(:,1),spectra{i}(:,2),'color',map_index(i,:),'linewidth',1.5);
    title('Emission Spectra','FontSize',14);
    axis ([350 700 0 100]);
    hold on;
end
 %% load images
%load groundtruth
inimg = cell(1,num_fluo);inimg_color = cell(1,num_fluo);
for i = 1: num_fluo
    inimg{i} = imread(fullfile(file_path0,gtfolder,gtinfo(i).name));
end


h = msgbox('Loading images......');
I = cell(1,channel_num);
I_cos = cell(1,channel_num);
I_color = cell(1,channel_num);
ave = zeros(1,channel_num);
aveun = zeros(1,num_fluo);
class = cell(1,channel_num);
count_fluo = zeros(1,channel_num);
p = zeros(1,num_fluo);
x = 1:channel_num;
H_cla = zeros(channel_num,channel_num);
    for ch = 1:channel_num
        I_sum = zeros(512,512);
        for n = 1:gap
           I_R = cell(1,gap);
           I_R{n} = im2double(imread(fullfile(file_path0,rawfolder,folderinfo0(gap*(ch-1)+n).name)));
        I_sum = I_sum+I_R{n};
        end
          I_cos{ch} = I_sum;
          I_color{ch} = ind2rgb(im2uint8(im2gray(I_cos{ch})),map{ch});
          I_cos{ch} = im2double(I_cos{ch});
        if  ch == 1
             I_pre = I_color{ch};
             I_pre_mip = I_color{ch};
        else
             I_pre = I_pre+I_color{ch};
             I_pre_mip= max(I_pre_mip,I_color{ch});
        end
    end
    I_pre = imadjust(I_pre,[0.0 1],[]);
    I_pre = imgaussfilt(I_pre);
    close(h);

Hm = [0.0268297809629270	0.0453738942755383	0.0519404268011286	0.0465575610827263	0.0321562815952728	0.0233915107134763	0.0141758191432272	0.00975115988778649	0.00493194502994981	0.00741200881643886	0.00476284977178010	0.00476284977178010	0.00476284977178010	0.00476284977178011	0.00467830214269525	0.00707381830009945
0.000133226567042800	0.0494936696564003	0.348021099757555	0.440247190792934	0.325772263061407	0.226585083898043	0.147448503074619	0.0974552337918084	0.0622501134507484	0.0448307398099023	0.0288768584065269	0.0191513190124025	0.0102917523040563	0.00945908626003881	0.00476284977178011	0.00942577961827811
0	0.00349275649930541	0.00433948534762188	0.0156644836938546	0.163842032149236	0.306621684196600	0.333928689554806	0.281960706489383	0.214963286366342	0.164688760997553	0.115208043924059	0.0749884236290269	0.0424422835218628	0.0286829397367202	0.0164053714361315	0.0143414698683601
0	0.00466515028928205	0.00161204146121788	0.00476284977178010	0.00476284977178010	0.00537347153739294	0.0475307982353030	0.0764498450547269	0.0634313890118612	0.0490939899552719	0.0363442074892759	0.0274291297113285	0.0173905078846535	0.0134336788434824	0.00849985497733065	0.00945242493168667
0.000425254443908939	0.00456439769795593	0.00209792192328410	0.00476284977178010	0.00476284977178010	0.00476284977178010	0.0130978368723953	0.225526606753040	0.315822300343038	0.291100842003799	0.270490176622345	0.202165962634309	0.124060896436368	0.0893034332208770	0.0490743628270915	0.0312703767754372
0.00199469776766860	0.00223894647391373	0.000732746118735402	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00492568224261020	0.0514550607823082	0.338691539326586	0.692038001027880	0.658901593213957	0.363767739834419	0.219783127502913	0.117157962762249	0.0650922802143281
6.61506912747239e-05	0.000793808295296686	0.00224912350334061	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00780578157041741	0.0428656479460210	0.710723027055634	1	0.642918568499041	0.381557187272607	0.202024211153007
0.000272162844101721	0.00333399484024608	0.00129277350948317	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00874323136676777	0.0663737136053071	0.426717319195986	0.360139483457602	0.210415898846143	0.123595951577694	0.0659994896946673
3.96904147648343e-05	0.00341337566977574	0.00107164119865053	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00476284977178011	0.00373089898789442	0.00476284977178011	0.00476284977178011	0.00734272673149433	0.0216312760468347	0.293193093867831	0.518594959317324	0.331414963286366
0	0.00376143007617506	0.00258903628619842	0.00476284977178010	0.00476284977178010	0.00476284977178010	0.00476284977178010	0.00476284977178010	0.00398125391179568	0.00476284977178010	0.00476284977178010	0.00476284977178010	0.00886622803669835	0.0320454302593615	0.128572518967438	0.530825713282550];
for i = 1:num_fluo
[~,p(i)] = max(Hm(i,:));
end

%% 
figure(2);
Hm = Hm./(max(max(Hm)));
for i = 1: num_fluo
    plot(x,Hm(i,:),'color',map_index(i,:),'linewidth',1.5);
    title('Measured Spectra','FontSize',14);
    hold on;
end
hold off;

pixel = size(I_cos{ch},1);

h = msgbox('Start Unmixing......');

%% flattening
     p(1) = 2;
     Vor = flattening(I_cos,pixel,channel_num);
     V = Vor;
     H0or = H0;Hmor = Hm;I_raw = cell(1,num_fluo);
     V = zeros(size(Vor,1),num_fluo);
     H0 = zeros(num_fluo,num_fluo);
     Hm = zeros(num_fluo,num_fluo);
     for i = 1:num_fluo
         V(:,i) = Vor(:,p(i));
         H0(:,i) = H0or(:,p(i));
         Hm(:,i) = Hmor(:,p(i));
         I_raw{i} = im2gray(I_cos{p(i)});
     end

     para.H0 = H0;
     para.Hm = Hm;
    
     %% caculate corrcoef
     corr_raw = zeros(size(V,2),size(V,2));
        for i = 1:size(V,2)
            for j = i:size(V,2)
                a = corr2(I_raw{i},I_raw{j});
                corr_raw(i,j) = a;
            end
        end
       
%% UNMIXING 
%nmf
 tstart = tic;
 [W_BINGO,H_BINGO] = BINGO(V,num_fluo,spH);
 para.time_BINGO = toc(tstart);
  [I_unmix_BINGO,fname_BINGO,tag_BINGO,para.SAD_BINGO,para.SID_BINGO,para.RMSE_BINGO,para.corr_BINGO,para.H_BINGO] = evluation(V,W_BINGO,H_BINGO,Hm,inimg,1);
  [para.SSIM_BINGO,para.CORR2_BINGO,para.hist_BINGO,para.Dice_BINGO,para.MIoU_BINGO] = similarity(I_unmix_BINGO,inimg);

 %nnmf
 tstart = tic;
 [W_nnmf,H_nnmf] = nnmf_unmixing(V,num_fluo);
 para.time_nnmf = toc(tstart);
 [I_unmix_nnmf,fname_nnmf,tag_nnmf,para.SAD_nNMF,para.SID_nNMF,para.RMSE_nNMF,para.corr_nNMF,para.H_nNMF] = evluation(V,W_nnmf,H_nnmf,Hm,inimg,1);
 [para.SSIM_nNMF,para.CORR2_nNMF,para.hist_nNMF,para.Dice_nNMF,para.MIoU_nNMF] = similarity(I_unmix_nnmf,inimg);
 %linear unmixing(r)
 tstart = tic;
 [W_pinv,H_pinv] = LU(V,H0);
 para.time_pinv = toc(tstart);
 [I_unmix_pinv,fname_pinv,tag_pinv,para.SAD_PINV,para.SID_PINV,para.RMSE_PINV,para.corr_PINV] = evluation(V,W_pinv,H_pinv,Hm,inimg,0);
[para.SSIM_PINV,para.CORR2_PINV,para.hist_PINV,para.Dice_PINV,para.MIoU_PINV] = similarity(I_unmix_pinv,inimg);
 %linear unmixing(m)
 tstart = tic;
 [W_pinvm,H_pinvm] = LU(V,Hm);
para.time_pinvm = toc(tstart);
 [I_unmix_pinvm,fname_pinvm,tag_pinvm,para.SAD_PINVm,para.SID_PINVm,para.RMSE_PINVm,para.corr_PINVm] = evluation(V,W_pinvm,H_pinvm,Hm,inimg,0);
[para.SSIM_PINVm,para.CORR2_PINVm,para.hist_PINVm,para.Dice_PINVm,para.MIoU_PINVm] = similarity(I_unmix_pinvm,inimg);
 %SIMI(r)
 tstart = tic;
 [W_simi,H_simi] = simi_unmixing(V,H0);
para.time_simi = toc(tstart);
 [I_unmix_simi,fname_simi,tag_simi,para.SAD_SIMI,para.SID_SIMI,para.RMSE_SIMI,para.corr_SIMI] = evluation(V,W_simi,H_simi,Hm,inimg,0);
 [para.SSIM_SIMI,para.CORR2_SIMI,para.hist_SIMI,para.Dice_SIMI,para.MIoU_SIMI] = similarity(I_unmix_simi,inimg);
 %SIMI(m)
 tstart = tic;
 [W_simim,H_simim] = simi_unmixing(V,Hm);
para.time_simim = toc(tstart);
 [I_unmix_simim,fname_simim,tag_simim,para.SAD_SIMIm,para.SID_SIMIm,para.RMSE_SIMIm,para.corr_SIMIm] = evluation(V,W_simim,H_simim,Hm,inimg,0);
 [para.SSIM_SIMIm,para.CORR2_SIMIm,para.hist_SIMIm,para.Dice_SIMIm,para.MIoU_SIMIm] = similarity(I_unmix_simim,inimg);
%%
%PICASSO
tstart = tic;
I_unmix_picasso = picasso_unmixing(I_cos,p,num_fluo);
para.time_picasso = toc(tstart);
I_unmix_picasso = subb(I_unmix_picasso);
[para.SSIM_picasso,para.CORR2_picasso,para.hist_picasso,para.Dice_picasso,para.MIoU_picasso] = similarity(I_unmix_picasso,inimg);
 %% caculate corrcoef
        for i = 1:num_fluo
            for j = i:num_fluo
                a = corr2(I_unmix_picasso{i},I_unmix_picasso{j});
                corr(i,j) = a;
            end
        end
       para.corr_picasso = corr;

%NMF-RI
tstart = tic;
[W_ri, H_ri] = nmfri_unmixing(V,H0);
para.time_ri = toc(tstart);
[I_unmix_ri,fname_ri,tag_ri,para.SAD_RI,para.SID_RI,para.RMSE_RI,para.corr_RI,para.H_RI] = evluation(V,W_ri,H_ri,Hm,inimg,0);
 [para.SSIM_RI,para.CORR2_RI,para.hist_RI,para.Dice_RI,para.MIoU_RI] = similarity(I_unmix_ri,inimg);

%% image display
close(h);
 I_raw = cell(1,num_fluo);
for n = 1:num_fluo
      I_raw{n} = uint8(0.5*round(im2gray(I_cos{p(n)})./max(max(im2gray(I_cos{p(n)})))*255));
      inimg{n} = uint8(0.5*round(inimg{n}./max(max(inimg{n}))*255));

end
  [para.SSIM_RAW,para.CORR2_RAW,para.hist_RAW,para.Dice_RAW,para.MIoU_RAW] = similarity(I_raw,inimg);

[BINGOunmix_SUM,BINGOunmix_MIP] = merge(I_unmix_BINGO,tag_BINGO,num_fluo,pixel,pixel);
[nNMFunmix_SUM,nNMFunmix_MIP] = merge(I_unmix_nnmf,tag_nnmf,num_fluo,pixel,pixel);
[PINVunmix_SUM,PINVunmix_MIP] = merge(I_unmix_pinv,tag_pinv,num_fluo,pixel,pixel);
[SIMIunmix_SUM,SIMIunmix_MIP] = merge(I_unmix_simi,tag_simi,num_fluo,pixel,pixel);
[PINVmunmix_SUM,PINVmunmix_MIP] = merge(I_unmix_pinvm,tag_pinvm,num_fluo,pixel,pixel);
[SIMImunmix_SUM,SIMImunmix_MIP] = merge(I_unmix_simim,tag_simim,num_fluo,pixel,pixel);
[PICASSOunmix_SUM,PICASSOunmix_MIP] = merge(I_unmix_picasso,tag_BINGO,num_fluo,pixel,pixel);
[RIunmix_SUM,RIunmix_MIP] = merge(I_unmix_ri,tag_ri,num_fluo,pixel,pixel);
[GT_SUM,GT_MIP] = merge(inimg,tag_ri,num_fluo,pixel,pixel);
[Ipre_SUM,Ipre_MIP] = merge(I_raw,tag_ri,num_fluo,pixel,pixel);
h1 = figure(3);
     set(h1,'position',[0 0 1920 1080]);
     subplot(2,4,1);
     imshow(Ipre_SUM);
     title('Raw','FontSize',14);
     subplot(2,4,2);
     imshow(PINVunmix_SUM);
     title('Linear unmixing(r)','FontSize',14);
     subplot(2,4,3);
     imshow(PICASSOunmix_SUM);
     title('PICASSO','FontSize',14);
     subplot(2,4,4);
     imshow(BINGOunmix_SUM);
     title('BINGO','FontSize',14);
     subplot(2,4,5);
     imshow(SIMIunmix_SUM);
     title('SIMI','FontSize',14);
     subplot(2,4,6);
     imshow(RIunmix_SUM);
     title('NMF-RI','FontSize',14);
     subplot(2,4,7);
     imshow(nNMFunmix_SUM);
     title('nNMF','FontSize',14);
     subplot(2,4,8);
     imshow(GT_MIP);
     title('Ground Truth','FontSize',14);
 %% save unmixing images
 ori_path=pwd;%remember current path
 file_path = file_path0;
 cd(file_path);%open file path
  newfolder='Results';mkdir(newfolder) %creat fold
   save_path=strcat(newfolder);
 cd(save_path);%open fold

  for n = 1:num_fluo
    [I_unmix_nnmf{n},~] = color_ge(I_unmix_nnmf{n},n); 
    [I_unmix_BINGO{n},~] = color_ge(I_unmix_BINGO{n},n);
    [I_unmix_pinv{n},~] = color_ge(I_unmix_pinv{n},n);
    [I_unmix_simi{n},~] = color_ge(I_unmix_simi{n},n);
    [I_unmix_pinvm{n},~] = color_ge(I_unmix_pinvm{n},n);
    [I_unmix_simim{n},~] = color_ge(I_unmix_simim{n},n);
    [I_unmix_picasso{n},~] = color_ge(I_unmix_picasso{n},n);
    [I_unmix_ri{n},~] = color_ge(I_unmix_ri{n},n);
    [I_raw{n},~] = color_ge(I_raw{n},n);    
    [inimg_color{n},~] = color_ge(inimg{n},n);
 end
 for n=1:num_fluo
    imwrite(I_unmix_nnmf{n},strcat('nNMF_',num2str(n),'.tif'));%save images in sequence
     imwrite(I_unmix_BINGO{n},strcat('BINGO_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_unmix_pinv{n},strcat('PINV_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_unmix_simi{n},strcat('SIMI_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_unmix_pinvm{n},strcat('PINVm_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_unmix_simim{n},strcat('SIMIm_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_unmix_picasso{n},strcat('PICASSO_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_unmix_ri{n},strcat('NMFRI_',num2str(n),'.tif'));%save images in sequence
    imwrite(I_raw{n},strcat('RAW_',num2str(n),'.tif'));%save images in sequence
    imwrite(inimg_color{n},strcat('GT_',num2str(n),'.tif'));%save images in sequence
 end
%  
%  

 %% save enchanced images

  imwrite(nNMFunmix_SUM,strcat('nNMF_SUM.tif'));%save images in sequence  
 imwrite(nNMFunmix_MIP,strcat('nNMF_MIP.tif'));%save images in sequence  
 imwrite(BINGOunmix_SUM,strcat('BINGO_SUM.tif'));%save images in sequence  
 imwrite(BINGOunmix_MIP,strcat('BINGO_MIP.tif'));%save images in sequence  
 imwrite(PINVunmix_SUM,strcat('PINV_SUM.tif'));%save images in sequence  
 imwrite(PINVunmix_MIP,strcat('PINV_MIP.tif'));%save images in sequence 
 imwrite(SIMIunmix_SUM,strcat('SIMI_SUM.tif'));%save images in sequence  
 imwrite(SIMIunmix_MIP,strcat('SIMI_MIP.tif'));%save images in sequence 
 imwrite(PINVmunmix_SUM,strcat('PINVm_SUM.tif'));%save images in sequence  
 imwrite(PINVmunmix_MIP,strcat('PINVm_MIP.tif'));%save images in sequence 
 imwrite(SIMImunmix_SUM,strcat('SIMIm_SUM.tif'));%save images in sequence  
 imwrite(SIMImunmix_MIP,strcat('SIMIm_MIP.tif'));%save images in sequence 
 imwrite(PICASSOunmix_SUM,strcat('PICASSO_SUM.tif'));%save images in sequence  
 imwrite(PICASSOunmix_MIP,strcat('PICASSO_MIP.tif'));%save images in sequence 
 imwrite(RIunmix_SUM,strcat('RI_SUM.tif'));%save images in sequence  
 imwrite(RIunmix_MIP,strcat('RI_MIP.tif'));%save images in sequence 
 imwrite(Ipre_SUM,strcat('RAW_SUM.tif'));%save images in sequence  
 imwrite(Ipre_MIP,strcat('RAW_MIP.tif'));%save images in sequence 
 imwrite(GT_SUM,strcat('GT_SUM.tif'));%save images in sequence  
 imwrite(GT_MIP,strcat('GT_MIP.tif'));%save images in sequence 

    end
end

save results;
cd(ori_path);

fprintf('###################################\n');
fprintf('###                             ###\n');
fprintf('###         Accomplish!!!       ###\n');
fprintf('###                             ###\n');
fprintf('###################################\n');












 