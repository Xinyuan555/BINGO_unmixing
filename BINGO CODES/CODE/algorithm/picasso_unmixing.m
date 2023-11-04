%%PICASSO UNMIXING

function I_unmix_picasso = picasso_unmixing(outimg,pe,num_fluo)
inpica =  outimg;
if num_fluo/2==fix(num_fluo/2)
    a = num_fluo/2;
    b = 0;
else
    a = (num_fluo-1)/2;b = 1;
end

 for i = 0:a-2
     for j = 1:3
     pi_input(:,:,j) = single(im2gray(inpica{pe(j+i*2)})); 
     end
  I_unmix_picasso_mid = picasso3(pi_input);
  if i == 0
      for j = 1:3
            I_unmix_picasso{i+j}= I_unmix_picasso_mid{j};
            inpica{j} = uint8(round((I_unmix_picasso_mid{j}./max(max(I_unmix_picasso_mid{j})))*255));
      end
  else
      for j = 1:2
      I_unmix_picasso{2*(i+1)-1+j}= I_unmix_picasso_mid{j+1};
      inpica{2*(i+1)-1+j} = uint8(round((I_unmix_picasso_mid{j+1}./max(max(I_unmix_picasso_mid{j+1})))*255));
      end
  end
 end
  
     for j = 1:3
     pi_input(:,:,j) = single(im2gray(inpica{pe(num_fluo-(3-j))})); 
     end
  I_unmix_picasso_mid = picasso3(pi_input);
  if b == 0
  I_unmix_picasso{num_fluo}= I_unmix_picasso_mid{3};
  else
  I_unmix_picasso{num_fluo}= I_unmix_picasso_mid{3};
  I_unmix_picasso{num_fluo-1}= I_unmix_picasso_mid{2};
  end
end
  function I_ouput = picasso3(sourceIMG)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% 4 color Unmixing %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qN = 100; maxIter = 200; step_size = 0.2; 
[imgDemixed, unmixing_log] = PICASSO_3C(sourceIMG, qN, maxIter, step_size, 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(sourceIMG,3)
    imgDemixed(:,:,i) = imgDemixed(:,:,i)./max(max(imgDemixed(:,:,i)));
    I_ouput{i} = uint8(round(imgDemixed(:,:,i)*255));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
