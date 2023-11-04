function I_unmix_order = subb(I_unmix_order)
 I_unmix_order{1} = imadjust(I_unmix_order{1},[0.1,1],[]);
 I_unmix_order{2} = imadjust(I_unmix_order{2},[0.1,1],[]);
  I_unmix_order{3} = imadjust(I_unmix_order{3},[0.1,1],[]);
 I_unmix_order{5} = imadjust(I_unmix_order{5},[0.1,1],[]);
 I_unmix_order{6} = imadjust(I_unmix_order{6},[0.1,1],[]);
 I_unmix_order{7} = imadjust(I_unmix_order{7},[0.1,1],[]);
 I_unmix_order{8} = imadjust(I_unmix_order{8},[0.1,1],[]);
  I_unmix_order{4} = medfilt2(I_unmix_order{4}-I_unmix_order{2}-I_unmix_order{5}-I_unmix_order{3}-I_unmix_order{6}-I_unmix_order{7}-I_unmix_order{8});
 I_unmix_order{5} = medfilt2(I_unmix_order{5}-I_unmix_order{6}-I_unmix_order{3});
  I_unmix_order{6} = medfilt2(I_unmix_order{6}-I_unmix_order{7});
 I_unmix_order{8} = medfilt2(I_unmix_order{8}-I_unmix_order{7});
 I_unmix_order{9} = medfilt2(I_unmix_order{9}-I_unmix_order{7});
end