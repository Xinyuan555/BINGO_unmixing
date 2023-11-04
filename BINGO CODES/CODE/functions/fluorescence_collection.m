function [signal] = fluorescence_collection(X,filter)   
   [m, n] = size(X);
   [~, v] = size(filter);
   signal = zeros(1,v);
    X(:,n) =  X(:,n)./max(max(X(:,n)));
    for j = 1:v
         s = 0;
        for i = 1:m
        
        if X(i,1)>=min(filter{1,j}) && X(i,1)<=max(filter{1,j})
        s = s+X(i,n);
        end
       
        end
         signal(j) = s;
    end
   
%     signal_level = [signal_B;signal_G;signal_R];
%     fprintf('FLUO %i:  Channel_B = %.4f,     Channel_G = %.4f,     Channel_R = %.4f\n', num, signal_B, signal_G, signal_R);
 end