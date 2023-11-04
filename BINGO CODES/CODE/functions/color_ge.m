function [I_COLOR,name] = color_ge(I_gray,n)
%% define the color of each fluo
channel_num = 10;
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

%% search fluo
switch n
    case 1
        name = 'blue';
    case 2
        name = 'blue-green';
    case 3
        name = 'green';
    case 4
        name = 'yellow-green';
    case 5
        name = 'yellow';
    case 6
        name = 'orange';
    case 7
        name = 'red-orange';
    case 8
        name = 'red';
    case 9
        name = 'carmine';
    case 10
        name = 'crimson';
    case 11
        name = 'scarlet';
end
 %% convert to rgb
I_a = gray2ind(I_gray,255);
I_COLOR = ind2rgb(I_a,map{n});
% imshow(I_COLOR);
% imagesc(I_gray),colormap(map);
% title(name);
% colorbar;
% colormap(map);
% axis off;

end