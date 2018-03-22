function [l_min r_max t_max b_min] = find_bb(structure)
% This function takes in a cell array containing GDS Structures and outputs
% the bounding box containing the entire strucutre.

l_min = 25000;
r_max = -25000;
t_max = -25000;
b_min = 25000;




for i=1:length(structure)
    temp = structure{i};
    %temp2 = temp.el;
    temp2 = structure{i}.el;
    temp3 = temp2{1}.xy;
    if temp2{1}.layer ~=6   %Skip the structure if the layer is not SOI
        continue;
    end
    temp4 = temp3{1};
    
    l_min = min([l_min;temp4(:,1)]);
    r_max = max([r_max;temp4(:,1)]);
    
    b_min = min([b_min;temp4(:,2)]);
    t_max = max([t_max;temp4(:,2)]);
end
