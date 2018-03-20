function gca_struct = make_comb_array(h_ca)
% Function to make a comb array


if ~isfield(h_ca,'rotation_theta')
    h_ca.rotation_theta = 0;
end

if ~isfield(h_ca,'rotation_center')
    h_ca.rotation_center = [0 0];
end


if ~isfield(h_ca,'VREP_group')
    h_ca.VREP_group = 'none';
end


% gs=gds_structure(['comb_arr_' num2str(round(h_ca.pos(1))) '_' num2str(round(h_ca.pos(2)))]);
x=h_ca.pos(1);
y=h_ca.pos(2);
length_tot=h_ca.length+h_ca.length_support;
% xy1=[x,y;x,y+h_ca.width;x+length_tot,y+h_ca.width;x+length_tot,y;x,y];
% 
% rot.pts = xy1;
% rot.theta = h_ca.rotation_theta;
% rot.p0 = h_ca.rotation_center;
% xy1=rotate_pts(rot);
%             
%             
% gs(end+1)=gds_element('boundary', 'xy',xy1, 'layer',h_ca.layer);

% Initial (bottom most) stator arm
% h_rect.x = x;
% h_rect.y = y;
% h_rect.w = length_tot;
% h_rect.l = h_ca.width;
% h_rect.p0 = h_ca.rotation_center;
% h_rect.theta = h_ca.rotation_theta;
% h_rect.layer = h_ca.layer;
% init = rect(h_rect);

if h_ca.side==0 % left side array
    h_rect.x = x;
    h_rect.y = y;
    h_rect.w = length_tot;
    h_rect.l = h_ca.width;
    h_rect.p0 = h_ca.rotation_center;
    h_rect.ynum = floor((h_ca.N)/2);
    h_rect.yspace = h_ca.gap1 + h_ca.gap2 + h_ca.width;
    h_rect.theta = h_ca.rotation_theta;
    h_rect.VREP_group = [h_ca.VREP_group  '_a'];
    h_rect.layer = h_ca.layer;
    left_side_fingers = rect(h_rect);
    
    h_rect.x = x+h_ca.length_support;
    h_rect.y = y+h_ca.width+h_ca.gap1;
    h_rect.w = length_tot;
    h_rect.l = h_ca.width;
    h_rect.p0 = h_ca.rotation_center;
    h_rect.ynum = floor((h_ca.N)/2);
    h_rect.VREP_group = [h_ca.VREP_group  '_r'];
    h_rect.yspace = h_ca.gap1 + h_ca.gap2 + h_ca.width;
    h_rect.theta = h_ca.rotation_theta;
    h_rect.layer = h_ca.layer;
    left_side_fingers_2 = rect(h_rect);

else              % right side array
    h_rect.x = x;
    h_rect.y = y;
    h_rect.w = length_tot;
    h_rect.l = h_ca.width;
    h_rect.p0 = h_ca.rotation_center;
    h_rect.ynum = floor((h_ca.N)/2);
    h_rect.VREP_group = [h_ca.VREP_group  '_a'];
    h_rect.yspace = h_ca.gap1 + h_ca.gap2 + h_ca.width;
    h_rect.theta = h_ca.rotation_theta;
    h_rect.layer = h_ca.layer;
    right_side_fingers = rect(h_rect);
   
   
    h_rect.x = x-h_ca.length_support;
    h_rect.y = y+h_ca.width+h_ca.gap1;
    h_rect.w = length_tot;
    h_rect.l = h_ca.width;
    h_rect.p0 = h_ca.rotation_center;
    h_rect.ynum = floor((h_ca.N)/2);
    h_rect.VREP_group = [h_ca.VREP_group  '_r'];
    h_rect.yspace = h_ca.gap1 + h_ca.gap2 + h_ca.width;
    h_rect.theta = h_ca.rotation_theta;
    h_rect.layer = h_ca.layer;
    right_side_fingers2 = rect(h_rect);
    
end


% for i=1:h_ca.N-1
%     if h_ca.side==0
%         if mod(i,2)==0
%             y=y+h_ca.width+h_ca.gap2;
%             x=x-h_ca.length_support;
%             xy1=[x,y;x,y+h_ca.width;x+length_tot,y+h_ca.width;x+length_tot,y;x,y];
%             
%             rot.pts = xy1;
%             rot.theta = h_ca.rotation_theta;
%             rot.p0 = h_ca.rotation_center;
%             xy1=rotate_pts(rot);
%             
%             gs(end+1) = gds_element('boundary', 'xy',xy1, 'layer',h_ca.layer);
%         else
%             y=y+h_ca.width+h_ca.gap1;
%             x=x+h_ca.length_support;
%             xy1=[x,y;x,y+h_ca.width;x+length_tot,y+h_ca.width;x+length_tot,y;x,y];
%             rot.pts = xy1;
%             rot.theta = h_ca.rotation_theta;
%             rot.p0 = h_ca.rotation_center;
%             xy1=rotate_pts(rot);
%             
%             gs(end+1) = gds_element('boundary', 'xy',xy1, 'layer',h_ca.layer);
%         end
%     elseif h_ca.side==1
%         if mod(i,2)==0
%             y=y+h_ca.width+h_ca.gap2;
%             x=x+h_ca.length_support;
%             xy1=[x,y;x,y+h_ca.width;x+length_tot,y+h_ca.width;x+length_tot,y;x,y];
%             rot.pts = xy1;
%             rot.theta = h_ca.rotation_theta;
%             rot.p0 = h_ca.rotation_center;
%             xy1=rotate_pts(rot);
%             
%             gs(end+1) = gds_element('boundary', 'xy',xy1, 'layer',h_ca.layer);
%         else
%             y=y+h_ca.width+h_ca.gap1;
%             x=x-h_ca.length_support;
%             xy1=[x,y;x,y+h_ca.width;x+length_tot,y+h_ca.width;x+length_tot,y;x,y];
%             rot.pts = xy1;
%             rot.theta = h_ca.rotation_theta;
%             rot.p0 = h_ca.rotation_center;
%             xy1=rotate_pts(rot);
%             
%             gs(end+1) = gds_element('boundary', 'xy',xy1, 'layer',h_ca.layer);
%         end
%     end
% end

% clear gs
%% Grab all the GDS structures and arrays of structures

a=whos();
b={};
c = 0;
for i=1:length(a)
    if(strcmp(a(i).class,'gds_structure'))
        c = c+1;
        str = sprintf('b{c} = %s;',a(i).name);
        eval(str);
    elseif(strcmp(a(i).class,'cell'))
        str = sprintf('temp = %s;',a(i).name);
        eval(str);
        if(isempty(temp))
            fprintf('Empty Cell! Something went wrong with %s!!\n',a(i).name)
            break;
        end
        str = sprintf('strcmp(class(%s{1}),''gds_structure'');',a(i).name);
        if(eval(str))
            str = sprintf('temp = %s;',a(i).name);
            eval(str)
            for i=1:length(temp)
                c = c+1;
                b{c} = temp{i};
            end
        end
    end
end

% Outputs a cell array of GDS Structures
gca_struct = b;







