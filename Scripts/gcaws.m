function out = gcaws(h_gca)
% Function to create the footprint of a gcaws mask
% p0 = bottom left coordinate of gcaws mask
% layer = layer to create the mask in

if ~isfield(h_gca,'layer')  % Place template into layer 6 by default
    default_layer_properties;
    h_gca.layer = [SOI SOIHOLE];
end             


h_rect.x = h_gca.p0(1);
h_rect.y = h_gca.p0(2);
h_rect.w = 2500;
h_rect.l = 25000;
h_rect.xnum = 2;
h_rect.xspace = 20000;
h_rect.layer = h_gca.layer(1);
left_and_right = rect(h_rect);

h_rect = [];
%Create top and bottom on the mask
h_rect.x = h_gca.p0(1);
h_rect.y = h_gca.p0(2);
h_rect.w = 25000;
h_rect.l = 2500;
h_rect.ynum = 2;
h_rect.yspace = 20000;
h_rect.layer = h_gca.layer(1);
t_and_b = rect(h_rect);

h_rect = [];

%Create bottom cutout which is offset
h_rect.x = h_gca.p0(1)+14500;
h_rect.y = h_gca.p0(2)+800;
h_rect.w = 400;
h_rect.l = 400;
h_rect.layer = h_gca.layer(2);
c_off = rect(h_rect);

%Create top cutout which is not offset
h_rect.x = h_gca.p0(1)+12300;
h_rect.y = h_gca.p0(2)+23800;
h_rect.w = 400;
h_rect.l = 400;
h_rect.layer = h_gca.layer(2);
c_on = rect(h_rect);

%Creates the triangular cutouts for the corners of the mask

x1 = h_gca.p0(1) + 2500;
x2 = x1 + 1850;
x3 = h_gca.p0(1) + 25000 - 2500 - 1850;
x4 = x3 + 1850;

y1 = h_gca.p0(2) + 2500;
y2 = y1 + 2780;
y3 = h_gca.p0(2) + 25000 - 2500 - 2780;
y4 =  y3 + 2780;
	

p1 = [x1 y1];
p2 = [x2 y1];
p3 = [x1 y2];
ele = gds_element('boundary', 'xy',[p1;p2;p3],'layer',h_gca.layer(1));
str_name = sprintf('gca_[%d,%d]',round(p1(1)),round(p1(2)));
gc_1 = gds_structure(str_name,ele);
    
p1 = [x3 y1];
p2 = [x4 y1];
p3 = [x4 y2];
ele = gds_element('boundary', 'xy',[p1;p2;p3],'layer',h_gca.layer(1));
str_name = sprintf('gca_[%d,%d]',round(p1(1)),round(p1(2)));
gc_2 = gds_structure(str_name,ele);


p1 = [x1 y4];
p2 = [x2 y4];
p3 = [x1 y3];
ele = gds_element('boundary', 'xy',[p1;p2;p3],'layer',h_gca.layer(1));
str_name = sprintf('gca_[%d,%d]',round(p1(1)),round(p1(2)));
gc_3 = gds_structure(str_name,ele);
    
p1 = [x3 y4];
p2 = [x4 y4];
p3 = [x4 y3];
ele = gds_element('boundary', 'xy',[p1;p2;p3],'layer',h_gca.layer(1));
str_name = sprintf('gca_[%d,%d]',round(p1(1)),round(p1(2)));
gc_4 = gds_structure(str_name,ele);
    

%% Grab all the GDS structures and arrays of structures

%Find all gds structures
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

% Outputs a cell array of
out = b;
	