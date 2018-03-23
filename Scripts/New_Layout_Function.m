function out = New_Layout_Function(h_nlf)
% Boilerplate to make a new function

%Default values
if ~isfield(h_nlf,'rect_size')
    h_nlf.rect_size = 10;
end

%Add a rectangle
h_rect.x = h_nlf.p0(1);
h_rect.y = h_nlf.p0(2);
h_rect.w = h_nlf.rect_size;
h_rect.l = h_nlf.rect_size;
h_rect.layer = 8;
rect_1 = rect(h_rect);


%Add a path 
p0 = [0 0];
p1 = [10 0];
p2 = [10 10];
p3 = [0 10];
p4 = [0 0];
width = 3;
layer = 8;

be = gds_element('path', 'xy',[p0;p1;p2;p3;p4],'width', width,'layer',layer);
str_name = sprintf('path_[%d,%d]',round(p2(1)),round(p2(2)));
rand_path = gds_structure(str_name,be);


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

% Outputs a cell array of GDS Structures
out = b;