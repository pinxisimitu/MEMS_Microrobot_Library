function out = m_shape(h_shape)
% Function to make an arbitrary path

%Default values
if ~isfield(h_shape,'smooth')
    h_shape.smooth = 0;
end

% Check if path should to be smoothed
if h_shape.smooth
    interm_pts = fnplt(cscvn(h_shape.pts'));
    h_shape.pts = interm_pts';  %To feed these points back into the path function, you will need to transpose again
end


be = gds_element('boundary', 'xy',h_shape.pts,'layer',h_shape.layer);
str_name = sprintf('shape_[%d,%d]_%d',round(h_shape.pts(2,1)),round(h_shape.pts(2,2)),h_shape.layer);
rand_shape = gds_structure(str_name,be);

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