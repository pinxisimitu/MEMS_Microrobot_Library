function out = m_routing(h_routing)
% Function to aid in routing. Draws both routing lines and dummy exclude
% lines

%Default values
if ~isfield(h_routing,'smooth')
    h_routing.smooth = 0;
end

if ~isfield(h_routing,'layer_dummy')
    default_layer_properties;
    h_routing.layer_dummy = NOTDUMMY;
end

if ~isfield(h_routing,'w')  % Total width of routing (including dummy exclude)
    h_routing.w = 50;
end

if ~isfield(h_routing,'gap')   % gap between routing line and dummy fill
    h_routing.gap = 5;
end

% Check if path should to be smoothed
if h_routing.smooth
    interm_pts = fnplt(cscvn(h_routing.pts'));
    h_routing.pts = interm_pts';  %To feed these points back into the path function, you will need to transpose again
end

% Create actual path
be = gds_element('path', 'xy',h_routing.pts,'width', h_routing.w-2*h_routing.gap,'layer',h_routing.layer);
str_name = sprintf('path_[%d,%d]_%d',round(h_routing.pts(2,1)),round(h_routing.pts(2,2)),h_routing.layer);
rand_path = gds_structure(str_name,be);


% Create dummy exclude around the path
be = gds_element('path', 'xy',h_routing.pts,'width', h_routing.w,'layer',h_routing.layer_dummy);
str_name = sprintf('path_[%d,%d]_%d',round(h_routing.pts(2,1)),round(h_routing.pts(2,2)),h_routing.layer_dummy);
rand_path_exclude = gds_structure(str_name,be);

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