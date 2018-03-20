function out = eh_window_array(h_ehwa)
% This function creates an array of etch holes with backside windows to see
% footing characteristics

if ~isfield(h_ehwa,'j')
    h_ehwa.j = 3;
end

if ~isfield(h_ehwa,'window_l')
    h_ehwa.window_l = 200;
end




init = h_ehwa.p0;
chunk_w = 20;

for i=1:10
    section.p0 = init + (i-1)*[chunk_w-4 2];
    h_etch.regions = cell(1,1);
    h_etch.r = 2;
    h_etch.undercut = 4;
    h_etch.circle_etch = 1;
    section.type = 'rect';
    section.w = chunk_w;
    section.l = h_ehwa.l;
    
    h_etch.regions{1,1} = section;
    
    str = sprintf('eh_%d = etch_hole(h_etch);',i);
    eval(str)
end



%Add a backside window into the middle of the array
h_rect.x = init(1);
h_rect.y = init(2)+h_ehwa.l/2 - h_ehwa.window_l/2;
h_rect.w = h_ehwa.w - 40;
h_rect.l = h_ehwa.window_l;
h_rect.layer = 5;
backside_wind = rect(h_rect);

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