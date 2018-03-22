function out = align_mark(h_align)
% Function to create an alignment mark
% p0 = central point of cross
% theta = angle of rotation of cross
% w = width of cross

if ~isfield(h_align,'theta')
    h_align.theta = 0;
end


if ~isfield(h_align,'num_align')
    h_align.num_align = 2;
end

h_align.l = 100;
h_align.w = [5 10 5 10 5 10 5 10 5 10];
layers = [h_align.layer 5];



for k = 1:length(layers)
    for i=1:h_align.num_align
        h_rect.x = h_align.p0(1) - h_align.l/2 + 1.5*(i-1)*h_align.l;
        h_rect.y = h_align.p0(2) - h_align.w(i)/2;
        h_rect.w = h_align.l;
        h_rect.l = h_align.w(i);
        h_rect.layer = layers(k);
        h_rect.theta = h_align.theta;
        h_rect.p0 = h_align.p0;
        str = sprintf('vertical_%d_%d = rect(h_rect);',i,k);
        eval(str)
        
        h_rect.x = h_align.p0(1) - h_align.w(i)/2+ 1.5*(i-1)*h_align.l;
        h_rect.y = h_align.p0(2) - h_align.l/2;
        h_rect.w = h_align.w(i);
        h_rect.l = h_align.l;
        h_rect.layer = layers(k);
        h_rect.theta = h_align.theta;
        h_rect.p0 = h_align.p0;
        str = sprintf('horiz_%d_%d = rect(h_rect);',i,k);
        eval(str)
    end
end

str_name = sprintf('AM_[%d,%d]_%d',round(h_align.p0(1)),round(h_align.p0(2)),h_align.layer);

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