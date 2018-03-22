function out = SS_dage(h_dage)
% This function creates a serpentine spring test structure designed to be
% used with the Dage

%% Default Values
if ~isfield(h_dage,'label')
    h_dage.label = 'Hopper';
end

if ~isfield(h_dage,'w')
    h_dage.w = 200;
end

if ~isfield(h_dage,'layer')
    default_layer_properties;
    h_dage.layer = [SOI SOIHOLE];
end

if ~isfield(h_dage,'NOTDUMMY')
    default_layer_properties;
    h_dage.NOTDUMMY = NOTDUMMY;
end


if ~isfield(h_dage,'noetch')
    h_dage.noetch = 0;
end

if ~isfield(h_dage,'shuttle_dx')
    h_dage.shuttle_dx = 1000;    %Total possible deflection on main shuttle
end

if ~isfield(h_dage,'N')
    h_dage.N = 2;
end


%% 

shuttle_w = 400;
spring_gap = 20; 

bar_width = 200; 

meander_gap = 3;

dage_gap = 300;         %Gap to leave to let Dage get in



% Top foot
h_rect.x = h_dage.p0(1);
h_rect.y = h_dage.p0(2);
h_rect.w = 2*(h_dage.h_ss.dpp + spring_gap) + shuttle_w;
h_rect.l = bar_width;
h_rect.layer = h_dage.layer(1);
h_rect.rounded = 0;
h_rect.etch = 1;
h_rect.ETCH_R = 5;
h_rect.SOI_HOLE = h_dage.layer(2);
h_rect.circle_etch = 1;
foot_s = rect(h_rect);

% Dummy fill
dummy_gap = spring_gap;

h_rect.x = h_dage.p0(1) - dummy_gap;
h_rect.y = h_dage.p0(2);
h_rect.w = 2*(h_dage.h_ss.dpp + spring_gap) + shuttle_w + 2*dummy_gap;
h_rect.l = bar_width + dage_gap;
h_rect.layer = h_dage.NOTDUMMY;
h_rect.etch = 0;
df_from_foot_up = rect(h_rect);

h_rect.y = h_rect.y + 10;
h_rect.l = - h_dage.shuttle_dx - 10;

df_from_foot_down = rect(h_rect);



h_ss = h_dage.h_ss;
% left side Serpentine Springs
for i = 1:h_dage.N/2
    % anchor
    h_rect.l = bar_width;
    h_rect.x = h_dage.p0(1)+h_ss.dpp/2 - (h_dage.N/2-1)*(h_ss.dpp + spring_gap) + (i-1)*(h_ss.dpp + spring_gap) - h_ss.dpp/2;
    h_rect.y = h_dage.p0(2) - bar_width - h_dage.shuttle_dx;
    h_rect.layer = h_dage.layer(1);
    h_rect.w = h_dage.h_ss.dpp;
    h_rect.etch = 0;
    str = sprintf('anc_%d = rect(h_rect);',i);
    eval(str);
    
    % spring
    h_ss.p1 = [h_dage.p0(1)+h_ss.dpp/2 - (h_dage.N/2-1)*(h_ss.dpp + spring_gap) h_rect.y] + [(i-1)*(h_ss.dpp + spring_gap) 0];
    h_ss.p2 = h_ss.p1 - [0 2*h_ss.n*(h_ss.w+meander_gap)];
    h_ss.layer = h_dage.layer(1);
    str = sprintf('ss_%d = s_spring(h_ss);',i);
    eval(str);
end

% right side Serpentine Springs
xshift = h_dage.N/2*(h_ss.dpp + spring_gap) + shuttle_w + spring_gap;
for i = 1:h_dage.N/2
    % anchor
    h_rect.l = bar_width;
    h_rect.x = xshift + h_dage.p0(1)+h_ss.dpp/2 - (h_dage.N/2-1)*(h_ss.dpp + spring_gap) + (i-1)*(h_ss.dpp + spring_gap) - h_ss.dpp/2;
    h_rect.y = h_dage.p0(2) - bar_width - h_dage.shuttle_dx;
    h_rect.w = h_dage.h_ss.dpp;
    h_rect.etch = 0;
    str = sprintf('anc_r_%d = rect(h_rect);',i);
    eval(str);
    
    % spring
    h_ss.p1 = [xshift + h_dage.p0(1)+h_ss.dpp/2 - (h_dage.N/2-1)*(h_ss.dpp + spring_gap) h_rect.y] + [(i-1)*(h_ss.dpp + spring_gap) 0];
    h_ss.p2 = h_ss.p1 - [0 2*h_ss.n*(h_ss.w+meander_gap)];
    str = sprintf('ss_r_%d = s_spring(h_ss);',i);
    eval(str);
end

% bottom cross beam
h_rect.l = bar_width;
h_rect.x = h_dage.p0(1) - (h_dage.N/2-1)*(h_dage.h_ss.dpp + spring_gap);
h_rect.y = h_ss.p2(2) - bar_width;
h_rect.w = h_dage.N*(h_dage.h_ss.dpp + spring_gap) + shuttle_w;
h_rect.layer = h_dage.layer(1);
h_rect.etch = 1;
h_rect.SOI_HOLE = h_dage.layer(2);
h_rect.circle_etch = 1;
cross_beam_1 = rect(h_rect);

% Dummy fill
h_rect.l = 2*bar_width + 2*h_ss.n*(h_ss.w+meander_gap);
h_rect.x = h_dage.p0(1) - (h_dage.N/2-1)*(h_dage.h_ss.dpp + spring_gap)- dummy_gap;
h_rect.y = h_ss.p2(2) - bar_width;
h_rect.w = h_dage.N*(h_dage.h_ss.dpp + spring_gap) + shuttle_w + 2*dummy_gap;
h_rect.layer = h_dage.NOTDUMMY;
h_rect.etch = 0;
df_up = rect(h_rect);

h_rect.y = h_rect.y + 10; 
h_rect.l = - h_dage.shuttle_dx - 10;
df_down = rect(h_rect);

% Add in large block to decrease exposed area (and use as a staple!)

h_rect.x = h_rect.x + 100;
h_rect.y = h_rect.y - 10 - 50;
h_rect.w = h_rect.w - 200;
h_rect.l = h_rect.l + 100;
h_rect.layer = h_dage.layer(1);
h_rect.rounded = 0;
h_rect.etch = 1;
h_rect.ETCH_R = 10;
h_rect.SOI_HOLE = h_dage.layer(2);
h_rect.circle_etch = 0;
df_staple = rect(h_rect);



% Central shuttle
h_rect.x = h_dage.p0(1) + h_dage.h_ss.dpp + spring_gap; 
h_rect.y = h_ss.p2(2) - h_dage.NOTDUMMY;
h_rect.w = shuttle_w;
h_rect.l = h_dage.shuttle_dx + 2*h_ss.n*(h_ss.w+meander_gap) + bar_width + 16;
h_rect.layer = h_dage.layer(1);
h_rect.rounded = 0;
h_rect.etch = 1;
h_rect.SOI_HOLE = h_dage.layer(2);
h_rect.circle_etch = 1;
central_shuttle = rect(h_rect);



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