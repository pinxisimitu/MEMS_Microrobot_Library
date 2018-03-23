function out = piezoR(h_piezoR)
% This function creates a greater inchworm motor. It uses two latching
% mechanisms to pull down a MEMS hammer body with long springs. 

if ~isfield(h_piezoR,'j')
    h_piezoR.j = 3;
end

% Make hammer body
h_hammer.num_springs = h_piezoR.num_springs;
h_hammer.theta_arm = -90;
h_hammer.displacement = 40;     %Microns of travel for the hammer
h_hammer.alumina = 0;
h_hammer.no_backstops = 1;
h_hammer.h_latch.ss.n = 10;
h_hammer.h_latch.ss.dpp = 85;
h_hammer.h_latch.ss.dist_from_rotor = 300;
h_hammer.large_dx = 0;
h_hammer.chiplet = 0;
h_hammer.no_head = 1;
h_hammer.p0 = h_piezoR.p0;
h_hammer.stages = 0;
h_hammer.spring_l = h_piezoR.spring_l;
h_hammer.no_bottom = 1;
hamout = make_hammer(h_hammer);

pad = 100;
gap = 50;

%Add more not dummy to the top
h_rect.x = h_piezoR.p0(1)-pad+50-gap;
h_rect.y = h_piezoR.p0(2)-gap;
h_rect.w = h_piezoR.spring_l + 400 + 2*pad + 2*gap;
h_rect.l = pad+2*gap;
h_rect.layer = 8;
dummy_add = rect(h_rect);

%Add more not dummy to the bottom
dummy2_len = 200;
h_rect.x = h_piezoR.p0(1)-pad+50-gap;
h_rect.y = h_piezoR.p0(2)-dummy2_len;
h_rect.w = h_piezoR.spring_l + 400 + 2*pad + 2*gap;
h_rect.l = dummy2_len;
h_rect.layer = 8;
dummy_add2 = rect(h_rect);

%Add pad on left and right to contact and measure
h_rect.x = h_piezoR.p0(1)-pad+50;
h_rect.y = h_piezoR.p0(2);
h_rect.w = pad;
h_rect.l = pad;
h_rect.layer = 6;
h_rect.xnum = 2;
h_rect.xspace = h_piezoR.spring_l + 400;
pads = rect(h_rect);


%Add path for handle
handle_dy = 100;
handle_dx = 100;
w = 6;
layer = 6;

p0 = h_piezoR.p0 + [50+h_piezoR.spring_l/2+200-handle_dx/2 1];
p1 = h_piezoR.p0 + [50+h_piezoR.spring_l/2+200-handle_dx/2 0];
p2 = h_piezoR.p0 + [50+h_piezoR.spring_l/2+200 -handle_dy];
p3 = h_piezoR.p0 + [50+h_piezoR.spring_l/2+200+handle_dx/2 0];
p4 = h_piezoR.p0 + [50+h_piezoR.spring_l/2+200+handle_dx/2 1];




be = gds_element('path', 'xy',[p0;p1;p2;p3;p4],'width', w,'layer',layer);
str_name = sprintf('pulltab_[%d,%d]',round(p2(1)),round(p2(2)));
pull_tab = gds_structure(str_name,be);






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