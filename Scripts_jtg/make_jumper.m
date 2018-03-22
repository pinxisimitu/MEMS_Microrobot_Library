function out=make_jumper(h_jumper)
% Function to create a jumping microrobot
% p0 = bottom left coordinate of the hammer
% w = width (in x direction)
% l = length (in y direction)
% layer = GDS layer to place circle into
% xnum = number of repeats in x direction
% xspace = gap between x repeats
% ynum = number of repeats in y direction
% yspace = gap between y repeats
% theta = angle by which to rotate rectangle
% p0 = point about which to rotate
% rp = return points, if rp==1, returns column vector of verticies

h_hammer.p0 = h_jumper.p0;
h_hammer.num_springs = h_jumper.num_springs;
h_hammer.theta_arm = -90;
h_hammer.chiplet = 2;
h_hammer.alumina = 1;
h_hammer.stages = 4;
h_hammer.no_backstops = 1;      %Removes unnecessary backstops
h_hammer.compact_latch = 1;
h_hammer.h_latch.ss.n = 10;
h_hammer.zif = 3;
h_hammer.spring_l = h_jumper.spring_l;
h_hammer.displacement = h_jumper.dx;            %Microns of travel for the hammer
h_hammer.hhead_r = 10;                          %No hammer head
h_hammer.no_head = 1;
h_hammer.closed = h_jumper.closed;
h_hammer.h_latch.ss.dpp = 85;
h_hammer.h_latch.ss.dist_from_rotor = 300;
h_hammer.h_latch_r = h_jumper.r;
h_hammer.r_stages = h_jumper.r_stages;
hamo = make_hammer(h_hammer);

h_hammer.rp = 1;
out = make_hammer(h_hammer);


%Add rectangular extension with etch holes with length dx
h_hbody.x = out.body_tl(1);
h_hbody.y = out.body_tl(2);
h_hbody.w = out.body_w;
h_hbody.l = h_jumper.dx;
h_hbody.layer = 6;
ham_extension = rect(h_hbody);

h_etch.regions = cell(1,1);
h_etch.r = 2;
h_etch.undercut = 6;
section.p0 = [h_hbody.x, h_hbody.y - h_etch.r];
section.type = 'rect';
section.w = h_hbody.w;
section.l = h_hbody.l + h_etch.r;
h_etch.regions{1,1} = section;
extension_holes = etch_hole(h_etch);

%Add rectangular cross bar with etch holes
cross_bar_w = h_jumper.cross_bar_w;
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2;
h_hbody.y = out.body_tl(2) + h_jumper.dx;
h_hbody.w = cross_bar_w;
h_hbody.l = h_jumper.dx;
h_hbody.layer = 6;
cross_bar = rect(h_hbody);

h_etch.regions = cell(1,1);
h_etch.r = 2;
h_etch.undercut = 6;
section.p0 = [h_hbody.x, h_hbody.y - h_etch.r];
section.type = 'rect';
section.w = h_hbody.w;
section.l = h_hbody.l + h_etch.r;
h_etch.regions{1,1} = section;
cross_bar_holes = etch_hole(h_etch);

%Add pad on left that backside Si will attach to
pad_w = 250;
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = pad_w;
h_hbody.l = pad_w;
h_hbody.layer = 6;
left_pad = rect(h_hbody);


%Adds layer 4 to re-trench
extra = 50;
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2-extra;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = pad_w+2*extra;
h_hbody.l = pad_w+extra;
h_hbody.layer = 4;
left_pad_backside = rect(h_hbody);


%Adds trench around the pad on the left
extra = extra + 150;
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2-extra;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = pad_w+2*extra;
h_hbody.l = pad_w+extra;
h_hbody.layer = 5;
left_pad_retrench= rect(h_hbody);
h_hbody.layer = 8;
left_pad_retrench_2 = rect(h_hbody);

%Subtract area right avobe crossbar so it's not bolted to substrate 
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = cross_bar_w;
h_hbody.l = 50;
h_hbody.layer = 8;
not_trench_unbolt = rect(h_hbody);






%Add pad on right that backside Si will attach to
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2+cross_bar_w - pad_w;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = pad_w;
h_hbody.l = pad_w;
h_hbody.layer = 6;
right_pad = rect(h_hbody);


%Adds layer 4 to re-trench
extra = 50;
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2+cross_bar_w - pad_w - extra;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = pad_w+2*extra;
h_hbody.l = pad_w+extra;
h_hbody.layer = 4;
right_pad_backside = rect(h_hbody);

%Adds trench around the pad
extra = extra + 150;
h_hbody.x = out.body_tl(1) - (cross_bar_w-out.body_w)/2+cross_bar_w - pad_w - extra;
h_hbody.y = out.body_tl(2) + 2*h_jumper.dx;
h_hbody.w = pad_w+2*extra;
h_hbody.l = pad_w+extra;
h_hbody.layer = 5;
right_pad_retrench = rect(h_hbody);
h_hbody.layer = 8;
right_pad_retrench_2 = rect(h_hbody);


%% Add label to jumper


h_label.p0 = h_hammer.p0 - [200 0];
h_label.text = h_jumper.label;
h_label.height = 30;
name_label = add_label(h_label);


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