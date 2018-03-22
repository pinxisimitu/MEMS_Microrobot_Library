function out = gim(h_gim)
% This function creates a greater inchworm motor. It uses two latching
% mechanisms to pull down a MEMS hammer body with long springs. 

if ~isfield(h_gim,'pin_gap')
    h_gim.pin_gap = 3;
end

if ~isfield(h_gim,'beak_l')
    h_gim.beak_l = 22;
end

if ~isfield(h_gim,'ss_contact_w')
    h_gim.ss_contact_w = 7;
end


if ~isfield(h_gim,'r')
    h_gim.r = [50 100];
end




%% Create the body of the hammer
% Make a MEMS Hammer
h_hammer.p0 = h_gim.p_init;            % Coordinate of bottom left of hammer
h_hammer.num_springs = h_gim.num_springs;      % Number of springs on hammer
h_hammer.spring_l = h_gim.spring_l;       % Total length of the spring (two sides added together)
h_hammer.stages = 0;
h_hammer.no_bottom = 1;
h_hammer.no_head = 1;
h_hammer.ham_w = 400;
h_hammer.anchor_w = 100;
h_hammer.displacement = 100;


h_hammer.closed = h_gim.closed;
h_hammer.noetch = h_gim.noetch;
g_hammer=make_hammer(h_hammer); % This is the function that generates the hammer gds_strucutre

%Points to make the dummy fill for the GIM
dfill_pts= zeros(4,2);
dfill_pts(1,:) = h_hammer.p0 + [h_hammer.anchor_w -2*h_hammer.displacement];
dfill_pts(2,:) = h_hammer.p0 + [h_hammer.anchor_w+h_hammer.spring_l+h_hammer.ham_w -2*h_hammer.displacement];

%% Create the Shuttle for GIM

h_gim.hw = 21;                  %Rough width of the hold (on the shuttle)

h_gim.l = 2*(h_gim.n+1)*h_gim.dx;

h_gim.p0 = h_hammer.p0 + [h_hammer.anchor_w+h_hammer.spring_l/2+h_hammer.ham_w/2-h_gim.w/2 -h_gim.l];

p1 = h_gim.p0 + [h_gim.w 0];
p2 = p1 + [0 h_gim.l];
p3 = p2 - [h_gim.w 0];

increase_dy_factor = 1.3; 

%Creates the cutouts on the right side
temp_points_r = [];
for i=1:h_gim.n
    p4 = p1 + [0 h_gim.dx] + 2*(i-1)*[0 h_gim.dx];
    np_p5 = 10;              %Number of points to put in p5
    p5_x = linspace(p4(1),p4(1)-h_gim.hw,np_p5);
    p5_y = p4(2)*ones(1,np_p5);
    p7 =  p4 + [0 increase_dy_factor*h_gim.dx];
    p6 = p7 - [h_gim.hw/2 h_gim.dx/3];
    init_points = [p4',[p5_x;p5_y],p6',p7'];
    
    interm_pts = fnplt(cscvn(init_points));
    
    final_pts = [interm_pts(1,1) downsample(interm_pts(1,2:end-1),2) interm_pts(1,end);
        interm_pts(2,1) downsample(interm_pts(2,2:end-1),2) interm_pts(2,end)];
    
    temp_points_r = [temp_points_r;final_pts'];
    
    %Add etch holes to the area that was just laid out
    close_to_edge = 2;                      %Moves start of etch holes closer to edges
    section2.p0 = [h_gim.p0(1) p4(2)-close_to_edge];
    section2.type = 'rcurve';
    section2.w = h_gim.w;
    section2.l = h_gim.dx+2*close_to_edge;
    section2.rcurve = final_pts';
    h_etch.regions = cell(1,1);
    h_etch.regions = {section2};
    h_etch.undercut = 5;
    el = sprintf('eh_rs_%d = etch_hole(h_etch);',i);
    eval(el);
end

%Creates the cutouts on the left side
temp_points_l = [];
for i=1:h_gim.n
    p4 = p3 - 2*i*[0 h_gim.dx];
    np_p5 = 10;                                     %Number of points to put in p5
    p5_x = linspace(p4(1)+h_gim.hw,p4(1),np_p5);
    p5_y = p4(2)*ones(1,np_p5);
    p7 =  p4 + [0 increase_dy_factor*h_gim.dx];
    p6 = p7 + [h_gim.hw/2 -h_gim.dx/3];
    init_points = [p4',[p5_x;p5_y],p6',p7'];
    init_points = [p7',p6',[p5_x;p5_y],p4'];
    
    interm_pts = fnplt(cscvn(init_points));
    
    final_pts = [interm_pts(1,1) downsample(interm_pts(1,2:end-1),2) interm_pts(1,end);
        interm_pts(2,1) downsample(interm_pts(2,2:end-1),2) interm_pts(2,end)];
    
    temp_points_l = [temp_points_l;final_pts'];
    
    %Add etch holes to the area that was just laid out
    close_to_edge = 2;                      %Moves start of etch holes closer to edges
    section2.p0 = [h_gim.p0(1)+h_gim.w p4(2)-close_to_edge];
    section2.type = 'lcurve';
    section2.w = h_gim.w;
    section2.l = h_gim.dx+2*close_to_edge;
    section2.rcurve = final_pts';
    h_etch.regions = cell(1,1);
    h_etch.regions = {section2};
    el = sprintf('eh_ls_%d = etch_hole(h_etch);',i);
    eval(el);
end

points = [h_gim.p0;p1;temp_points_r;p2;p3;temp_points_l];

temp = gds_element('boundary', 'xy',points,'layer',6);
str_name = sprintf('GIM_S_[%d,%d]',round(h_gim.p0(1)),round(h_gim.p0(2)));
serpentine_spline_l = gds_structure(str_name,temp);

% Add etch holes to the top and bottom rectangles on the shuttle

% Add etch holes to the bottom rectangle of the shuttle
closer_to_edge = 1;                     %Pushing etch holes closer to edges.
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = h_gim.p0 - [0 closer_to_edge];
section.type = 'rect';
section.w = h_gim.w;
section.l = h_gim.dx+2*closer_to_edge;
h_etch.regions{1,1} = section;

g_shuttle_bottom = etch_hole(h_etch);

% Add etch holes to the top rectangle of the shuttle
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = h_gim.p0 + [0 2*h_gim.dx*(h_gim.n+.5)-closer_to_edge];
section.type = 'rect';
section.w = h_gim.w;
section.l = h_gim.dx+2*closer_to_edge+4;
h_etch.regions{1,1} = section;

g_shuttle_top = etch_hole(h_etch);

% Add guide to left of the shuttle
gap = 2;
guide_w = 500;
guide_l = 200;
h_rect.x = h_gim.p0(1) - gap - guide_w;
h_rect.y = h_gim.p0(2) + 250;
h_rect.w = guide_w;
h_rect.l = guide_l;
h_rect.layer = 6;
h_rect.rounded = 20;
left_GIM_guide = rect(h_rect);

% Add guide to right of the shuttle
h_rect.x = p1(1) + gap;
h_rect.y = p1(2) + 250;
h_rect.w = guide_w;
h_rect.l = guide_l;
h_rect.layer = 6;
h_rect.rounded = 20;
right_GIM_guide = rect(h_rect);



%% Add rotors and lever arms

%Right side lever arm
h_latch.r = h_gim.r;
%h_latch.p0 = h_gim.p0 + [h_gim.w+h_latch.r(2)+2 0] + [0 h_gim.dx]; 
h_latch.p0 = h_gim.p0 + [h_gim.w+h_latch.r(2)+2 0];
dfill_pts(3,:) = h_latch.p0;

%Beak length
h_latch.blength = h_gim.beak_l;

h_latch.orientation = 1;

h_latch.h_joint = 0;
h_latch.inchworm = 0;

h_latch.ss_contact_w = h_gim.ss_contact_w;

h_latch.ss.n = 7;
h_latch.ss.dpp = 70;                          %was 70 in the designs that worked
h_latch.ss.dist_from_rotor = 200;

h_latch.bbeak_fillet_length = 5;
h_latch.n = 100;
h_latch.arml =900;
h_latch.armw = 2*h_latch.r(1);

h_latch.hhead_r = 50;
h_latch.closed = h_gim.closed;
h_latch.alumina = 0;

h_latch.actuation_angle = 40;
h_latch.init_angle = 40;

h_latch.theta = h_latch.init_angle;
h_latch.theta_arm = -90;
h_latch.cstage = 1;

h_latch.chiplet = 0;
h_latch.layer = 6;
h_latch.noetch = h_hammer.noetch;
h_latch.stages = 1;
h_latch.mech_latch = 0;
h_latch.compact_latch = 0;
h_latch.no_backstops = 0;
h_latch.pin_gap = h_gim.pin_gap;
[g_latch latch_ref_pts]= ham_latch(h_latch);

%Left side lever arm
h_latch.p0 = h_gim.p0 - [h_latch.r(2)+2 0];
dfill_pts(4,:) = h_latch.p0;
h_latch.orientation = -1;

h_latch.actuation_angle = 40;
h_latch.init_angle = 140;
h_latch.theta = h_latch.init_angle;
h_latch.theta_arm = -90;
h_latch.gim = 1;

[g_ham_latch latch_ref_pts]= ham_latch(h_latch);

% Add fillets to the right side of the shuttle
h_fillet.d = .7;
h_fillet.layer = 6;
h_fillet.x_dist = (h_hammer.ham_w - h_gim.w)/2;
%h_fillet.y_dist = -2*(h_gim.dx+.5*(1-increase_dy_factor)*h_gim.dx);
h_fillet.y_dist = -2*h_gim.dx;
h_fillet.p0 = p2; 

h_fillet.p1 = p2 + [(h_hammer.ham_w - h_gim.w)/2 0];
h_fillet.p2 = p2 + [0 -2*(h_gim.dx+.5*(1-increase_dy_factor)*h_gim.dx)];
h_fillet.downsample = 5;

%For rotating the fillet, (just copied and pasted, will need editing to
%actually rotate these points
%h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%h_rot.theta =  -(latch_arm_initial_angle + pi/2);
%h_rot.p0 = h_fillet.p0;
%out=rotate_pts(h_rot);

%h_fillet.p0 = out(1,:);
%h_fillet.p1 = out(2,:);
%h_fillet.p2 = out(3,:);
right_fillet = fillet(h_fillet);
h_fillet.rp = 1;
rfill_pts = fillet(h_fillet);
h_fillet.rp = 0;


%Add etch holes to the right fillet
close_to_edge_x = 6;                      %Moves start of etch holes closer to edges
close_to_edge_y = 5;                      %Moves start of etch holes closer to edges
section2.p0 = h_fillet.p2 - [close_to_edge_x close_to_edge_y];
section2.type = 'rcurve';
section2.w = (h_hammer.ham_w - h_gim.w)/2 + 2*close_to_edge_x;
section2.l = 2*h_gim.dx + 2*close_to_edge_y + (1-increase_dy_factor)*h_gim.dx;
section2.rcurve = rfill_pts';
h_etch.regions = cell(1,1);
h_etch.regions = {section2};
h_etch.undercut = 4;
g_etch_r_fill = etch_hole(h_etch);

% Add fillets to the left of the shuttle
h_fillet.d = .7;
h_fillet.layer = 6;
h_fillet.x_dist = -(h_hammer.ham_w - h_gim.w)/2;
h_fillet.y_dist = -2*h_gim.dx;
h_fillet.p0 = p3; 

h_fillet.p1 = p3 - [(h_hammer.ham_w - h_gim.w)/2 0];
h_fillet.p2 = p3 - [0 (h_gim.dx+(1-increase_dy_factor)*h_gim.dx)];
h_fillet.downsample = 5;

%For rotating the fillet, (just copied and pasted, will need editing to
%actually rotate these points
%h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%h_rot.theta =  -(latch_arm_initial_angle + pi/2);
%h_rot.p0 = h_fillet.p0;
%out=rotate_pts(h_rot);

%h_fillet.p0 = out(1,:);
%h_fillet.p1 = out(2,:);
%h_fillet.p2 = out(3,:);
left_fillet = fillet(h_fillet);
h_fillet.rp = 1;
lfill_pts = fillet(h_fillet);
h_fillet.rp = 0;


%Add etch holes to the left fillet
close_to_edge_x = -6;                      %Moves start of etch holes closer to edges
close_to_edge_y = 5;                      %Moves start of etch holes closer to edges
section2.p0 = h_fillet.p2 - [close_to_edge_x close_to_edge_y];
section2.type = 'lcurve';
section2.w = (h_hammer.ham_w - h_gim.w)/2 + 2*close_to_edge_x;
section2.l = h_gim.dx + 2*close_to_edge_y + (1-increase_dy_factor)*h_gim.dx;
section2.rcurve = lfill_pts';
h_etch.regions = cell(1,1);
h_etch.regions = {section2};
h_etch.undercut = 4;
g_etch_l_fill = etch_hole(h_etch);

%% Add dummy fill

%Add dummy fill for the shuttle
extra = 50;
h_df.x = h_gim.p0(1)-extra;
h_df.y = h_gim.p0(2)-extra;
h_df.w = h_gim.w+2*extra;
h_df.l = h_gim.l+2*extra;
h_df.layer = 8;
hammer_springs = rect(h_df);

%Dummy fill for dage
extra = 150;
h_df.x = h_gim.p0(1)-extra;
h_df.y = h_gim.p0(2)-extra;
h_df.w = h_gim.w+2*extra;
h_df.l = h_gim.l+2*extra + 100;
h_df.layer = 8;
dage_DF = rect(h_df);



be = gds_element('boundary', 'xy',dfill_pts,'layer',8);
str_name = sprintf('dummff_[%d,%d]%d',round(dfill_pts(1,1)),round(dfill_pts(1,2)),8);
df_trap = gds_structure(str_name,be);


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