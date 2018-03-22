function out = manual_jumper(h_gim)
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

if ~isfield(h_gim,'chiplet')
    h_gim.chiplet = 0;
end


if ~isfield(h_gim,'manual')
    h_gim.manual = 0;
end



tot_desired_displacement = 1400;    %On the central shuttle



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
%g_hammer=make_hammer(h_hammer); % This is the function that generates the hammer gds_strucutre

%Points to make the dummy fill for the GIM
dfill_pts= zeros(4,2);
dfill_pts(1,:) = h_hammer.p0 + [h_hammer.anchor_w -2*h_hammer.displacement];
dfill_pts(2,:) = h_hammer.p0 + [h_hammer.anchor_w+h_hammer.spring_l+h_hammer.ham_w -2*h_hammer.displacement];

%% Create the Shuttle for GIM

%h_gim.hw = 25;                  %Rough width of the hold (on the shuttle)

h_gim.l = 2*(h_gim.n+1)*h_gim.dx;

h_gim.p0 = h_hammer.p0 + [h_hammer.anchor_w+h_hammer.spring_l/2+h_hammer.ham_w/2-h_gim.w/2 -h_gim.l];

p1 = h_gim.p0 + [h_gim.w 0];
p2 = p1 + [0 h_gim.l];
p3 = p2 - [h_gim.w 0];

increase_dy_factor = h_gim.increase_factor; 

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

%For hopper foot
foot_y = max(points(:,2));
foot_w = 2000;
foot_h = 200;

%Add foot
h_rect.x = h_gim.p0(1) - foot_w/2 + h_gim.w/2;
h_rect.y = foot_y;
h_rect.w = foot_w;
h_rect.l =  foot_h;
h_rect.layer = 6;
h_rect.rounded = 0;
foot_s = rect(h_rect);

foot_tl = [h_rect.x h_rect.y+foot_h]; %Top left point of the foot





%Add foot etch holes
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = [h_rect.x h_rect.y-2*(h_etch.r)];
section.type = 'rect';
section.w = h_rect.w;
section.l = h_rect.l+4*(h_etch.r);
h_etch.regions{1,1} = section;
foot_eh = etch_hole(h_etch);

%Add dummy for the foot
h_rect.x = h_gim.p0(1) - foot_w/2 + h_gim.w/2 - 50;
h_rect.y = foot_y - tot_desired_displacement;
h_rect.w = foot_w + 100;
h_rect.l =  foot_h + tot_desired_displacement+250;
h_rect.layer = 8;
h_rect.rounded = 0;
foot_s_df = rect(h_rect);


if h_gim.manual == 0
    %Add backside foot re-trench
    edge_gap = 20;
    sub_foot = 200;
    h_rect.x = h_gim.p0(1) - foot_w/2 + h_gim.w/2 + edge_gap;
    h_rect.y = foot_y + edge_gap;
    h_rect.xnum = 2;
    h_rect.xspace = foot_w - 2*(edge_gap + sub_foot);
    h_rect.w = sub_foot;
    h_rect.l =  foot_h - 2*edge_gap;
    h_rect.layer = 4;
    h_rect.rounded = 0;
    foot_substr = rect(h_rect);
    
    h_rect.xnum = 1;
    
end




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
%on top
gap = 2;
guide_w = 200;
guide_l = 350;
h_rect.x = h_gim.p0(1) - gap - guide_w;
h_rect.y = h_gim.p0(2) + 150;
h_rect.w = guide_w;
h_rect.l = guide_l;
h_rect.layer = 6;
h_rect.rounded = 20;
left_GIM_guide = rect(h_rect);

%on bottom
gap = 2;
guide_ws = 80;
h_rect.x = h_gim.p0(1) - gap - guide_ws;
h_rect.y = h_gim.p0(2) - 500;
h_rect.w = guide_ws;
h_rect.l = guide_l;
h_rect.layer = 6;
h_rect.rounded = 20;
left_GIM_guide_bot = rect(h_rect);



% Add guide to right of the shuttle
%on top
h_rect.x = p1(1) + gap;
h_rect.y = p1(2) + 150;
h_rect.w = guide_w;
h_rect.l = guide_l;
h_rect.layer = 6;
h_rect.rounded = 20;
right_GIM_guide = rect(h_rect);

%on bottom
h_rect.x = p1(1) + gap;
h_rect.y = p1(2) - 500;
h_rect.w = guide_ws;
h_rect.l = guide_l;
h_rect.layer = 6;
h_rect.rounded = 20;
right_GIM_guide_bot = rect(h_rect);

%Add dummy fill near the guides
h_rect.x = h_gim.p0(1) - 1.5*h_gim.w;
h_rect.y = h_gim.p0(2) - 750;
h_rect.w = 4*h_gim.w;
h_rect.l = 1500;
h_rect.layer = 8;
h_rect.rounded = 0;
df_guide_area = rect(h_rect);








%% Add rotors and lever arms

%Right side lever arm
h_latch.r = h_gim.r;
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
h_latch.arml = 900;
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
[g_latch latch_ref_pts]= mj_latch(h_latch);

%Add right geared lever arm
h_arm.n = 100;
h_arm.opening_theta = 360 - h_latch.actuation_angle;
h_arm.inner_arm_theta = 180 - h_latch.actuation_angle/2;
h_arm.p0 = h_latch.p0;
h_arm.r = [h_gim.r(2)-6 1000];      %6 is the etch hole undercut
h_arm.orientation = h_latch.orientation;
h_arm.layer = 6;
h_arm.shuttle_w = 50;
h_arm.rack_teeth = h_gim.rack_teeth;
h_arm.manual = h_gim.manual;

[right_ga right_ga_points] = gear_arm(h_arm);


%Add right side serpentine spring
%Generate Serpentine Spring
h_ss.p1 = h_arm.p0 - [-150 500];         % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 800];          % Second point that the SS will span to
ss_r_p2 = h_ss.p2;
h_ss.n = 20;                          % Number of meanders 
h_ss.w = 7;                           % Width of beams 
h_ss.dpp = 450;                       % Peak to peak distance of meanders
h_ss.layer = 6;                       % Layer 
serpentine_h = s_spring(h_ss);

h_ss.rp = 1;
points = s_spring(h_ss);        %matrix that has each vertex of SS
h_ss.rp = 0;



%Spring protect line right
protection_w = 50;
protection_l = 4650;

h_rect.w = protection_w;
h_rect.l = protection_l;
h_rect.x = h_ss.p1(1)+h_ss.dpp/2+15;
h_rect.y = h_ss.p1(2)-h_rect.l+50;
h_rect.layer = 6;
h_rect.rounded = 0;
h_rect.xnum = 1;
h_rect.xspace = 10;
prot1 = rect(h_rect);

h_rect.xnum = 1;


%Add anchor
h_rect.w = 400;
h_rect.l = 100;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
r_ss_anc = rect(h_rect);

%Second copy
temp = h_ss.p1(2) - h_ss.p2(2);
h_ss.p1 = h_ss.p1 - [0 temp+50+100+tot_desired_displacement];            
h_ss.p2 = h_ss.p2 - [0 temp+50+100+tot_desired_displacement];          
serpentine_h_rs_2 = s_spring(h_ss);


right_ss_pf = h_ss.p2;

%Add anchor
h_rect.w = 400;
h_rect.l = 100;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
r_ss_anc2 = rect(h_rect);

%Reset for future use
h_ss.p1 = h_arm.p0 - [-250 500];         % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 800];          % Second point that the SS will span to






%Left side lever arm
h_latch.p0 = h_gim.p0 - [h_latch.r(2)+2 0];
dfill_pts(4,:) = h_latch.p0;
h_latch.orientation = -1;

h_latch.actuation_angle = 40;
h_latch.init_angle = 140;
h_latch.theta = h_latch.init_angle;
h_latch.theta_arm = -90;
h_latch.gim = 1;

[g_latch3 latch_ref_pts]= mj_latch(h_latch);

%Add left geared lever arm
h_arm.n = 100;
h_arm.opening_theta = 360 - h_latch.actuation_angle;
h_arm.inner_arm_theta = 0 + h_latch.actuation_angle/2;
h_arm.p0 = h_latch.p0;
h_arm.r = [h_gim.r(2)-6 1000];      %6 is the etch hole undercut
h_arm.layer = 6;
h_arm.orientation = h_latch.orientation;
h_arm.shuttle_w = 50;

[left_ga left_ga_points]= gear_arm(h_arm);



%Add left side serpentine spring
%Generate Serpentine Spring
h_ss.p1 = h_arm.p0 - [150 500];         % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 800];          % Second point that the SS will span to
ss_l_p2 = h_ss.p2;
serpentine_h_ls = s_spring(h_ss);

%left side protection spring
h_rect.w = protection_w;
h_rect.l = protection_l;
h_rect.x = h_ss.p1(1)-h_ss.dpp/2-35 - protection_w + 10;
h_rect.y = h_ss.p1(2)-h_rect.l+50;
h_rect.layer = 6;
h_rect.xnum = 1;
h_rect.xspace = 10;
h_rect.rounded = 0;
prot2 = rect(h_rect);
h_rect.xnum  = 1;


%Add dummy fill for the spring system
h_rect.w = 1400;
h_rect.l = protection_l;
h_rect.x = h_ss.p1(1)-h_ss.dpp/2-35 - protection_w + 10-50;
h_rect.y = h_ss.p1(2)-h_rect.l+50;
h_rect.layer = 8;
h_rect.rounded = 0;
prot2_df = rect(h_rect);




%Add anchor
h_rect.w = 400;
h_rect.l = 100;
h_rect.layer = 6;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
l_ss_anc2 = rect(h_rect);



%Second copy
temp = h_ss.p1(2) - h_ss.p2(2);
h_ss.p1 = h_ss.p1 - [0 temp+50+100+tot_desired_displacement];         
h_ss.p2 = h_ss.p2 - [0 temp+50+100+tot_desired_displacement];          
serpentine_h_ls_2 = s_spring(h_ss);

left_ss_pf = h_ss.p2;

h_ss.rp = 1;
points = s_spring(h_ss);        %matrix that has each vertex of SS
h_ss.rp = 0;

%Add anchor
h_rect.w = 400;
h_rect.l = 100;
h_rect.layer = 6;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
l_ss_anc1 = rect(h_rect);

%Rest vars for future use
h_ss.p1 = h_arm.p0 - [250 500];         % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 800];          % Second point that the SS will span to




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

%right_fillet = fillet(h_fillet);
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
%g_etch_r_fill = etch_hole(h_etch);

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
%left_fillet = fillet(h_fillet);
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
%g_etch_l_fill = etch_hole(h_etch);

%Extend the shuttle downwards
ext_length = h_gim.p0(2) - h_ss.p2(2); 
h_rect.x = h_gim.p0(1);
h_rect.y = h_gim.p0(2)-ext_length;
h_rect.w = h_gim.w;
h_rect.l = ext_length;
h_rect.layer = 6;
h_rect.rounded = 0;
extended_shuttle = rect(h_rect);

%Add etch holes to the shuttle extension 
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = [h_rect.x h_rect.y-2*(h_etch.r)];
section.type = 'rect';
section.w = h_rect.w;
section.l = h_rect.l+4*(h_etch.r);
h_etch.regions{1,1} = section;

SE_Eh = etch_hole(h_etch);

%Put top T onto shuttle to attach SS
T_width = 50;
h_rect.x = ss_l_p2(1)-h_ss.w/2;
h_rect.y = ss_l_p2(2)-T_width;
h_rect.w = ss_r_p2(1) - ss_l_p2(1)+h_ss.w;
h_rect.l = T_width;
h_rect.layer = 6;
h_rect.rounded = 0;
extended_shuttle_t = rect(h_rect);

%Add etch holes to the top T
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = [h_rect.x h_rect.y];
section.type = 'rect';
section.w = h_rect.w;
section.l = h_rect.l;
h_etch.regions{1,1} = section;

T_Eh = etch_hole(h_etch);


%Add bottom T to shuttle


h_rect.x = left_ss_pf(1) - h_ss.w/2;
h_rect.y = left_ss_pf(2)-T_width;
h_rect.w = right_ss_pf(1) - left_ss_pf(1) + h_ss.w;
h_rect.l = T_width;
h_rect.layer = 6;
h_rect.rounded = 0;
extended_shuttle_t_2 = rect(h_rect);

%Add etch holes to the top T
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = [h_rect.x h_rect.y];
section.type = 'rect';
section.w = h_rect.w;
section.l = h_rect.l;
h_etch.regions{1,1} = section;

T_Eh_2 = etch_hole(h_etch);

%Second Extension of the shuttle downwards
h_rect.x = h_gim.p0(1);
h_rect.y = left_ss_pf(2);
h_rect.w = h_gim.w;
h_rect.l =  ss_l_p2(2)-T_width - left_ss_pf(2);
h_rect.layer = 6;
h_rect.rounded = 0;
extended_shuttle_2 = rect(h_rect);

%Add etch holes to the shuttle extension 
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = [h_rect.x h_rect.y-2*(h_etch.r)];
section.type = 'rect';
section.w = h_rect.w;
section.l = h_rect.l+4*(h_etch.r);
h_etch.regions{1,1} = section;

SE_Eh_2 = etch_hole(h_etch);



%% Add the motors to the layout


if h_gim.manual == 0        %Should we draw a motor?
    %Right side motor
    load Motor_TT
    
    % h_motor.toothW = toothRad*1e6+overDraw;
    % h_motor.toothL = toothH*1e6;
    % h_motor.toothS = toothS*1e6-overDraw;
    % h_motor.shuttleG = gshut*1e6;
    % h_motor.fingerW = out.wffm*1e6+overDraw;
    % h_motor.fingerL = out.Lf*1e6;
    % h_motor.fingerSupportL = Ls*1e6;
    % h_motor.gap1 = out.g0*1e6-overDraw;
    % h_motor.gap2 = out.gb*1e6-overDraw;
    % h_motor.armW = warm*1e6+overDraw;
    % h_motor.armL =  out.Larm*1e6;
    % h_motor.armAngle = out.phifm;
    % h_motor.springW = springW*1e6+overDraw;
    % h_motor.springL = springL*1e6;
    % h_motor.gstopGap =  gstopGap*1e6-overDraw;
    
    
    h_motor.shuttle_w = 40;
    h_motor.pos = right_ga_points(1,:) + [h_motor.shuttle_w/2  0];         %set motor at bottom of gear rack
    h_motor.N = 25;
    h_motor.travel = 3000;
    h_motor.angle = 90;                                                     %Angle of shuttle in degrees
    h_motor.label = 'moto11';
    h_motor.num_inch_sets = 3;
    [m1 moto_pts]= motor(h_motor);
    
    %Add serpentine spring at the end of the motor
    make_bar_offset = 4;                        %make Bar funciton can change length of bar...
    h_ss.p1 = h_motor.pos - [0 moto_pts(1,1)-make_bar_offset];             % First point that the SS will span from
    h_ss.p2 = h_ss.p1 - [0 250];          % Second point that the SS will span to
    ss_l_p2 = h_ss.p2;
    h_ss.n = 8;                          % Number of meanders
    h_ss.w = 4;                           % Width of beams
    h_ss.dpp = 500;                       % Peak to peak distance of meanders
    h_ss.layer = 6;                       % Layer
    serpentine_moto_s = s_spring(h_ss);
    
    h_ss.rp = 1;
    points = s_spring(h_ss);        %matrix that has each vertex of SS
    h_ss.rp = 0;
    
    
    %Add anchor
    h_rect.w = 100;
    h_rect.l = 100;
    h_rect.x = h_ss.p2(1) - h_rect.w/2;
    h_rect.y = h_ss.p2(2) - h_rect.l;
    h_rect.rounded = 0;
    moto_ss_anc = rect(h_rect);
    
    
    
    % add left motor
    h_motor.pos = left_ga_points(1,:) - [h_motor.shuttle_w/2  0];         %set motor at bottom of gear rack
    h_motor.label = 'M2';
    [m2 moto_pts]= motor(h_motor);
    
    %Add serpentine spring at the end of the motor
    h_ss.p1 = h_motor.pos - [0 moto_pts(1,1)-make_bar_offset];             % First point that the SS will span from
    h_ss.p2 = h_ss.p1 - [0 250];          % Second point that the SS will span to
    ss_l_p2 = h_ss.p2;
    serpentine_moto_s_l = s_spring(h_ss);
    
    h_ss.rp = 1;
    points = s_spring(h_ss);        %matrix that has each vertex of SS
    h_ss.rp = 0;
    
    
    %Add anchor
    h_rect.w = 100;
    h_rect.l = 100;
    h_rect.x = h_ss.p2(1) - h_rect.w/2;
    h_rect.y = h_ss.p2(2) - h_rect.l;
    h_rect.rounded = 0;
    moto_ss_anc_l = rect(h_rect);
end


%% Add trench mask for chiplet
if h_gim.chiplet == 1
    %Add trench along entire robot
    chip_w = 4500;
    chip_l = 8000;
    trench_w = 400;
    y_travel = 1000;
    p0 = foot_tl - [trench_w/2  -trench_w/2];
    p1 = p0 + [foot_w+trench_w/2 0];
    p11 = p1 - [0 y_travel+trench_w/2];
    p12 = p11 + [(chip_w-foot_w)/2 0];
    p2 = p12 + [0 -chip_l+y_travel+trench_w/2];
    p3 = p2 + [-chip_w 0];
    p31 = p3 - [0 -chip_l+y_travel+trench_w/2];
    p32 = p31 + [(chip_w-foot_w)/2 0];
    p4 = p32 + [0 y_travel+trench_w/2];
    
    points = [p0;p1;p11;p12;p2;p3;p31;p32;p4];
    
    
    temp = gds_element('path', 'xy',points,'width',trench_w,'layer',5);
    str_name = sprintf('Trench_l_[%d,%d]',round(foot_tl(1)),round(foot_tl(2)));
    trench_cutout = gds_structure(str_name,temp);
    
    %add trench rectangle
    h_rect.w = foot_w;
    h_rect.l = y_travel + trench_w/2;
    h_rect.x = p32(1);
    h_rect.y = p32(2)-trench_w/2;
    h_rect.layer = 5;
    trench_rect_hop1 = rect(h_rect);
end

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