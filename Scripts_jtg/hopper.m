function out = hopper(h_gim)
% This function creates a fully synthesized MEMS Hopper 

if ~isfield(h_gim,'pin_gap')
    h_gim.pin_gap = 3;
end

if ~isfield(h_gim,'beak_l')
    h_gim.beak_l = 22;
end


if ~isfield(h_gim,'label')
    h_gim.label = 'Hopper';
end

if ~isfield(h_gim,'closed')
    h_gim.closed = 0;
end

if ~isfield(h_gim,'ss_contact_w')
    h_gim.ss_contact_w = 7;
end

if ~isfield(h_gim,'rack_teeth')
    h_gim.rack_teeth = 12;
end


if ~isfield(h_gim,'w')
    h_gim.w = 200;
end

if ~isfield(h_gim,'foot_etch_holes')
    h_gim.foot_etch_holes = 0;
end

if ~isfield(h_gim,'dx')
    h_gim.dx = 40;
end

if ~isfield(h_gim,'r')
    h_gim.r = [50 100];
end

if ~isfield(h_gim,'chiplet')
    h_gim.chiplet = 1;
end

if ~isfield(h_gim,'lever_arm_l')
    h_gim.lever_arm_l = 1200;
end


if ~isfield(h_gim,'manual')
    h_gim.manual = 0;
end

if ~isfield(h_gim,'noetch')
    h_gim.noetch = 0;
end

if ~isfield(h_gim,'bond_pad_gripper')
    h_gim.bond_pad_gripper = 0;
end

if ~isfield(h_gim,'shuttle_dx')
    h_gim.shuttle_dx = 1400;    %Total possible deflection on main shuttle
end


if ~isfield(h_gim,'trench_fill')
    h_gim.trench_fill = 0;
end

if ~isfield(h_gim,'motor_shuttle_w')
    h_gim.motor_shuttle_w = 40;
end

if ~isfield(h_gim,'internal_guide')
    h_gim.internal_guide = 0;
end



%% Create the Shuttle for GIM

% Number of iterations here will be a function of shuttle_dx
% Should be ~84um per h_gim.n, there is about 700um of dead space
% until the edge of the chip
h_gim.n = 8 + round(h_gim.shuttle_dx/84); 

h_gim.n = round(h_gim.shuttle_dx/80);

%h_gim.hw = 25;                  %Rough width of the hold (on the shuttle)
h_etch.noetch = h_gim.noetch;

h_gim.l = 2*(h_gim.n+1)*h_gim.dx;

p1 = h_gim.p0 + [h_gim.w 0];
p2 = p1 + [0 h_gim.l];
p3 = p2 - [h_gim.w 0];

%increase_dy_factor = h_gim.increase_factor; 
increase_dy_factor = 1.7; 
h_gim.hw = 30;      % This used to be 20, in original Lighthouse Design

%Creates the cutouts on the right side
temp_points_r = [];
central_push = 3;

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
    section2.p0 = [h_gim.p0(1)+h_gim.w/2 - central_push p4(2)-close_to_edge];
    section2.type = 'rcurve';
    section2.w = h_gim.w/2 + central_push;
    section2.l = increase_dy_factor*h_gim.dx+2*close_to_edge;
    section2.l = 2*h_gim.dx + 2*close_to_edge;
    
    section2.rcurve = final_pts';
    h_etch.regions = cell(1,1);
    h_etch.regions = {section2};
    h_etch.undercut = 5;
    h_etch.circle_etch = 1;
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
    section2.p0 = [h_gim.p0(1)+h_gim.w/2+central_push p4(2)-close_to_edge];
    section2.type = 'lcurve';
    section2.w = h_gim.w/2 + central_push;
    section2.l = 2*h_gim.dx+2*close_to_edge;
    section2.rcurve = final_pts';
    h_etch.regions = cell(1,1);
    h_etch.regions = {section2};
    h_etch.circle_etch = 1;
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

%Add top segment of main shuttle (no cutouts)
ms_w = h_gim.w;
ms_l = 750;

shift = 6;
h_rect.x = h_gim.p0(1) - ms_w/2 + h_gim.w/2;                  
h_rect.y = foot_y - shift;                  
h_rect.w = ms_w;                  
h_rect.l = ms_l + shift;
h_rect.etch = 1;           
main_shuttle_upper = rect(h_rect);    % Function to create a second rectangle GDS structure
h_rect.etch = 0;

%Add foot
foot_y = foot_y + ms_l;

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
h_etch.circle_etch = 1;
if h_gim.foot_etch_holes
    foot_eh = etch_hole(h_etch);
end

%Add dummy for the foot
h_rect.x = h_gim.p0(1) - foot_w/2 + h_gim.w/2;
h_rect.y = foot_y - h_gim.shuttle_dx;
h_rect.w = foot_w;

h_rect.l =  foot_h + h_gim.shuttle_dx;

h_rect.layer = 8;
h_rect.rounded = 0;
foot_s_df = rect(h_rect);

%Allow enough room in dummy for testing if not a chiplet
if h_gim.chiplet == 0
    dummy_gap = 50;
    
    h_rect.x = h_gim.p0(1) - foot_w/2 + h_gim.w/2-dummy_gap;
    h_rect.y = foot_y - h_gim.shuttle_dx;
    h_rect.w = foot_w + 2*dummy_gap;
    
    h_rect.l =  foot_h + h_gim.shuttle_dx + 400;
    
    h_rect.layer = 8;
    h_rect.rounded = 0;
    foot_s_df_2 = rect(h_rect);
end


if h_gim.manual == 0
    %Add backside foot re-trench
    edge_gap = 0;
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

% Add etch holes to the bottom rectangle of the shuttle with cutouts
closer_to_edge = 1;                     %Pushing etch holes closer to edges.
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = h_gim.p0 - [0 5*closer_to_edge];
section.type = 'rect';
section.w = h_gim.w;
section.l = h_gim.dx+8*closer_to_edge;
h_etch.regions{1,1} = section;
g_shuttle_bottom = etch_hole(h_etch);

% Add etch holes to the bottom rectangle of the shuttle with cutouts
% half width rectangle on bottom left of cutout shuttle
closer_to_edge = 1;                     
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = h_gim.p0 - [0 10*closer_to_edge] + [0 h_gim.dx+2*closer_to_edge];
section.type = 'rect';
section.w = h_gim.w/2;
section.l = h_gim.dx + 9*closer_to_edge;
h_etch.regions{1,1} = section;
g_shuttle_bottom2 = etch_hole(h_etch);


% Add etch holes to the top rectangle of the shuttle
h_etch.regions = cell(1,1);
h_etch.r = 2;
section.p0 = h_gim.p0 + [h_gim.w/2-6*closer_to_edge 2*h_gim.dx*(h_gim.n+.5)-6*closer_to_edge];
section.type = 'rect';
section.w = h_gim.w/2+8*closer_to_edge;
section.l = h_gim.dx+6*closer_to_edge+4;
h_etch.regions{1,1} = section;
g_shuttle_top = etch_hole(h_etch);

% Add guide on top for main shuttle (right and left shuttle guide)
gap = 2;
guide_w = 50;
guide_off = 120;
guide_l = ms_l + h_gim.l - h_gim.shuttle_dx-guide_off;
h_rect.x = h_gim.p0(1) - gap - guide_w + h_gim.w/2 - ms_w/2;
h_rect.y = h_gim.p0(2) + guide_off;
h_rect.w = guide_w;
h_rect.l = guide_l;
h_rect.xnum = 2;
h_rect.xspace = ms_w + 2*gap;
h_rect.layer = 6;
h_rect.rounded = 20;
left_GIM_guide = rect(h_rect);

% Add guide on bottom (left and right) of shuttle (under grab point)
gap = 2;
guide_ws = 80;
h_rect.x = h_gim.p0(1) - gap - guide_ws;
h_rect.y = h_gim.p0(2) - 500;
h_rect.w = guide_ws;
h_rect.l = 380;
h_rect.layer = 6;
h_rect.xspace = h_gim.w + 2*gap;
h_rect.rounded = 20;
left_GIM_guide_bot = rect(h_rect);
h_rect.xnum = 1;

%Add bottom dummy fill near the guides
h_rect.x = h_gim.p0(1) - (ms_w - h_gim.w)/2 - 25;
h_rect.y = h_gim.p0(2);
h_rect.w = ms_w + 50;
h_rect.l = h_gim.shuttle_dx + ms_l;
h_rect.layer = 8;
h_rect.rounded = 0;
df_guide_area = rect(h_rect);

h_rect.x = h_gim.p0(1) - 4*gap;
h_rect.y = h_gim.p0(2) - 750;
h_rect.w = h_gim.w + 8*gap;
h_rect.l = 800;
h_rect.layer = 8;
h_rect.rounded = 0;
df_guide_area_bottom = rect(h_rect);


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

h_latch.ss.n = 7;                             % Number meanders in serpentine spring on rotor
h_latch.ss.dpp = 70*2;                          %was 70 in the designs that worked
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
h_latch.noetch = h_gim.noetch;
h_latch.stages = 1;
h_latch.mech_latch = 0;
h_latch.compact_latch = 0;
h_latch.no_backstops = 0;
h_latch.pin_gap = h_gim.pin_gap;
[g_latch latch_ref_pts]= latch(h_latch);

%Add right geared lever arm
h_arm.n = 100;
h_arm.opening_theta = 360 - h_latch.actuation_angle;
h_arm.inner_arm_theta = 180 - h_latch.actuation_angle/2;
h_arm.p0 = h_latch.p0;
h_arm.r = [h_gim.r(2)-6 h_gim.lever_arm_l];      %6 is the etch hole undercut
h_arm.orientation = h_latch.orientation;
h_arm.layer = 6;
h_arm.shuttle_w = 50;
h_arm.rack_teeth = h_gim.rack_teeth;
h_arm.manual = h_gim.manual;
h_arm.noetch = h_gim.noetch;
h_arm.rack_offset = 172;

[right_ga right_ga_points] = gear_arm(h_arm);

%%
%Add right side serpentine springs
% Generate Top right Serpentine Spring
num_meanders = h_gim.ss_n;
spring_l = h_gim.ss_l;
spring_w = h_gim.ss_w;
meander_gap = h_gim.ss_gap;
spring_buffer = 20;         %Gap between springs and stationary things

h_ss.p1 = h_arm.p0 - [(h_gim.r(2) - spring_l/2 - spring_buffer) 500];                       % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 2*num_meanders*(spring_w+meander_gap)];   % Second point that the SS will span to
ss_r_p2 = h_ss.p2;
h_ss.n = num_meanders;                                 % Number of meanders 
h_ss.w = spring_w;                                     % Width of beams 
h_ss.dpp = spring_l;                                   % Peak to peak distance of meanders
h_ss.layer = 6;                                        % Layer 
serpentine_h = s_spring(h_ss);

h_ss.rp = 1;
points = s_spring(h_ss);        %matrix that has each vertex of SS
h_ss.rp = 0;



%Spring protect line right
protection_w = 50;
%protection_l = 3150;
protection_l = 2*(abs(h_ss.p2(2) - h_ss.p1(2)) + h_gim.shuttle_dx) + 300;   %300 from anchors etc


h_rect.w = protection_w;
h_rect.l = protection_l;
h_rect.x = h_ss.p1(1)+h_ss.dpp/2+15;
protectoin_point = h_rect.x;
h_rect.y = h_ss.p1(2)-h_rect.l+50;

%For motor routing later. The left side routing boundary (x coordinate) for
%the right side motor
left_routing_boundary_r = h_rect.x + h_rect.w;


h_rect.layer = 6;
h_rect.rounded = 0;
h_rect.xnum = 1;
h_rect.xspace = 10;
prot1 = rect(h_rect);

%Add anchor
h_rect.w = spring_l;
h_rect.l = 100;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
r_ss_anc = rect(h_rect);

%Second copy of right hand spring (on bottom)
temp = h_ss.p1(2) - h_ss.p2(2);
h_ss.p1 = h_ss.p1 - [0 temp+50+100+h_gim.shuttle_dx];            
h_ss.p2 = h_ss.p2 - [0 temp+50+100+h_gim.shuttle_dx];          
serpentine_h_rs_2 = s_spring(h_ss);


right_ss_pf = h_ss.p2;

%Add anchor
h_rect.w = spring_l;
h_rect.l = 100;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
r_ss_anc2 = rect(h_rect);

%Reset ffor future use
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

[g_latch2 latch_ref_pts]= latch(h_latch);

%Add left geared lever arm
h_arm.n = 100;
h_arm.opening_theta = 360 - h_latch.actuation_angle;
h_arm.inner_arm_theta = 0 + h_latch.actuation_angle/2;
h_arm.p0 = h_latch.p0;
h_arm.r = [h_gim.r(2)-6 h_gim.lever_arm_l];      %6 is the etch hole undercut
h_arm.layer = 6;
h_arm.orientation = h_latch.orientation;
h_arm.shuttle_w = 50;
h_arm.rack_offset = 166;

[left_ga left_ga_points]= gear_arm(h_arm);

%Add left side serpentine spring
%Generate Serpentine Spring
h_ss.p1 = h_arm.p0 - [-(h_gim.r(2) - spring_l/2 - spring_buffer) 500];         % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 2*num_meanders*(spring_w+meander_gap)];          % Second point that the SS will span to
ss_l_p2 = h_ss.p2;
serpentine_h_ls = s_spring(h_ss);

%left side protection spring
h_rect.w = protection_w;
h_rect.l = protection_l;
h_rect.x = h_ss.p1(1)-h_ss.dpp/2-35 - protection_w + 20;
h_rect.y = h_ss.p1(2)-h_rect.l+50;
h_rect.layer = 6;
h_rect.xnum = 1;
h_rect.xspace = 10;
h_rect.rounded = 0;
prot2 = rect(h_rect);
h_rect.xnum  = 1;

right_routing_boundary_l = h_rect.x;


%Add dummy fill for the spring system
h_rect.w = protectoin_point - h_rect.x + protection_w;
h_rect.l = protection_l;
h_rect.x = h_ss.p1(1)-h_ss.dpp/2-35 - protection_w + 20;
h_rect.y = h_ss.p1(2)-h_rect.l+50;
routing_SS_barrier = h_rect.y-200;
h_rect.layer = 8;
h_rect.rounded = 0;
prot2_df = rect(h_rect);


%Add anchor
h_rect.w = spring_l;
h_rect.l = 100;
h_rect.layer = 6;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
l_ss_anc2 = rect(h_rect);


%Second copy
temp = h_ss.p1(2) - h_ss.p2(2);
h_ss.p1 = h_ss.p1 - [0 temp+50+100+h_gim.shuttle_dx];         
h_ss.p2 = h_ss.p2 - [0 temp+50+100+h_gim.shuttle_dx];          
serpentine_h_ls_2 = s_spring(h_ss);

left_ss_pf = h_ss.p2;

h_ss.rp = 1;
points = s_spring(h_ss);        %matrix that has each vertex of SS
h_ss.rp = 0;

%Add anchor
h_rect.w = spring_l;
h_rect.l = 100;
h_rect.layer = 6;
h_rect.x = h_ss.p1(1) - h_rect.w/2;
h_rect.y = h_ss.p1(2);
h_rect.rounded = 0;
l_ss_anc1 = rect(h_rect);

%Rest vars for future use
h_ss.p1 = h_arm.p0 - [250 500];         % First point that the SS will span from
h_ss.p2 = h_ss.p1 - [0 2*num_meanders*(spring_w+meander_gap)];          % Second point that the SS will span to



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
h_etch.circle_etch = 1;

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
h_etch.circle_etch = 1;
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
h_etch.circle_etch = 1;
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

% Add an internal guide to the bottom most part of the shuttle
if h_gim.internal_guide
    internal_guide_w = 100;
    y_gap = 58;             % offset from bottom of main shuttle
    slot_guide_gap = 50;
    
    tolerance = 50;
    
    guide_l = 200;
    
    h_rect.x = h_gim.p0(1) + h_gim.w/2 - internal_guide_w/2;
    h_rect.y = left_ss_pf(2)+y_gap;
    h_rect.w = internal_guide_w;
    h_rect.l =  h_gim.shuttle_dx + guide_l + slot_guide_gap + tolerance + 1;
    h_rect.layer = 2;
    h_rect.rounded = 20;
    internal_guide_slot = rect(h_rect);
    
    h_rect.x = h_gim.p0(1) + h_gim.w/2 - internal_guide_w/2 + gap;
    h_rect.y = left_ss_pf(2)+y_gap + slot_guide_gap;
    h_rect.w = internal_guide_w - 2*gap;
    h_rect.l = guide_l;
    h_rect.layer = 9;
    h_rect.rounded = 20;
    internal_guide = rect(h_rect);
    h_rect.rounded = 0;
end
%% Add trench mask for chiplet

%Add trench along entire robot
chip_w = h_gim.chip_w; % 5500
chip_l = h_gim.chip_l; % 8000
trench_w = 400;
tether_w = 200;
y_travel = h_gim.shuttle_dx;
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

routing_edge_gap = 100;     % Space to leave around edge of chip to protect routing
right_routing_boundary_r = p2(1) - trench_w/2 - routing_edge_gap;
left_routing_boundary_l = p3(1) + trench_w/2 + routing_edge_gap;

if h_gim.chiplet == 1
    temp = gds_element('path', 'xy',points,'width',trench_w,'layer',5);
    str_name = sprintf('Trench_l_[%d,%d]',round(foot_tl(1)),round(foot_tl(2)));
    trench_cutout = gds_structure(str_name,temp);
    
    temp = gds_element('path', 'xy',points,'width',trench_w,'layer',8);
    str_name = sprintf('Trench_l_df_[%d,%d]',round(foot_tl(1)),round(foot_tl(2)));
    trench_cutout_df = gds_structure(str_name,temp);
    
    %add trench rectangle
    h_rect.w = foot_w;
    h_rect.l = y_travel + trench_w/2;
    h_rect.x = p32(1);
    h_rect.y = p32(2)-trench_w/2;
    h_rect.layer = 5;
    trench_rect_hop1 = rect(h_rect);
    
    % Add tether
    h_rect.w = tether_w;
    h_rect.l = 400;
    h_rect.x = mean([p2(1) p3(1)]) - 100;
    h_rect.y = p2(2)-trench_w/2;
    h_rect.layer = 4;
    tether_robo = rect(h_rect);
    
    % Add Label
    h_label.text = h_gim.label;
    h_label.layer = 2;
    h_label.height = 50;
    h_label.p0 = p31 + [400 -400];
    robo_lab = add_label(h_label);  % Add label underneath rectangle
    
    h_rect.etch = 0;
end


% This code executes if the trench should be filled with dummy silicon to
% attempt to prevent breakthrough during the backside etch
if h_gim.trench_fill
    s_gap = 20;                             % Gap between dummy fill and structures
    dev_tether_center_offset = 500;         % Offset from center of chip to tethers in center
    
    t0 = p11 + [trench_w/2-s_gap 0];
    t1 = p12;
    t2 = p2;
    t3 = t2 - [chip_w/2-s_gap-tether_w/2 0];
    
    points = [t0;t1;t2;t3];
    temp = gds_element('path', 'xy',points,'width',trench_w-2*s_gap,'layer',9);
    str_name = sprintf('SOI_protection_[%d,%d]',round(t0(1)),round(t0(2)));
    soi_pro_1 = gds_structure(str_name,temp);
    
    t4 = t3 - [tether_w+2*s_gap 0];
    t5 = p3;
    t6 = p31;
    t7 = p32 - [trench_w/2-s_gap 0];
    
    points = [t4;t5;t6;t7];
    temp = gds_element('path', 'xy',points,'width',trench_w-2*s_gap,'layer',9);
    str_name = sprintf('SOI_protection_[%d,%d]',round(t4(1)),round(t4(2)));
    soi_pro_2 = gds_structure(str_name,temp);
    
    % Add releasable buckle over the tether
    h_rect.w = tether_w + 2*s_gap;
    h_rect.l = trench_w - 2*s_gap;
    h_rect.x = t4(1);
    h_rect.y = t4(2) - (trench_w/2-s_gap);
    h_rect.layer = 6;
    h_rect.etch = 1;
    soi_pro_buckle = rect(h_rect);
    h_rect.etch = 0;
    
    % Add small device side tether at top right of chip
    square_w = 50;
    path_w = 10;
    path_g = 2;
    
    p0 = t1;
    p1 = t1 - .7*[trench_w trench_w];
    h_path.pts = [p0;p1];
    h_path.w = path_w;
    h_path.layer = 6;
    dev_tether_tr = m_path(h_path);
    
    h_path.w = path_w+2*path_g;
    h_path.layer = 8;
    dev_tether_tr_df = m_path(h_path);
    
    
    h_rect.w = square_w;
    h_rect.l = square_w;
    h_rect.x = p1(1) - square_w/2;
    h_rect.y = p1(2) - square_w/2;
    h_rect.layer = 6;
    h_rect.etch = 1;
    h_rect.UNDERCUT = 4;
    soi_pro_square_1 = rect(h_rect);
    h_rect.etch = 0;
    h_rect.UNDERCUT = 6;
    
    h_rect.w = square_w+2*path_g;
    h_rect.l = square_w+2*path_g;
    h_rect.x = p1(1) - square_w/2-path_g;
    h_rect.y = p1(2) - square_w/2-path_g;
    h_rect.layer = 8;
    h_rect.etch = 0;
    soi_pro_square_1_df = rect(h_rect);
    
    % Add small device side tether at top center right of chip    
    p0 = t1 - [chip_w/2-dev_tether_center_offset 0];
    p1 = p0 - .7*[0 trench_w];
    h_path.pts = [p0;p1];
    h_path.w = path_w;
    h_path.layer = 6;
    dev_tether_tr_c = m_path(h_path);
    
    h_path.w = path_w+2*path_g;
    h_path.layer = 8;
    dev_tether_tr_df_c = m_path(h_path);
    
    
    h_rect.w = square_w;
    h_rect.l = square_w;
    h_rect.x = p1(1) - square_w/2;
    h_rect.y = p1(2) - square_w/2;
    h_rect.layer = 6;
    h_rect.etch = 1;
    h_rect.UNDERCUT = 4;
    soi_pro_square_1_c = rect(h_rect);
    h_rect.etch = 0;
    h_rect.UNDERCUT = 6;
    
    h_rect.w = square_w+2*path_g;
    h_rect.l = square_w+2*path_g;
    h_rect.x = p1(1) - square_w/2-path_g;
    h_rect.y = p1(2) - square_w/2-path_g;
    h_rect.layer = 8;
    h_rect.etch = 0;
    soi_pro_square_1_df_c = rect(h_rect);
    
    % Add small device side tether at top left of chip
    p0 = t6;
    p1 = t6 + .7*[trench_w -trench_w];
    h_path.pts = [p0;p1];
    h_path.w = 10;
    h_path.layer = 6;
    dev_tether_tl = m_path(h_path);
    
    h_path.w = path_w+2*path_g;
    h_path.layer = 8;
    dev_tether_tl_df = m_path(h_path);
    
    
    h_rect.w = square_w;
    h_rect.l = square_w;
    h_rect.x = p1(1) - square_w/2;
    h_rect.y = p1(2) - square_w/2;
    h_rect.layer = 6;
    h_rect.etch = 1;
    h_rect.UNDERCUT = 4;
    soi_pro_square_1l = rect(h_rect);
    h_rect.etch = 0;
    h_rect.UNDERCUT = 6;
    
    h_rect.w = square_w+2*path_g;
    h_rect.l = square_w+2*path_g;
    h_rect.x = p1(1) - square_w/2-path_g;
    h_rect.y = p1(2) - square_w/2-path_g;
    h_rect.layer = 8;
    h_rect.etch = 0;
    soi_pro_square_1_df2 = rect(h_rect);
    
    
    % Add small device side tether at top center left of chip
    p0 = t6 + [chip_w/2-dev_tether_center_offset 0];
    p1 = p0 + .7*[0 -trench_w];
    h_path.pts = [p0;p1];
    h_path.w = 10;
    h_path.layer = 6;
    dev_tether_tl_ce = m_path(h_path);
    
    h_path.w = path_w+2*path_g;
    h_path.layer = 8;
    dev_tether_tl_df_ce = m_path(h_path);
    
    h_rect.w = square_w;
    h_rect.l = square_w;
    h_rect.x = p1(1) - square_w/2;
    h_rect.y = p1(2) - square_w/2;
    h_rect.layer = 6;
    h_rect.etch = 1;
    h_rect.UNDERCUT = 4;
    soi_pro_square_1l_ce = rect(h_rect);
    h_rect.etch = 0;
    h_rect.UNDERCUT = 6;
    
    h_rect.w = square_w+2*path_g;
    h_rect.l = square_w+2*path_g;
    h_rect.x = p1(1) - square_w/2-path_g;
    h_rect.y = p1(2) - square_w/2-path_g;
    h_rect.layer = 8;
    h_rect.etch = 0;
    soi_pro_square_1_df2_ce = rect(h_rect);
    
    
    
    
    
    % Long skinny protection bars at left and right of foot
    h_rect.w = trench_w/2 - 2*s_gap;
    h_rect.l = h_gim.shuttle_dx + trench_w + foot_h - 2*s_gap;
    h_rect.x = t7(1);
    h_rect.y = t7(2) - (trench_w/2-s_gap);
    h_rect.layer = 9;
    h_rect.xnum = 2;
    h_rect.xspace = foot_w + 2*s_gap;
    soi_pro3 = rect(h_rect);
    h_rect.xnum = 1;
    
    % Top protection bar
    h_rect.w = foot_w + 2*s_gap;
    h_rect.l = trench_w - 2*s_gap;
    h_rect.x = t7(1) + trench_w/2 - 2*s_gap;
    h_rect.y = t7(2) - (trench_w/2-s_gap) + (h_gim.shuttle_dx + trench_w + foot_h - 2*s_gap) - trench_w + 2*s_gap;
    h_rect.layer = 9;
    soi_pro4 = rect(h_rect);
    
    % Protection under the foot (2 rects)
    h_rect.w = (foot_w - h_gim.w)/2;
    h_rect.l = h_gim.shuttle_dx - 2*s_gap;
    h_rect.x = t7(1) + trench_w/2 - 2*s_gap;
    h_rect.y = t7(2) - (trench_w/2-s_gap);
    h_rect.layer = 9;
    h_rect.xnum = 2;
    h_rect.xspace = h_gim.w + 2*s_gap;
    soi_pro5 = rect(h_rect);
    h_rect.xnum = 1;
    
    
end


%% Add the motors to the layout
if h_gim.manual == 0        %Should we draw a motor?
    %Right side motor
    load Motor_TT    
    h_motor.layer = [6 2];
    h_motor.shuttle_w = h_gim.motor_shuttle_w;
    h_motor.pos = right_ga_points(1,:) + [h_motor.shuttle_w/2  0];         %set motor at bottom of gear rack
    h_motor.N = h_gim.motorN;
    h_motor.travel = 2000;
    h_motor.angle = 90;                                                     %Angle of shuttle in degrees
    h_motor.hopper_shuttle_extend = 1;                                      % Extend shuttle for etch hole alignment
    h_motor.label = 'moto11';
    h_motor.circle_etch_holes = 1;
    h_motor.num_inch_sets = 2;
    h_motor.ground_serpentine = 1;
    [m1 moto_pts]= motor(h_motor);
    
    % Add manual-ish dummy fill for motors
    %cent_to_cent = 1830;            % Distance from center to center of backstops
    cent_to_cent = 1020 + (h_gim.motorN-20)*22.5;            % Distance from center to center of backstops
    array_w = 330;                  % Width of the GCA array
    array_w_decrease = 55;          % Amount to pull width of total array (for routing)
                                    % Effectively pulls routing contact
                                    % closer to center of stator anchor
       
    for ii = 2:2:length(moto_pts)
        %Find the bottom left point
        px = mean([moto_pts(ii,1) moto_pts(ii+1,1)])-cent_to_cent/2;
        py = moto_pts(ii,2) - array_w/2;
        h_rect.x = px;                   
        h_rect.y = py;                   
        h_rect.w = cent_to_cent;                   
        h_rect.l = array_w;                 
        h_rect.layer = 8;             
        str = sprintf('df_%d = rect(h_rect);',ii);         
        eval(str);
    end
    
    % Routing - for the motors on the right hand side
    rps = zeros(4,2); 
    lrps = zeros(4,2);
    
    for jj = 1:2                % 1 is the top motors, 2 is the bottom motors
        for kk = 1:-2:-1        % positive KK is on left side of shuttle, negative is on right side of shuttle
        route_pts = zeros(2*h_motor.num_inch_sets,2);
        for ii = 2:2:length(route_pts)
            route_pts(ii-1,1) = mean([moto_pts(ii+(jj-1)*4,1) moto_pts(ii+(jj-1)*4+1,1)])-cent_to_cent*kk/4;
            route_pts(ii-1,2) = moto_pts(ii+(jj-1)*4,2) + (array_w-array_w_decrease)/2;
            route_pts(ii,1) = mean([moto_pts(ii+(jj-1)*4,1) moto_pts(ii+(jj-1)*4+1,1)])-cent_to_cent*kk/4;
            route_pts(ii,2) = moto_pts(ii+(jj-1)*4,2) - (array_w-array_w_decrease)/2;
        end
        
        
        h_routing.w = 50;
        
        % Routing - Hard coded for 2 sets of motors
        pitch = route_pts(1,2) - route_pts(3,2);
        
        p1 = route_pts(1,:);
        if kk>0
            p2 = [left_routing_boundary_r + kk*((jj-1)*h_routing.w  + h_routing.w/2) route_pts(1,2)];
        else
            p2 = [right_routing_boundary_r + kk*((jj-1)*h_routing.w  + h_routing.w/2) route_pts(1,2)];
        end
        
        p3 = [p2(1) route_pts(3,2)];
        p4 = route_pts(3,:);
        p5 = route_pts(2,:);
        p6 = [p2(1)  route_pts(3,2)-pitch+h_routing.w-2];           %2 is from phase offset between pawls
        if kk>0
            if jj==1
                rps(1,:) = p6;
            else
                rps(2,:) = p6;
            end
        else
            if jj==1
                rps(4,:) = p6;
            else
                rps(3,:) = p6;
            end
        end
        
        p7 = [route_pts(4,1) p6(2)];
        p8 = route_pts(4,:);
        
        pts = [p1;p2;p3;p4;p5];
        pts2 = [p3;p6;p7;p8];
        
        h_routing.pts = pts;
        
        h_routing.layer = 6;
        str = sprintf('route1_%d_%d = m_routing(h_routing);',jj,kk+2);
        eval(str);
        
        h_routing.pts = pts2;
        str = sprintf('route2_%d_%d = m_routing(h_routing);',jj,kk+2);
        eval(str);
        end
    end
    
    % Routing - Motor have been wired initially, need to connect sides
    % Bring oustide routing lines down to level with the inside ones
    % Right side motor
    pts = [rps(1,:); [rps(1,1) rps(2,2)]];
    h_routing.pts = pts;
    rps_1_left = m_routing(h_routing);
    rps(1,:) = [rps(1,1) rps(2,2)];
    
    pts = [rps(4,:); [rps(4,1) rps(3,2)]];
    h_routing.pts = pts;
    rps_1_right = m_routing(h_routing);
    rps(4,:) = [rps(4,1) rps(3,2)];
    
    % Left side motor
    
    
    
    
    
    
    
    %Add serpentine spring at the end of the motors
    SS_total_deflection = 850;                % Total deflection of these SS (for dummy fill)
    
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
        
    
    %Add Dummy exclude over serpentine spring
    dummy_gap = 30;
    h_rect.w = h_ss.dpp+dummy_gap;
    h_rect.l = abs(h_ss.p2(2) - h_ss.p1(2)) + SS_total_deflection;
    h_rect.layer = 8;  
    h_rect.x = h_ss.p1(1) - (h_ss.dpp + dummy_gap)/2;
    h_rect.y = h_ss.p2(2);
    h_rect.rounded = 0;
    df_SS_right = rect(h_rect);
    
    % Round the corners where the shuttle enters the guides
    h_rect.w = (h_ss.dpp+dummy_gap-h_motor.shuttle_w - 2*4)/2;
    h_rect.l = 20;
    h_rect.layer = 6;  
    h_rect.x = h_ss.p1(1) - (h_ss.dpp + dummy_gap)/2;
    h_rect.y = h_ss.p2(2)+abs(h_ss.p2(2) - h_ss.p1(2)) + SS_total_deflection-10;
    h_rect.rounded = 10;
    h_rect.xnum = 2;
    h_rect.xspace = h_motor.shuttle_w + 2*4;
    df_corner_round = rect(h_rect);
    h_rect.xnum = 1;
    
    %Add anchor
    h_rect.w = h_ss.dpp+300;
    h_rect.l = 100;
    h_rect.layer = 6;  
    h_rect.x = h_ss.p2(1) - h_rect.w/2;
    h_rect.y = h_ss.p2(2) - h_rect.l;
    h_rect.rounded = 0;
    moto_ss_anc = rect(h_rect);
    
    
    
    % add left motor
    h_motor.pos = left_ga_points(1,:) - [h_motor.shuttle_w/2  0];         %set motor at bottom of gear rack
    h_motor.label = 'M2';
    [m2 moto_pts]= motor(h_motor);
     
    for ii = 2:2:length(moto_pts)
        %Find the bottom left point
        px = mean([moto_pts(ii,1) moto_pts(ii+1,1)])-cent_to_cent/2;
        py = moto_pts(ii,2) - array_w/2;
        h_rect.x = px;                   
        h_rect.y = py;                   
        h_rect.w = cent_to_cent;                   
        h_rect.l = array_w;                 
        h_rect.layer = 8;             
        str = sprintf('dff_%d = rect(h_rect);',ii);         
        eval(str);
    end
    
    % Routing - for the left motors on the left hand side
    for jj = 1:2
        for kk = 1:-2:-1
        route_pts = zeros(2*h_motor.num_inch_sets,2);
        for ii = 2:2:length(route_pts)
            route_pts(ii-1,1) = mean([moto_pts(ii+(jj-1)*4,1) moto_pts(ii+(jj-1)*4+1,1)])-cent_to_cent*kk/4;
            route_pts(ii-1,2) = moto_pts(ii+(jj-1)*4,2) + (array_w-array_w_decrease)/2;
            route_pts(ii,1) = mean([moto_pts(ii+(jj-1)*4,1) moto_pts(ii+(jj-1)*4+1,1)])-cent_to_cent*kk/4;
            route_pts(ii,2) = moto_pts(ii+(jj-1)*4,2) - (array_w-array_w_decrease)/2;
        end
        
        
        h_routing.w = 50;
        
        % Routing - Hard coded for 2 sets of motors
        pitch = route_pts(1,2) - route_pts(3,2);
        
        p1 = route_pts(1,:);
        if kk>0
            p2 = [left_routing_boundary_l + kk*((jj-1)*h_routing.w  + h_routing.w/2) route_pts(1,2)];
        else
            p2 = [right_routing_boundary_l + kk*((jj-1)*h_routing.w  + h_routing.w/2) route_pts(1,2)];
        end
        
        p3 = [p2(1) route_pts(3,2)];
        p4 = route_pts(3,:);
        p5 = route_pts(2,:);
        p6 = [p2(1)  route_pts(3,2)-pitch+h_routing.w-2];           %2 is from phase offset between pawls
        
        if kk>0
            if jj==1
                lrps(1,:) = p6;
            else
                lrps(2,:) = p6;
            end
        else
            if jj==1
                lrps(4,:) = p6;
            else
                lrps(3,:) = p6;
            end
        end
        
        
        p7 = [route_pts(4,1) p6(2)];
        p8 = route_pts(4,:);
        
        pts = [p1;p2;p3;p4;p5];
        pts2 = [p3;p6;p7;p8];
        
        h_routing.pts = pts;
        
        h_routing.layer = 6;
        str = sprintf('route3_%d_%d = m_routing(h_routing);',jj,kk+2);
        eval(str);
        
        h_routing.pts = pts2;
        str = sprintf('route4_%d_%d = m_routing(h_routing);',jj,kk+2);
        eval(str);
        end
    end
    
    % Routing - Motor have been wired initially, need to connect sides
    % Bring oustide routing lines down to level with the inside ones
    % Left side motor
    pts = [lrps(1,:); [lrps(1,1) lrps(2,2)]];
    h_routing.pts = pts;
    lrps_1_left = m_routing(h_routing);
    lrps(1,:) = [lrps(1,1) lrps(2,2)];
    
    pts = [lrps(4,:); [lrps(4,1) lrps(3,2)]];
    h_routing.pts = pts;
    lrps_1_right = m_routing(h_routing);
    lrps(4,:) = [lrps(4,1) lrps(3,2)];
        
    %Add serpentine spring at the end of the motor
    h_ss.p1 = h_motor.pos - [0 moto_pts(1,1)-make_bar_offset];             % First point that the SS will span from
    h_ss.p2 = h_ss.p1 - [0 250];          % Second point that the SS will span to
    ss_l_p2 = h_ss.p2;
    serpentine_moto_s_l = s_spring(h_ss);
    
    h_ss.rp = 1;
    points = s_spring(h_ss);        %matrix that has each vertex of SS
    h_ss.rp = 0;
    
    %Add Dummy exclude over serpentine spring
    h_rect.w = h_ss.dpp+dummy_gap;
    h_rect.l = abs(h_ss.p2(2) - h_ss.p1(2)) + SS_total_deflection;
    h_rect.layer = 8;  
    h_rect.x = h_ss.p1(1) - (h_ss.dpp + dummy_gap)/2;
    h_rect.y = h_ss.p2(2);
    h_rect.rounded = 0;
    df_SS_left = rect(h_rect);
    
    hx = h_rect.x;
    hy = h_rect.y;
    hw = h_rect.w;
    
    % Round the corners where the shuttle enters the guides
    h_rect.w = (h_ss.dpp+dummy_gap-h_motor.shuttle_w - 2*4)/2;
    h_rect.l = 20;
    h_rect.layer = 6;  
    h_rect.x = h_ss.p1(1) - (h_ss.dpp + dummy_gap)/2;
    h_rect.y = h_ss.p2(2)+abs(h_ss.p2(2) - h_ss.p1(2)) + SS_total_deflection-10;
    h_rect.rounded = 10;
    h_rect.xnum = 2;
    h_rect.xspace = h_motor.shuttle_w + 2*4;
    df_corner_round_l = rect(h_rect);
    h_rect.xnum = 1;
    
    
    % Create pads to route signals to
    num_pads = 5;
    pad_w = 500;
    pad_l = 700;
    gap = 40;
    offset = 64;            % Offset in x direction
    
    %Bottom left corner where bond pads will start
    init = [hx+hw-230 hy - pad_l - 300];
    
    %Contains coordinates of center of all pads
    pad_pts = zeros(num_pads,2);
    
    for i=1:num_pads
        % Make outline of not dummy
        h_rect.x = offset + init(1) + (i-1)*(pad_w+2*gap) - 2*gap;
        h_rect.y = init(2) - 2*gap;
        h_rect.w = pad_w + 2*gap;
        h_rect.l = pad_l + 2*gap;
        h_rect.layer = 8;
        h_rect.rounded = 0;
        str = sprintf('df_pads_%d = rect(h_rect);',i);
        eval(str);
        pad_pts(i,1) = h_rect.x + h_rect.w/2;
        pad_pts(i,2) = h_rect.y + h_rect.l/2;
        top_bond_pads = h_rect.y + h_rect.l;
        
        ry = h_rect.y;
        
        if h_gim.bond_pad_gripper
            % Made a bond pad that has a gripper
            h_bpg.p0 = [h_rect.x+gap h_rect.y+gap];
            h_bpg.w = pad_w;
            h_bpg.l = pad_l;
            h_bpg.ring = 100;
            h_bpg.dx = 15;
            h_bpg.r = 25;
            
            str = sprintf('pad_gripper_%d = bond_pad_gripper(h_bpg);',i);
            eval(str);
            
            % Add trench behind the bond pads
            h_rect.x = h_bpg.p0(1)-gap;
            h_rect.y = h_bpg.p0(2) + (pad_l-400-100);
            h_rect.w = pad_w+2*gap;
            h_rect.l = 400;
            h_rect.layer = 5;
            h_rect.rounded = 0;
            str = sprintf('trench_pad_%d = rect(h_rect);',i);
            eval(str);
            
            
        else
            % Standard SOI only bond pad
            h_rect.x = h_rect.x + gap;
            h_rect.y = h_rect.y + gap;
            h_rect.w = pad_w;
            h_rect.l = pad_l;
            h_rect.layer = 6;
            h_rect.rounded = 0;
            str = sprintf('padss_%d = rect(h_rect);',i);
            eval(str);
        end
        

        

        

    end

    
    
    % Right side motor routing
    routing_lower_most_barrier = ry - 150;
    
    %Set the routing boundary to the top of the bond pads
    routing_SS_barrier = top_bond_pads + h_routing.w/2;
    
    % Routing top right HV lines to their pads for right side motor
    
    % Pad routing offset
    route_offset = [0 pad_l/2];
    
    % Top left HV line (farthest to the left)
    p1 = rps(1,:);
    p2 = [rps(1,1) routing_SS_barrier];
    p3 = [pad_pts(4,1) p2(2)];
    p4 = pad_pts(4,:) + route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    route_HV_right_top_left = m_routing(h_routing);
    
    % Bottom left HV line (2nd from left)
    p1 = rps(2,:);
    p2 = [rps(2,1) routing_SS_barrier+h_routing.w];
    p3 = [pad_pts(5,1) p2(2)];
    p4 = pad_pts(5,:) + route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    route_HV_right_top_left_2 = m_routing(h_routing);
    
    % Top Right HV line (fartherst to the right)
    p1 = rps(4,:);
    p2 = [rps(4,1) routing_lower_most_barrier];
    p3 = [pad_pts(4,1) p2(2)];
    p4 = pad_pts(4,:) - route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    route_HV_right_top_left_3 = m_routing(h_routing);
    
    % Bottom right HV line (2nd from the right)
    p1 = rps(3,:);
    p2 = [rps(3,1) routing_lower_most_barrier+h_routing.w];
    p3 = [pad_pts(5,1) p2(2)];
    p4 = pad_pts(5,:) - route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    route_HV_right_top_left_4 = m_routing(h_routing);
    
    
    
    % Left motor routing
    
    rps = lrps;
    % Routing bottom left HV lines to their pads for left side motor
    
    % Top left HV line (farthest to the left)
    p1 = rps(1,:);
    p2 = [rps(1,1) routing_lower_most_barrier];
    p3 = [pad_pts(2,1) p2(2)];
    p4 = pad_pts(2,:) - route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    lroute_HV_right_top_left = m_routing(h_routing);
    
    % Bottom left HV line (2nd from left)
    p1 = rps(2,:);
    p2 = [rps(2,1) routing_lower_most_barrier+h_routing.w];
    p3 = [pad_pts(1,1) p2(2)];
    p4 = pad_pts(1,:) - route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    lroute_HV_right_top_left_2 = m_routing(h_routing);
    
    % Top right HV line (fartherst to the right)
    p1 = rps(4,:);
    p2 = [rps(4,1) routing_SS_barrier];
    p3 = [pad_pts(2,1) p2(2)];
    p4 = pad_pts(2,:) + route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    lroute_HV_right_top_left_3 = m_routing(h_routing);
    
    % Bottom right HV line (2nd from the right)
    p1 = rps(3,:);
    p2 = [rps(3,1) routing_SS_barrier+h_routing.w];
    p3 = [pad_pts(1,1) p2(2)];
    p4 = pad_pts(1,:) + route_offset;
    h_routing.pts = [p1;p2;p3;p4];
    h_routing.layer = 6;
    lroute_HV_right_top_left_4 = m_routing(h_routing);
    
    %Add anchor
    h_rect.w = h_ss.dpp+300;
    h_rect.l = 100;
    h_rect.x = h_ss.p2(1) - h_rect.w/2;
    h_rect.y = h_ss.p2(2) - h_rect.l;
    h_rect.layer = 6;
    h_rect.rounded = 0;
    moto_ss_anc_l = rect(h_rect);
    
    % Add connection for middle (ground) pad on top
    h_rect.w = h_routing.w;
    h_rect.l = gap;
    h_rect.x = pad_pts(3,1)-h_routing.w/2;
    h_rect.y = pad_pts(3,2) + route_offset(2);
    h_rect.layer = 6;
    h_rect.rounded = 0;
    gnd_cont_top = rect(h_rect);
    
    h_rect.w = h_routing.w;
    h_rect.l = gap;
    h_rect.x = pad_pts(3,1)-h_routing.w/2;
    h_rect.y = pad_pts(3,2) - route_offset(2)-gap;
    h_rect.layer = 6;
    h_rect.rounded = 0;
    gnd_cont_bot = rect(h_rect);
    
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