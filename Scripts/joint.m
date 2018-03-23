function out = joint(h_joint)

%%
if ~isfield(h_joint,'EH_layer') % Layer for the etch holes
    default_layer_properties;
    h_joint.EH_layer = SOIHOLE;
end

if ~isfield(h_joint,'n')        % Number of points in pin
    h_joint.n = round(h_joint.r(2)*2*pi); 
end

if ~isfield(h_joint,'gap')      % Gap between the pin and cuff 
    h_joint.gap = 3; 
end

if ~isfield(h_joint,'overlap')      % Overlap of the inner arm and the pin 
    h_joint.overlap = 6; 
end

if ~isfield(h_joint,'opening_theta_center')  % angle around which the opening in the joint will be centered 
    h_joint.opening_theta_center = h_joint.inner_arm_theta; 
end

if ~isfield(h_joint,'ss_cw')      % Generates the serpentine spring clockwise from the inner arm  
    h_joint.ss_cw = 1; 
end

if ~isfield(h_joint,'UNDERCUT')     
    default_etch_properties
    h_joint.UNDERCUT = UNDERCUT; 
end

if ~isfield(h_joint,'ETCH_R')     
    default_etch_properties
    h_joint.ETCH_R = ETCH_R; 
end

if ~isfield(h_joint,'CIRCULAR_ETCH')     
    default_etch_properties
    h_joint.CIRCULAR_ETCH = CIRCULAR_ETCH; 
end


%%

% Create the C that goes around the ball of the joint
points = zeros(2*(h_joint.n+1),2);
theta_gap = h_joint.opening_theta*pi/180;
%theta_init = h_joint.inner_arm_theta*pi/180 + theta_gap/2; 
theta_init = h_joint.opening_theta_center*pi/180 + theta_gap/2; 


for j=1:2
    if j==2
        mult = -1;
        %theta_init = h_joint.inner_arm_theta*pi/180 - theta_gap/2;
        theta_init = h_joint.opening_theta_center*pi/180 - theta_gap/2; 
    else
        mult = 1;
    end
    for i=1:h_joint.n+1
        points(i + (j-1)*(h_joint.n+1),1) = h_joint.p0(1) + (h_joint.r(j)+.05)*cos(theta_init+mult*(i-1)*(2*pi-theta_gap)/h_joint.n);
        points(i + (j-1)*(h_joint.n+1),2) = h_joint.p0(2) + (h_joint.r(j)+.05)*sin(theta_init+mult*(i-1)*(2*pi-theta_gap)/h_joint.n);
    end
end
rotor_points = [points(end,:);points];

% Calculate points for etch holes (only the bottom half of the outer circle
rotor_points2 = zeros(h_joint.n,2);
for i=1:h_joint.n+1
    rotor_points2(i,1) = h_joint.p0(1) + (h_joint.r(2)+.05)*cos(pi+i*pi/h_joint.n);
    rotor_points2(i,2) = h_joint.p0(2) + (h_joint.r(2)+.05)*sin(pi+i*pi/h_joint.n);
end


outer_rotor_ring = rotor_points(h_joint.n+1:end,:);

be = gds_element('boundary', 'xy',rotor_points,'layer',h_joint.layer);
str_name = sprintf('Joint_C_[%d,%d],[%d]',round(h_joint.p0(1)),round(h_joint.p0(2)),round(h_joint.r(1)));
c_joint = gds_structure(str_name,be);

%Make etch holes in rotor C
h_etch.regions = cell(1,1);
h_etch.r =  h_joint.ETCH_R;
h_etch.undercut = h_joint.UNDERCUT;
h_etch.layer = h_joint.EH_layer;
section.p0 = [h_joint.p0(1), h_joint.p0(2)];
section.type = 'partial_annulus';
section.r = h_joint.r;

%section.theta = [h_joint.inner_arm_theta*pi/180 + theta_gap/2 + (h_etch.undercut+h_etch.r)/h_joint.r(2), 2*pi-theta_gap-(h_etch.undercut+h_etch.r)/h_joint.r(2)]; 
section.theta = [h_joint.opening_theta_center*pi/180 + theta_gap/2 + (h_etch.undercut+h_etch.r)/h_joint.r(2), 2*pi-theta_gap-(h_etch.undercut+h_etch.r)/h_joint.r(2)]; 


h_etch.regions{1,1} = section;
C_holes = etch_hole(h_etch);

%Creates Pin
h_cir.x = h_joint.p0(1);
h_cir.y = h_joint.p0(2);
h_cir.r = h_joint.r(1)-h_joint.gap;
h_cir.layer = h_joint.layer;
h_cir.n = 100;
joint_ball = circle(h_cir);
h_cir.rp = 1;
pin_points = circle(h_cir);
pin_points_n = h_cir.n;
h_cir.rp = 0;


%Add etch holes to pin
h_etch.regions = cell(1,1);
h_etch.r = h_joint.ETCH_R;
h_etch.undercut = h_joint.UNDERCUT;
h_etch.layer = h_joint.EH_layer;
section.p0 = [h_joint.p0(1), h_joint.p0(2)];
section.type = 'annulus';
section.r = [0 h_joint.r(1)-h_joint.gap];
section.theta = [h_joint.inner_arm_theta*pi/180 + theta_gap/2 + (h_etch.undercut+h_etch.r)/h_joint.r(2), 2*pi-theta_gap-(h_etch.undercut+h_etch.r)/h_joint.r(2)]; 
h_etch.regions{1,1} = section;
ball_holes = etch_hole(h_etch);

% Create inner arm
h_rect.x = h_joint.p0(1);
h_rect.y = h_joint.p0(2)-h_joint.inner_arm_w/2;
h_rect.w = h_joint.inner_arm_l;
h_rect.l = h_joint.inner_arm_w;
h_rect.p0 = h_joint.p0;
h_rect.theta = h_joint.inner_arm_theta*pi/180;
h_rect.layer = h_joint.layer;
inner_arm = rect(h_rect);

% Create dummy fill for inner arm
% df_x_l = 1000;
% df_x_r = 200;
% df_y_b = 200;
% df_y_t = 500;
% h_rect.x = h_joint.p0(1) - df_x_l;
% h_rect.y = h_joint.p0(2)-h_joint.inner_arm_w/2 - df_y_b;
% h_rect.w = h_joint.inner_arm_l + df_x_l + df_x_r;
% h_rect.l = h_joint.inner_arm_w + df_y_b + df_y_t;
% h_rect.p0 = h_joint.p0;
% h_rect.layer = 8;
% h_rect.theta = h_joint.inner_arm_theta*pi/180;

%df_joint = rect(h_rect);

% % Serpentine Spring and Contacts
% Add serpentine spring contact onto inner arm
% ss_contact_w = 6;
% ss_w = 4;
% ss_contact_l_s = 60;
% ss_contact_l_l = 110;
% ss_contact_dist_from_c = 70;
% 
% h_joint.inner_arm_theta = 0;
% 
% if h_joint.ss_cw    %Generates the serpentine spring clockwise from the inner arm
%     
%     h_rect.x = h_joint.p0(1) + h_joint.r(2) + ss_contact_dist_from_c;
%     h_rect.y = h_joint.p0(2) - h_joint.inner_arm_w/2;
%     h_rect.w = ss_contact_w;
%     h_rect.layer = h_joint.layer;
%     h_rect.l = -ss_contact_l_l;
%     h_rect.p0 = h_joint.p0;
%     h_rect.theta = h_joint.inner_arm_theta*pi/180;
%     h_rect.theta = 0;
%     inner_arm_ss_contact = rect(h_rect);
%     
%     This is the contact (and base) for one side of the spring
%     cont2 = [h_rect.x+ss_contact_w/2 h_rect.y-ss_contact_l_l];
%     h_rot.pts = cont2;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont2=rotate_pts(h_rot);
%     
%     A shifted version of contact 2 so the serpentine overlaps exactly
%     cont2_shifted = [h_rect.x+ss_contact_w h_rect.y-ss_contact_l_l];
%     h_rot.pts = cont2_shifted;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont2_shifted=rotate_pts(h_rot);
%     
%     base2 = [h_rect.x+ss_contact_w/2 h_rect.y];
%     h_rot.pts = base2;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     base2=rotate_pts(h_rot);
%     
%     
%     Add fillets to the contact on inner arm
%     Right fillet
%     snap = .01;
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x+ss_contact_w h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 -y_fill_len];
%     h_fillet.p2 = h_fillet.p0 + [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_joint_r = fillet(h_fillet);
%     
%     Left fillet
%     snap = .01;
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 -y_fill_len];
%     h_fillet.p2 = h_fillet.p0 - [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_joint_l = fillet(h_fillet);
%     
%     Add serpentine spring contact onto C joint
%     snap = .5;
%     Contact on C
%     h_rect.x = h_joint.p0(1) - ss_contact_w/2;
%     h_rect.y = h_joint.p0(2) - h_joint.r(2) + snap;
%     h_rect.w = ss_contact_w;
%     h_rect.l = -ss_contact_l_s;
%     h_rect.p0 = h_joint.p0;
%     h_rect.theta = h_joint.inner_arm_theta*pi/180;
%     h_rect.theta = 0;
%     C_ss_contact = rect(h_rect);
%     
%     
%     This is the contact for the other side of the spring
%     cont1 = [h_rect.x+ss_contact_w/2 h_rect.y-ss_contact_l_s];
%     h_rot.pts = cont1;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont1=rotate_pts(h_rot);
%     
%     Shifted version to align the serpentine springs exactly
%     cont1_shifted = [h_rect.x h_rect.y-ss_contact_l_s];
%     h_rot.pts = cont1_shifted;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont1_shifted=rotate_pts(h_rot);
%     
%     base1 = [h_rect.x+ss_contact_w/2 h_rect.y];
%     h_rot.pts = base1;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     base1=rotate_pts(h_rot);
%     
%     Add fillets to the contact on the C
%     Right fillet
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x+ss_contact_w h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 -y_fill_len];
%     h_fillet.p2 = h_fillet.p0 + [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_C_joint_r = fillet(h_fillet);
%     
%     Left fillet
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 -y_fill_len];
%     h_fillet.p2 = h_fillet.p0 - [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_C_joint_l = fillet(h_fillet);
%     
%     Add serpentine spring conncting the two
%     h_ss.p1 = cont1_shifted;
%     h_ss.p2 = cont2_shifted;
%     h_ss.p1 = cont1 + [ss_contact_w/2 ss_w];
%     h_ss.p2 = cont2 + [-ss_contact_w/2 -ss_w/2];
%     h_ss.n = 4;
%     h_ss.w = ss_w;
%     h_ss.dpp = 70; %was 40 in the designs that didn't work
%     h_ss.layer = h_joint.layer;
%     serpentine = s_spring(h_ss);
%     h_ss.rp = 1;
%     ss_points = s_spring(h_ss);
%     
%     Now smooth out all the points
%     final_ss_points = [base1;ss_points;base2]';
%     
%     init_points = final_ss_points;
%     
%     interm_pts = fnplt(cscvn(init_points));
%     
%     Don't do any downsampling
%     final_pts = [interm_pts(1,1) downsample(interm_pts(1,2:end-1),2) interm_pts(1,end);
%        interm_pts(2,1) downsample(interm_pts(2,2:end-1),2) interm_pts(2,end)];
%     
%     points = [interm_pts'];
%     temp = gds_element('path', 'xy',points,'width',ss_contact_w,'layer',h_fillet.layer);
%     str_name = sprintf('SS_Spline_[%d,%d]',round(cont1(1)),round(cont1(2)));
%     serpentine_spline_l = gds_structure(str_name,temp);
%     
%     Create DC rotary springs
%     syms x y
%         
%     Set springs to be generated on correct side of inner/outer arms
%     if h_joint.inner_arm_theta > h_joint.outer_arm_theta
%         p0_inner = h_joint.p0 - h_joint.inner_arm_w/2*[-sind(h_joint.inner_arm_theta)  cosd(h_joint.inner_arm_theta)];
%         p0_outer = h_joint.p0 - h_joint.outer_arm_w/2*[sind(h_joint.outer_arm_theta)  -cosd(h_joint.outer_arm_theta)];
%     else
%         p0_inner = h_joint.p0 + h_joint.inner_arm_w/2*[-sind(h_joint.inner_arm_theta)  cosd(h_joint.inner_arm_theta)];
%         p0_outer = h_joint.p0 + h_joint.outer_arm_w/2*[sind(h_joint.outer_arm_theta)  -cosd(h_joint.outer_arm_theta)];      
%     end
% 
%     if h_joint.outer_arm_theta == 90
%        [sx sy] = solve(y-p0_inner(2) == tand(h_joint.inner_arm_theta)*(x-p0_inner(1)),x == p0_outer(1));
%     
%     elseif h_joint.outer_arm_theta == 270
%        [sx sy] = solve(y-p0_inner(2) == tand(h_joint.inner_arm_theta)*(x-p0_inner(1)),x == p0_outer(1));
%     
%     elseif h_joint.inner_arm_theta == 90
%            [sx sy] = solve(x == p0_inner(1),y-p0_outer(2) == tand(h_joint.outer_arm_theta)*(x-p0_outer(1)));
%       
%     elseif h_joint.inner_arm_theta == 270
%            [sx sy] = solve(x == p0_inner(1),y-p0_outer(2) == tand(h_joint.outer_arm_theta)*(x-p0_outer(1)));
%     else
%         Neither arm is vertical
%         [sx sy] = solve(y-p0_inner(2) == tand(h_joint.inner_arm_theta)*(x-p0_inner(1)),y-p0_outer(2) == tand(h_joint.outer_arm_theta)*(x-p0_outer(1)));
%     end
%     
%     h_spring.pos = [double(sx) double(sy)];
%     h_spring.N=3;
%     h_spring.l=200;
%     h_spring.radius=140;
%     h_spring.theta1 = h_joint.inner_arm_theta;
%     h_spring.theta2 = h_joint.outer_arm_theta;
%     
%     Make sure all angles are positive
%     if h_spring.theta1 < 0
%         h_spring.theta1 = h_spring.theta1+360;
%     end
%     if h_spring.theta2 < 0
%         h_spring.theta2 = h_spring.theta1+360;
%     end
%     
%     rotary_spring=make_rot_spring(h_spring);
% 
% else
%     Inner arm serpentine spring contact
%     h_rect.x = h_joint.p0(1) + h_joint.r(2) + ss_contact_dist_from_c;
%     h_rect.y = h_joint.p0(2) + h_joint.inner_arm_w/2;
%     h_rect.w = ss_contact_w;
%     h_rect.layer = h_joint.layer;
%     h_rect.l = ss_contact_l_l;
%     h_rect.p0 = h_joint.p0;
%     h_rect.theta = h_joint.inner_arm_theta*pi/180;
%     h_rect.theta = 0;
%     inner_arm_ss_contact = rect(h_rect);
%     
%     Inner arm contact point for serpentine spring
%     cont2_shifted = [h_rect.x+ss_contact_w h_rect.y+ss_contact_l_l];
%     h_rot.pts = cont2_shifted;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont2_shifted=rotate_pts(h_rot);
%     
%     base2 = [h_rect.x+ss_contact_w/2 h_rect.y];
%     h_rot.pts = base2;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     base2=rotate_pts(h_rot);
%     
%     
%     Add fillets to the contact on inner arm
%     Right fillet
%     snap = .01;
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x+ss_contact_w h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 y_fill_len];
%     h_fillet.p2 = h_fillet.p0 + [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_joint_r = fillet(h_fillet);
%     
%     Left fillet
%     snap = .01;
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 y_fill_len];
%     h_fillet.p2 = h_fillet.p0 - [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_joint_l = fillet(h_fillet);
%     
%     Add serpentine spring contact onto C joint
%     snap = .5;
%     Contact on C
%     h_rect.x = h_joint.p0(1) - ss_contact_w/2;
%     h_rect.y = h_joint.p0(2) + h_joint.r(2) - snap;
%     h_rect.w = ss_contact_w;
%     h_rect.l = ss_contact_l_s;
%     h_rect.p0 = h_joint.p0;
%     h_rect.theta = h_joint.inner_arm_theta*pi/180;
%     h_rect.theta = 0;
%     C_ss_contact = rect(h_rect);
%     
%     
%     This is the contact for the other side of the spring
%     cont1 = [h_rect.x+ss_contact_w/2 h_rect.y-ss_contact_l_s];
%     h_rot.pts = cont1;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont1=rotate_pts(h_rot);
%     
%     Shifted version to align the serpentine springs exactly
%     cont1_shifted = [h_rect.x h_rect.y+ss_contact_l_s];
%     h_rot.pts = cont1_shifted;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     cont1_shifted=rotate_pts(h_rot);
%     
%     base1 = [h_rect.x+ss_contact_w/2 h_rect.y];
%     h_rot.pts = base1;
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     base1=rotate_pts(h_rot);
%     
%     Add fillets to the contact on the C
%     Right fillet
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x+ss_contact_w h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 y_fill_len];
%     h_fillet.p2 = h_fillet.p0 + [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_C_joint_r = fillet(h_fillet);
%     
%     Left fillet
%     y_fill_len = 10;
%     x_fill_len = 6;
%     h_fillet.d = .6;
%     h_fillet.layer = h_joint.layer;
%     h_fillet.p0 = [h_rect.x h_rect.y];
%     
%     h_fillet.p1 = h_fillet.p0 + [0 y_fill_len];
%     h_fillet.p2 = h_fillet.p0 - [x_fill_len 0];
%     
%     Rotate P0,P1, and P2
%     h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
%     h_rot.theta = h_joint.inner_arm_theta*pi/180;
%     h_rot.p0 = h_joint.p0;
%     out=rotate_pts(h_rot);
%     
%     h_fillet.p0 = out(1,:);
%     h_fillet.p1 = out(2,:);
%     h_fillet.p2 = out(3,:);
%     
%     fill_C_joint_l = fillet(h_fillet);
%     
%     Add serpentine spring conncting the two
%     h_ss.p1 = cont1_shifted;
%     h_ss.p2 = cont2_shifted;
%     h_ss.p1 = cont1 + [ss_contact_w/2 ss_w];
%     h_ss.p2 = cont2 + [-ss_contact_w/2 -ss_w/2];
%     h_ss.n = 4;
%     h_ss.w = ss_w;
%     h_ss.dpp = 70; %was 40 in the designs that didn't work
%     h_ss.layer = h_joint.layer;
%     serpentine = s_spring(h_ss);
%     h_ss.rp = 1;
%     ss_points = s_spring(h_ss);
%     
%     Now smooth out all the points
%     final_ss_points = [base1;ss_points;base2]';
%     
%     init_points = final_ss_points;
%     
%     interm_pts = fnplt(cscvn(init_points));
%     
%     Don't do any downsampling
%     final_pts = [interm_pts(1,1) downsample(interm_pts(1,2:end-1),2) interm_pts(1,end);
%        interm_pts(2,1) downsample(interm_pts(2,2:end-1),2) interm_pts(2,end)];
%     
%     points = [interm_pts'];
%     temp = gds_element('path', 'xy',points,'width',ss_contact_w,'layer',h_fillet.layer);
%     str_name = sprintf('SS_Spline_[%d,%d]',round(cont1(1)),round(cont1(2)));
%     serpentine_spline_l = gds_structure(str_name,temp);
% 
% end

%% Add DC's serpentine springs
% Create DC rotary springs
syms x y

% Set springs to be generated on correct side of inner/outer arms
if h_joint.inner_arm_theta > h_joint.outer_arm_theta
    p0_inner = h_joint.p0 - h_joint.inner_arm_w/2*[-sind(h_joint.inner_arm_theta)  cosd(h_joint.inner_arm_theta)];
    p0_outer = h_joint.p0 - h_joint.outer_arm_w/2*[sind(h_joint.outer_arm_theta)  -cosd(h_joint.outer_arm_theta)];
else
    p0_inner = h_joint.p0 + h_joint.inner_arm_w/2*[-sind(h_joint.inner_arm_theta)  cosd(h_joint.inner_arm_theta)];
    p0_outer = h_joint.p0 + h_joint.outer_arm_w/2*[sind(h_joint.outer_arm_theta)  -cosd(h_joint.outer_arm_theta)];
end

if h_joint.outer_arm_theta == 90
    [sx sy] = solve(y-p0_inner(2) == tand(h_joint.inner_arm_theta)*(x-p0_inner(1)),x == p0_outer(1));
    
elseif h_joint.outer_arm_theta == 270
    [sx sy] = solve(y-p0_inner(2) == tand(h_joint.inner_arm_theta)*(x-p0_inner(1)),x == p0_outer(1));
    
elseif h_joint.inner_arm_theta == 90
    [sx sy] = solve(x == p0_inner(1),y-p0_outer(2) == tand(h_joint.outer_arm_theta)*(x-p0_outer(1)));
    
elseif h_joint.inner_arm_theta == 270
    [sx sy] = solve(x == p0_inner(1),y-p0_outer(2) == tand(h_joint.outer_arm_theta)*(x-p0_outer(1)));
else
    % Neither arm is vertical
    [sx sy] = solve(y-p0_inner(2) == tand(h_joint.inner_arm_theta)*(x-p0_inner(1)),y-p0_outer(2) == tand(h_joint.outer_arm_theta)*(x-p0_outer(1)));
end

h_spring.pos = [double(sx) double(sy)];
h_spring.N=3;
h_spring.l=200;
h_spring.radius=140;
h_spring.layer = h_joint.layer;
h_spring.theta1 = h_joint.inner_arm_theta;
h_spring.theta2 = h_joint.outer_arm_theta;

% Make sure all angles are positive
if h_spring.theta1 < 0
    h_spring.theta1 = h_spring.theta1+360;
end
if h_spring.theta2 < 0
    h_spring.theta2 = h_spring.theta1+360;
end

rotary_spring=make_rot_spring(h_spring);
    

%%

% Add etch holes to inner arm
h_etch.regions = cell(1,1);
section.type = 'tcurve';                %Pretend the circular curve is at the top and rotate
section.p0 = h_joint.p0 - [h_joint.inner_arm_w/2  h_joint.inner_arm_l];
section.tcurve = pin_points(round(pin_points_n/2):end,:);

h_etch.r = h_joint.ETCH_R;
h_etch.undercut = h_joint.UNDERCUT;
h_etch.w = h_joint.inner_arm_w;
h_etch.layer = h_joint.EH_layer;
h_etch.regions{1,1} = section;
h_etch.rotation_angle = -90 - h_joint.inner_arm_theta;
h_etch.rotation_point = h_joint.p0;
inner_arm_etch_holes = etch_hole(h_etch);

% Create outer arm on the C
h_rect.x = h_joint.p0(1) + h_joint.r(1);
h_rect.y = h_joint.p0(2)-h_joint.outer_arm_w/2;
h_rect.w = h_joint.outer_arm_l;
h_rect.l = h_joint.outer_arm_w;
h_rect.p0 = h_joint.p0;
h_rect.theta = h_joint.outer_arm_theta*pi/180;

outer_arm = rect(h_rect);

% Add etch holes to outer arm
h_etch.regions = cell(1,1);
section.type = 'tcurve';                %Pretend the circular curve is at the top and rotate
section.p0 = h_joint.p0 - [h_joint.outer_arm_w/2  h_joint.outer_arm_l+h_joint.r(1)];
section.tcurve = rotor_points(round(1.5*h_joint.n):end,:);
section.tcurve = rotor_points(h_joint.n:round(1.5*h_joint.n),:);
section.tcurve = rotor_points2;
%section.tcurve = outer_rotor_ring;
%section.tcurve(:,2) = section.tcurve(:,2) - h_joint.r(1); 

h_etch.r = h_joint.ETCH_R;
h_etch.undercut = h_joint.UNDERCUT;
h_etch.w = h_joint.outer_arm_w;
h_etch.regions{1,1} = section;
h_etch.rotation_angle = -90 - h_joint.outer_arm_theta;
h_etch.rotation_point = h_joint.p0;
outer_arm_etch_holes = etch_hole(h_etch);


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