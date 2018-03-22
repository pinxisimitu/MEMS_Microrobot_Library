function [out moto_pts]= motor(h_motor)
% Function to create an inchworm motor that can be rotated in any
% orientation.

% moto_pts(1,:) = end of shuttle (attachment point for serpentine spring
% moto_pts(2,:)

%% Default variables

% Number of sets of inchworms to break full motor into
if ~isfield(h_motor,'num_inch_sets')
    h_motor.num_inch_sets = 1;
end

% Width of the inchworm shuttle
if ~isfield(h_motor,'shuttle_w')
    h_motor.shuttle_w = 34;
end

% Ground connection Serpentines (1=they exist, 0=they are missing)
if ~isfield(h_motor,'ground_serpentine')
    h_motor.ground_serpentine = 0;
end

if ~isfield(h_motor,'ETCH_R')   %Default etch hole size (half length or radius)
    default_etch_properties;
    h_motor.ETCH_R = ETCH_R;
end

if ~isfield(h_motor,'UNDERCUT')   %Default etch hole undercut
    default_etch_properties;
    h_motor.UNDERCUT = UNDERCUT;
end

if ~isfield(h_motor,'layer')   %Default layers for structures and etch holes
    default_layer_properties;
    h_motor.layer = [SOI SOIHOLE];
end

if ~isfield(h_motor,'NOT_DUMMY')   %Default layers for structures and etch holes
    default_layer_properties;
    h_motor.NOTDUMMY = NOTDUMMY;
end

if ~isfield(h_motor,'hopper_shuttle_extend')
    h_motor.hopper_shuttle_extend = 0;
end


if ~isfield(h_motor,'circle_etch_holes')
    h_motor.circle_etch_holes = 0;
end

%%

% This global variable holds
global VREP_ignore s_motor joints

if VREP_ignore
    s_motor.motor_count = 0;
end

%Set etch hole flags
h_gca.ETCH_R = h_motor.ETCH_R;
h_gca.SOI_HOLE = h_motor.layer(2);
h_gca.UNDERCUT = h_motor.UNDERCUT;

% Set standard motor parameters for Tin Toy Motor
h_motor.toothW = 2.5000;
h_motor.toothL = 1.5000;
h_motor.toothS =  1.5000;
h_motor.shuttleG = 2;
h_motor.fingerW = 5;
h_motor.fingerL = 76.4000;
h_motor.fingerSupportL = 10;
h_motor.gap1 = 4.8000;
h_motor.gap2 = 7.7500;
h_motor.armW = 3;
h_motor.armL = 122;
h_motor.armAngle = 67.4000;
h_motor.springW = 3;
h_motor.springL = 240.8000;
h_motor.gstopGap = 3.8000;


% Output points that might be useful for downstream layout functions
moto_pts = [];          %motor_pts(1,:) = length of shuttle


%Need to rotate this point around the origin by negative h_motor.angle
%degrees
h_rot.pts =  h_motor.pos;
h_rot.theta = -pi/180*h_motor.angle;
h_rot.p0 = [0 0];
pos2 = rotate_pts(h_rot);                           % Position for teeth
pos =  h_motor.pos;

temp_label = ['Motor_' num2str(round(pos(1))) '_' num2str(round(pos(2)))];

toothW = h_motor.toothW;
toothL = h_motor.toothL;
toothS = h_motor.toothS;
shuttleG = h_motor.shuttleG;

fingerW = h_motor.fingerW;
fingerL = h_motor.fingerL;
fingerSupportL = h_motor.fingerSupportL;
gap1 = h_motor.gap1;
gap2 = h_motor.gap2;
N = h_motor.N;
travel = h_motor.travel;
armW = h_motor.armW;
armL = h_motor.armL;
armAngle = h_motor.armAngle;
springW = h_motor.springW;
springL = h_motor.springL;
gstopGap = h_motor.gstopGap;

% minimum anchor width
anchorW=50;

% default stucture spacing (spacing between structures, spacing around springs, etc.
space=5;

% routing width
routeW=30;

% stator anchor length
anchorL=N/2*(2*fingerW+gap1+gap2)-gap2+space+fingerW;

% gap stop parameters
% gapstopL, stop length
% gstopW, stop width
% gstopBumpW, stop bump width
% gstopBumpL, stop bump length
% gstopSpace, pitch between stop bumps
gstopL=110;
gstopW=10;
gstopBumpW=2;
gstopBumpL=5;
gstopSpace=15;

% mid-bar length that exceeds the anchor length, to give breadth for
% support springs on GCA
midBarOL=50;

% width and length of the rotor midbar, holding the rotor fingers
% default array has 2 rows of fingers
barW=20;
barL=N/2*(2*fingerW+gap1+gap2)+midBarOL+gstopW+gstopBumpW+anchorW;

% width and length for shuttle
% gap between pawls and shuttle
shuttleW=h_motor.shuttle_w;
shuttleL=4*anchorW+2*(fingerL+2*fingerSupportL)+2*(barW+springL)+4*space+3*routeW+travel;

moto_pts = [moto_pts; shuttleL,shuttleL];
%% Make the inchworm shuttle and teeth

% Create the shuttle
if h_motor.hopper_shuttle_extend
    extend_shuttle = 10;
else
    extend_shuttle = 0;
end

h_rect.x = pos(1) + extend_shuttle;
h_rect.y = pos(2) - shuttleW/2;
h_rect.w = -shuttleL-extend_shuttle;
h_rect.l = shuttleW;
h_rect.p0 = pos;
h_rect.theta = h_motor.angle *pi/180;
h_rect.layer = h_motor.layer;
h_rect.etch = 1;

if h_motor.circle_etch_holes
    h_rect.ETCH_R = 2;
    h_rect.circle_etch = 1; 
else
    h_rect.ETCH_R = 4;
    h_rect.circle_etch = 0;
end

h_rect.etch_layer = h_gca.SOI_HOLE;
h_rect.UNDERCUT = h_gca.UNDERCUT;
h_rect.VREP_group = sprintf('M%d_shuttle',s_motor.motor_count);
g1 = h_rect.VREP_group;
shuttle_re = rect(h_rect);
dimShuttle = [shuttleW shuttleL];

h_rect = rmfield(h_rect,'VREP_group');

if ~VREP_ignore
    % Add anchor for shuttle (for VREP)
    anch = 300;
    h_rect.x = pos(1)-shuttleL-anch - 50;
    h_rect.y = pos(2) - anch/2;
    h_rect.w = anch;
    h_rect.l = anch;
    h_rect.p0 = pos;
    h_rect.theta = h_motor.angle *pi/180;
    h_rect.layer = h_motor.layer;
    h_rect.etch = 0;
    h_rect.VREP_name = sprintf('M%d_shuttle_anchor',s_motor.motor_count);
    g2 = h_rect.VREP_name;
    shuttle_re_anchor = rect(h_rect);
    
    
    h_rect.rp = 1;
    ppts = rect(h_rect);
    h_rect.rp = 0;
    
    j_pos = midpt(ppts(2,:),ppts(3,:),.5);
    
    
    % Connect the shuttle
    name = sprintf('M%d_shuttle_j',s_motor.motor_count);
    joints = VREP_add_joint(joints,g1,g2,'prismatic',j_pos,h_motor.angle*pi/180,[0 1],name);
    
    h_rect = rmfield(h_rect,'VREP_name');
end

% Add dummy exclude layer around the shuttle
shuttle_dummy_gap = 5;
h_rect.x = pos(1);
h_rect.y = pos(2) - shuttleW/2 - shuttle_dummy_gap;
h_rect.w = -shuttleL;
h_rect.l = shuttleW+2*shuttle_dummy_gap;
h_rect.p0 = pos;
h_rect.theta = h_motor.angle *pi/180;
h_rect.layer = h_motor.NOTDUMMY;
h_rect.etch = 0;
shuttle_re_dummy = rect(h_rect);
dimShuttle = [shuttleW shuttleL];

% Add teeth to shuttle
nTooth = floor(shuttleL/(toothS+toothW));
h_tooth.pos = pos2 + [0 shuttleW/2];
h_tooth.toothW = toothW;
h_tooth.toothL = toothL;
h_tooth.toothS = toothS;
h_tooth.Ntooth = nTooth;
h_tooth.orientation = 0;
h_tooth.VREP_group = sprintf('M%d_shuttle',s_motor.motor_count);
h_tooth.theta = h_motor.angle *pi/180;
h_tooth.layer = h_motor.layer;

top_teeth = make_tooth_array(h_tooth);

h_tooth.pos = pos2 - [0 shuttleW/2];
h_tooth.orientation = 1;
bottom_teeth = make_tooth_array(h_tooth);


%% Define pawl parameters

nToothPawl = 2;                                    % Number of teeth on the pawl
pawlW=nToothPawl*toothW+toothS;                    % Width of pawl
pawlL=toothW*nToothPawl+toothS*(nToothPawl-1)+0.5; % Length of pawl

% define pawl origins
pawlX=anchorW+fingerL+fingerSupportL+barW/2-armW/sind(armAngle)/2-armL*cosd(armAngle)-pawlL-(pawlL-armW/sind(armAngle));
pawlX=pawlX-mod(pawlX,toothW+toothS)-toothW/2-toothS/2;
pawlY=dimShuttle(1)/2+2*toothL+shuttleG;

% define offset in the y-dimension for rear pawls
rearPawlX=anchorW+fingerL+fingerSupportL+barW+2*(springL+anchorW)+3*routeW+4*space+barW/2-armW/sind(armAngle)/2-armL*cosd(armAngle)-pawlL-(pawlL-armW/sind(armAngle));
sepErr=mod(rearPawlX,toothW+toothS);
spaceRoute=space-sepErr/4;
rearPawlX=rearPawlX-sepErr;

%% Generate GCAs and angled arms

% define GCA origin
gcaX=pawlX+armL*cosd(armAngle)+pawlL-armW/sind(armAngle)/2;
gcaY=pawlY+pawlW+armL*sind(armAngle);
reargcaX=rearPawlX+armL*cosd(armAngle)+pawlL-armW/sind(armAngle)/2;


%Make first set of drive (left/top side of shuttle)
x_shift = 548;      %How much space to put between motors (fix this todo)
gcaPos=pos-[gcaX-x_shift -gcaY];

num_inch_sets = h_motor.num_inch_sets;
%For jamming structures
jammer_len = 50;
jammer_dist = 4;

pts_top = [];
pts_bot = [];
for k =1:2
    % Generates front GCAs and angled arms (both above and below the shuttle)
    if k ==1
        for j = 1:num_inch_sets
            % Set initial position of the pawl
            pawlPos = pos - [pawlX+x_shift*(j-1) -pawlY];
            
            % Make the front top GCA's angled arm (jtg)
            h_arm.pawlPos = pawlPos;
            h_arm.width = armW;
            h_arm.length = armL;
            h_arm.nToothPawl = nToothPawl;
            h_arm.pawlW = pawlW;
            h_arm.pawlL = pawlL;
            h_arm.toothW = toothW;
            h_arm.toothL = toothL;
            h_arm.toothS = toothS;
            h_arm.theta = armAngle;
            h_arm.orientation = 0;
            h_arm.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j-1));
            h_arm.layer = h_motor.layer;
            h_arm.rotation_center = pos;
            h_arm.rotation_theta = h_motor.angle*pi/180;
            aa1 = make_angled_arm(h_arm);
            
            str = sprintf('aa1_%d_%d = aa1;',j,k);
            eval(str);
            
            % NEED TO ADD THE JOINTS IN!
            
            % Add joint at the angles arm attachment
            %g1 = [h_arm.VREP_group '_rotor'];
            %g2 = sprintf('M%d_ft_%d_rotor',s_motor.motor_count,j);
            %joints = VREP_add_joint(joints,g1,g2,'prismatic',h_arm.pawlPos,h_motor.angle*pi/180);
            
            
            %g1 = [h_arm.VREP_group '_r'];
            %g2 = h_rect.VREP_group;
            %joints = VREP_add_joint(joints,g1,g2,'prismatic',xy1(3,:),h_arm.rotation_theta,[0 h_tooth.toothW]);
            
            % Guide along shuttle
            guide_gap = 50;
            guide_w = 90;
            guide_round = 10;
            
            h_rect.x = h_arm.pawlPos(1)+guide_gap/2;
            h_rect.y = h_arm.pawlPos(2)-h_tooth.toothL;
            h_rect.w = x_shift - 3/2*guide_gap;
            h_rect.l = guide_w;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.VREP_group = sprintf('M%d_guides_anchor',s_motor.motor_count);
            h_rect.rounded = guide_round;
            
            if j>1  %Dont want a guide on the first pawl
                str = sprintf('guide_top_%d_%d = rect(h_rect);',j,k);
                eval(str);
            end
            
            % Add backstop to stop jamming of pawls
            h_rect.x = h_arm.pawlPos(1) - 7 - jammer_len - jammer_dist;
            h_rect.y = h_arm.pawlPos(2);
            h_rect.w = jammer_len;
            h_rect.l = 10;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.rounded = 0;
            h_rect.VREP_group = sprintf('M%d_guides_anchor',s_motor.motor_count);
            str = sprintf('backstop_anti_jammer_1_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = jammer_len-20;
            h_rect.l = 40;
            str = sprintf('backstop_anti_jammer_2_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = -200;
            h_rect.l = guide_w;
            str = sprintf('backstop_anti_jammer_12_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            % Front top GCA
            gcaPos = gcaPos - [x_shift 0];
            h_gca.pos = gcaPos;
            h_gca.width =fingerW;
            h_gca.length =fingerL;
            h_gca.lengthSupport =fingerSupportL;
            h_gca.gap1 =gap1;
            h_gca.gap2 = gap2;
            h_gca.nFingers = N;
            h_gca.springW = springW;
            h_gca.springL = springL;
            h_gca.gstopGap = gstopGap;
            h_gca.anchorW = anchorW;
            h_gca.space = space;
            h_gca.anchorL = anchorL;
            h_gca.ang = 90;
            h_gca.top = 1;
            h_gca.springOrient = 0;
            h_gca.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j-1));
            h_gca.rotation_center = pos;
            h_gca.rotation_theta = h_motor.angle *pi/180;
            h_gca.layers = [h_motor.layer h_motor.layer(2)];
            [GCA11 pts] = make_GCA_array_v2(h_gca);
            str = sprintf('GCA_%d_%d = GCA11;',j,k);
            eval(str);
            
            %Grab points from GCA Array
            rot.pts = h_gca.pos;
            rot.theta = h_gca.rotation_theta;
            rot.p0 = pos;
            gca_pos=rotate_pts(rot);
            moto_pts = [moto_pts; gca_pos];
            
            pts_top = [pts_top; pts];
            pawlPos=pos-[pawlX+x_shift*(j-1) pawlY];
            
            % Make the front bottom angled arm (jtg)
            h_arm.pawlPos = pawlPos;
            h_arm.orientation = 1;
            h_arm.VREP_group = h_gca.VREP_group;
            h_arm.rotation_center = pos;
            h_arm.rotation_theta = h_motor.angle *pi/180;
            h_arm.layer = h_motor.layer;
            %h_arm.VREP_group = sprintf('M%d_fb_%d',s_motor.motor_count,j);
            h_arm.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j));
            aarm2=make_angled_arm(h_arm);
            str = sprintf('aarm2_%d_%d = aarm2;',j,k);
            eval(str);
            
            % Guide along shuttle  (bottom)
            h_rect.x = h_arm.pawlPos(1)+guide_gap/2;
            h_rect.y = h_arm.pawlPos(2)+h_tooth.toothL-guide_w;
            h_rect.w = x_shift - 3/2*guide_gap;
            h_rect.l = guide_w;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.rounded = guide_round;
            
            if j>1  %Dont want a guide on the first pawl
                str = sprintf('guide_bot_%d_%d = rect(h_rect);',j,k);
                eval(str);
            end
            
            % Add backstop to stop jamming of pawls
            h_rect.x = h_arm.pawlPos(1) - 7 - jammer_len - jammer_dist;
            h_rect.y = h_arm.pawlPos(2);
            h_rect.w = jammer_len;
            h_rect.l = -10;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.rounded = 0;
            str = sprintf('backstop_anti_jammer_3_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = jammer_len-20;
            h_rect.l = -40;
            str = sprintf('backstop_anti_jammer_4_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = -200;
            h_rect.l = -guide_w;
            str = sprintf('backstop_anti_jammer_14_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            % Front bottom GCA
            gcaPos2=pos-[gcaX+x_shift*(j-1) gcaY];
            h_gca.pos = gcaPos2;
            h_gca.ang = 270;
            h_gca.top = 0;
            h_gca.rotation_center = pos;
            h_gca.rotation_theta = h_motor.angle *pi/180;
            h_gca.layers = [h_motor.layer h_motor.layer(2)];
            h_gca.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j));
            [GCAA22 pts] = make_GCA_array_v2(h_gca);
            str = sprintf('GCA22_%d_%d = GCAA22;',j,k);
            eval(str);
            
            
            %Grab points from GCA Array
            rot.pts = h_gca.pos;
            rot.theta = h_gca.rotation_theta;
            rot.p0 = pos;
            gca_pos=rotate_pts(rot);
            moto_pts = [moto_pts; gca_pos];
            
            
            h_gca.layers = h_motor.layer;
            h_gca.rotation_theta = 0;
            pts_top = [pts_top; pts];
            
        end
        x_shift = x_shift + 2;
        
    else
        % Generates the second (set of) top and bottom GCAs and arms
        %pawlPos=pos-[pawlX+x_shift*(num_inch_sets-1) -pawlY];   %Reset pawl position
        pawlPos = pawlPos + [0 2*pawlY];
        for j = 1:num_inch_sets
            if j==2
                x_shift = x_shift - 2;
            end
            
            %Make top rear angled arm (jtg)
            pawlPos=pawlPos - [x_shift 0];
            h_arm.pawlPos = pawlPos;
            h_arm.orientation = 0;
            h_arm.rotation_center = pos;
            %h_arm.VREP_group = sprintf('M%d_rt_%d',s_motor.motor_count,j);
            h_arm.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j-1));
            h_arm.rotation_theta = h_motor.angle *pi/180;
            h_arm.layer = h_motor.layer;
            aarm3=make_angled_arm(h_arm);
            str = sprintf('aarm3_%d_%d = aarm3;',j,k);
            eval(str);
            
            % Guide along shuttle
            h_rect.x = h_arm.pawlPos(1)+guide_gap/2;
            h_rect.y = h_arm.pawlPos(2)-h_tooth.toothL;
            h_rect.w = x_shift - 3/2*guide_gap;
            h_rect.l = guide_w;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.rounded = guide_round;
            str = sprintf('guide2_top_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            % Add backstop to stop jamming of pawls
            h_rect.x = h_arm.pawlPos(1) - 7 - jammer_len - jammer_dist;
            h_rect.y = h_arm.pawlPos(2);
            h_rect.w = jammer_len;
            h_rect.l = 10;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.rounded = 0;
            str = sprintf('backstop_anti_jammer_7_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = jammer_len-20;
            h_rect.l = 40;
            str = sprintf('backstop_anti_jammer_8_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = -200;
            h_rect.l = guide_w;
            str = sprintf('backstop_anti_jammer_18_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            
            
            gcaPos=gcaPos-[x_shift 0];
            
            % Rear top GCA
            h_gca.pos = gcaPos;
            h_gca.ang = 90;
            h_gca.top = 1;
            h_gca.springOrient = 0;
            h_gca.rotation_center = pos;
            h_gca.rotation_theta = h_motor.angle *pi/180;
            h_gca.layers = [h_motor.layer h_motor.layer(2)];
            h_gca.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j-1));
            [gca2322 pts] = make_GCA_array_v2(h_gca);
            str = sprintf('GCA23_%d_%d = gca2322;',j,k);
            eval(str);
            
            
            %Grab points from GCA Array
            rot.pts = h_gca.pos;
            rot.theta = h_gca.rotation_theta;
            rot.p0 = pos;
            gca_pos=rotate_pts(rot);
            moto_pts = [moto_pts; gca_pos];
            
            
            
            pts_bot = [pts_bot; pts];
            
            pp2y = pos(2) - pawlY;
            
            % Second angled arm (jtg)
            h_arm.pawlPos = [pawlPos(1) pp2y];
            h_arm.orientation = 1;
            h_arm.layer = h_motor.layer;
            h_arm.rotation_center = pos;
            %h_arm.VREP_group = sprintf('M%d_rb_%d',s_motor.motor_count,j);
            h_arm.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j));
            h_arm.rotation_theta = h_motor.angle *pi/180;
            aaaarm2j=make_angled_arm(h_arm);
            str = sprintf('aaaarm2j_%d_%d = aaaarm2j;',j,k);
            eval(str);
            
            % Guide along shuttle  (bottom)
            h_rect.x = h_arm.pawlPos(1)+guide_gap/2;
            h_rect.y = h_arm.pawlPos(2)+h_tooth.toothL-guide_w;
            h_rect.w = x_shift - 3/2*guide_gap;
            h_rect.l = guide_w;
            h_rect.etch = 0;
            h_rect.rounded = guide_round;
            h_rect.layer = h_motor.layer(1);
            
            str = sprintf('guide_bot22_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            
            % Add backstop to stop jamming of pawls
            h_rect.x = h_arm.pawlPos(1) - 7 - jammer_len - jammer_dist;
            h_rect.y = h_arm.pawlPos(2);
            h_rect.w = jammer_len;
            h_rect.l = -10;
            h_rect.etch = 0;
            h_rect.layer = h_motor.layer(1);
            h_rect.rounded = 0;
            str = sprintf('backstop_anti_jammer_5_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = jammer_len-20;
            h_rect.l = -40;
            str = sprintf('backstop_anti_jammer_6_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            h_rect.w = -200;
            h_rect.l = -guide_w;
            str = sprintf('backstop_anti_jammer_16_%d_%d = rect(h_rect);',j,k);
            eval(str);
            
            
            
            gcaPos2=gcaPos2-[x_shift 0];
            
            % Rear bottom GCA
            h_gca.pos = gcaPos2;
            h_gca.ang = 270;
            h_gca.rotation_theta = h_motor.angle *pi/180;
            h_gca.layers = [h_motor.layer h_motor.layer(2)];
            h_gca.top = 0;
            h_gca.springOrient = 0;
            h_gca.VREP_group = sprintf('M%d_%d',s_motor.motor_count,(k-1)*(2*num_inch_sets)+(2*j));
            [gca22ds2 pts] = make_GCA_array_v2(h_gca);
            str = sprintf('GCA23ds_%d_%d = gca22ds2;',j,k);
            eval(str);
            pts_bot = [pts_bot; pts];
            
            %Grab points from GCA Array
            rot.pts = h_gca.pos;
            rot.theta = h_gca.rotation_theta;
            rot.p0 = pos;
            gca_pos=rotate_pts(rot);
            moto_pts = [moto_pts; gca_pos];
            
            
        end
    end
end
clear gca22ds2 gca2322 GCAA22 GCA11
clear aaaarm2j aarm3 aarm2 aa1

s_motor.motor_count = s_motor.motor_count + 1;

%% Add in serpentine springs along the ground contacts

%temp fix

% Add top most serpentine to route grounds together
add_serp = 1;


additional_serpentine_top = [pts_top(1:4,:)-[pts_top(5,1)-pts_top(1,1) 0]; pts_top];
additional_serpentine_bot = [pts_bot(1:4,:)-[pts_bot(5,1)-pts_bot(1,1) 0];pts_bot(1:4,:)-[pts_bot(5,1)-pts_bot(1,1) 0]; pts_bot];

pts_top = [pts_top;pts_bot];

count = 0;
if h_motor.ground_serpentine
    ss_ground = gds_structure(['Serpentines_' num2str(round(pts_bot(1))) '_' num2str(round(pts_bot(2)))]);
    for serp = 0:add_serp
        for kk = 0:1
            for i = 1:2*num_inch_sets-1
                if(serp == 1 && kk==0 && i == 1)
                    pts_top = additional_serpentine_top;
                end
                %Make rectangels on both ends of serpentine beam
                ss_w = 3;
                p0 = pts_top((i-1)*4+2+kk*2,:);
                p1 = midpt(pts_top((i-1)*4+2+kk*2,:),pts_top((i-1)*4+5+kk*2,:),.1);
                p2 = midpt(pts_top((i-1)*4+2+kk*2,:),pts_top((i-1)*4+5+kk*2,:),.5);
                p3 = pts_top((i-1)*4+5+kk*2,:);
                
                %Rotate all of those points
                rot.pts = [p0;p1;p2;p3];
                rot.theta = h_motor.angle *pi/180;
                rot.p0 = pos;
                pts=rotate_pts(rot);
                p0 = pts(1,:);
                p1 = pts(2,:);
                p2 = pts(3,:);
                p3 = pts(4,:);
                
                
                %Connect to top of serpentine ground contact
                be = gds_element('path', 'xy',[p0;p1],'width', ss_w,'layer',h_motor.layer);
                str_name = sprintf('SS_a_[%d,%d],[%d,%d]',round(p1(1)),round(p1(2)),round(p2(1)),round(p2(2)));
                conn1 = gds_structure(str_name,be);
                
                be = gds_element('path', 'xy',[p2;p3],'width', ss_w,'layer',h_motor.layer);
                str_name = sprintf('SS_b_[%d,%d],[%d,%d]',round(p1(1)),round(p1(2)),round(p2(1)),round(p2(2)));
                conn2 = gds_structure(str_name,be);
                
                ss_ground=join_gds_structures(ss_ground,conn1);
                ss_ground=join_gds_structures(ss_ground,conn2);
                
                h_ss.p1 = p1;     % First point that the SS will span from
                h_ss.p2 = p2;    % 20 is the width of the shuttle
                h_ss.n = 3;                                       % Number of meanders
                h_ss.w = 3;                                      % Width of beams
                h_ss.dpp = 40;                                   % Peak to peak distance of meanders
                h_ss.layer = h_motor.layer;                                   % Layer
                spring_temp = s_spring(h_ss);
                ss_ground=join_gds_structures(ss_ground,spring_temp);
                
                %Add dummy fill for these springs
                h_path.pts = [p0;p3];       % Make sure each row is one point [x,y]
                
                h_path.w = h_ss.dpp*1.4;                           % Width of path
                h_path.layer = h_motor.NOTDUMMY;                     % Set layer of path
                dummy_path = m_path(h_path);               % Function to create the path GDS structure
                
                ss_ground=join_gds_structures(ss_ground,dummy_path{1});
                
                clear spring_temp conn1 conn2 dummy_path
                
                if(serp == 1 && i==2)
                    break;
                end
            end
        end
    end
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