function gs = make_angled_arm(h_arm)
% Function to create the angled arm of the inchworm motors.
% The origin will be at the rightmost corner of the pawl, above the teeth

global joints

if ~isfield(h_arm,'rotation_center')
    h_arm.rotation_center = [0 0];
end

if ~isfield(h_arm,'theta')
    h_arm.rotation_theta = 0;
end


if ~isfield(h_arm,'VREP_group')
    h_arm.VREP_group = 'none';
end


% set tooth h_arm.orientation
toothOrient=1;

if h_arm.orientation ==1
    h_arm.pawlW=-h_arm.pawlW;
    h_arm.length=-h_arm.length;
    toothOrient=0;
    h_arm.theta=180-h_arm.theta;
end

% begin gds structure
gs=gds_structure(['angled_arm_' num2str(round(h_arm.pawlPos(1))) '_' num2str(round(h_arm.pawlPos(2)))]);

% form pawl block with teeth
xy1=[h_arm.pawlPos ; h_arm.pawlPos+[-h_arm.pawlL 0] ; h_arm.pawlPos+[-h_arm.pawlL h_arm.pawlW] ; h_arm.pawlPos+[0,h_arm.pawlW] ; h_arm.pawlPos];

% Add rotation
rot.pts = xy1;
rot.theta = h_arm.rotation_theta;
rot.p0 = h_arm.rotation_center;
xy1=rotate_pts(rot);

h_rect.x = h_arm.pawlPos(1);
h_rect.y = h_arm.pawlPos(2);
h_rect.w = -h_arm.pawlL;
h_rect.l = h_arm.pawlW;
h_rect.p0 =  h_arm.rotation_center;
h_rect.VREP_group = [h_arm.VREP_group  '_r'];
h_rect.theta = h_arm.rotation_theta;
h_rect.layer = h_arm.layer;
pawl_block = rect(h_rect);



%gs(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_arm.layer);

h_tooth.pos = h_arm.pawlPos;
h_tooth.toothW = h_arm.toothW;
h_tooth.toothL = h_arm.toothL;
h_tooth.toothS = h_arm.toothS;
h_tooth.Ntooth = h_arm.nToothPawl;
h_tooth.orientation  = toothOrient;
h_tooth.theta = h_arm.rotation_theta;
h_tooth.rotation_center = h_arm.rotation_center;
h_tooth.layer = h_arm.layer;
h_tooth.VREP_group = [h_arm.VREP_group  '_r'];
t_array=make_tooth_array(h_tooth);

%gs=join_gds_structures(gs,gsTemp);



% define points of the arm, angled by theta
pt1=h_arm.pawlPos+[-h_arm.pawlL,h_arm.pawlW];
pt2=h_arm.pawlPos+[-h_arm.pawlL+h_arm.width/sind(h_arm.theta),h_arm.pawlW];
pt3=h_arm.pawlPos+[-h_arm.pawlL+h_arm.width/sind(h_arm.theta)-h_arm.length*cosd(h_arm.theta),h_arm.pawlW+h_arm.length*sind(h_arm.theta)];
pt4=h_arm.pawlPos+[-h_arm.pawlL-h_arm.length*cosd(h_arm.theta),h_arm.pawlW+h_arm.length*sind(h_arm.theta)];
xy1=[pt1 ; pt2 ; pt3 ; pt4 ; pt1];

% Add rotation
rot.pts = xy1;
rot.theta = h_arm.rotation_theta;
rot.p0 = h_arm.rotation_center;
xy1=rotate_pts(rot);

%gs(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_arm.layer);

h_rect.x = xy1(2,1);
h_rect.y = xy1(2,2);
h_rect.w = -h_arm.width;
h_rect.l = h_arm.length + sign(h_arm.length)*h_arm.width;
h_rect.p0 =  xy1(2,:);
h_rect.VREP_group = [h_arm.VREP_group  '_r'];
h_rect.theta = atan2(xy1(2,2)-xy1(1,2),xy1(2,1)-xy1(1,1))+(pi/2-h_arm.theta*pi/180);
h_rect.layer = h_arm.layer;
angled_arm = rect(h_rect);
            
%g1 = [h_arm.VREP_group '_r'];
%g2 = h_rect.VREP_group;
%joints = VREP_add_joint(joints,g1,g2,'prismatic',xy1(3,:),h_arm.rotation_theta,[0 h_tooth.toothW]);


% Set respondable mask for the angled arm structure
%VREP_SetParam(g2,3019,2^0+2^1+2^2+2^3,0);
            
%% Grab all the GDS structures and arrays of structures

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
gs = b;


