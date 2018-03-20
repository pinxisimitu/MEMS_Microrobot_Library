function out = iono_hole(h_ihole)
% function to make a new ionocraft hole for capillary tubes

% Set defaults
if ~isfield(h_ihole,'SOI')
    h_ihole.SOI = 6;
end

if ~isfield(h_ihole,'SOI_HOLE')
    h_ihole.SOI_HOLE = 2;
end

if ~isfield(h_ihole,'EH_UNDERCUT')
    h_ihole.EH_UNDERCUT = 6;
end

if ~isfield(h_ihole,'rot_ang')
    h_ihole.rot_ang = 0;
end


if ~isfield(h_ihole,'etch_r')
    h_ihole.etch_r = 2;
end


if ~isfield(h_ihole,'theta_1')
    h_ihole.theta_1 = 0;                        % Initial angle of annulus
end

if ~isfield(h_ihole,'theta_2')
    h_ihole.theta_2 = pi;                       % Final angle of annulus
end



%To rotate points
h_rot.p0 = h_ihole.p0;
h_rot.theta = h_ihole.rot_ang-pi/2;

% Make half annulus
num_pts = round(h_ihole.radii(1)/2);               % How many points to be in arc
points = zeros(2*(num_pts+1),2);    % Matrix to store annulus pts
for j=1:2
    if j==2
        mult = -1;
        theta_temp = h_ihole.theta_2;
        h_ihole.theta_2 = h_ihole.theta_1;
        h_ihole.theta_1 = theta_temp;
    else
        mult = 1;
    end
    for i=1:num_pts+1
        points(i + (j-1)*(num_pts+1),1) = h_ihole.p0(1) + (h_ihole.radii(j)+.05)*cos(h_ihole.theta_1+mult*(i-1)*abs(h_ihole.theta_2-h_ihole.theta_1)/num_pts);
        points(i + (j-1)*(num_pts+1),2) = h_ihole.p0(2) + (h_ihole.radii(j)+.05)*sin(h_ihole.theta_1+mult*(i-1)*abs(h_ihole.theta_2-h_ihole.theta_1)/num_pts);
    end
end
rotor_points = [points(end,:);points];
outer_rotor_ring = rotor_points(num_pts+1:end,:);

%Rotate points
h_rot.pts = rotor_points;
rotor_points = rotate_pts(h_rot);

be = gds_element('boundary', 'xy',rotor_points,'layer',h_ihole.SOI);
str_name = sprintf('Partial_arm_[%d,%d],[%d]',round(h_ihole.p0(1)),round(h_ihole.p0(2)),round(h_ihole.radii(1)));
c_joint = gds_structure(str_name,be);

% Add etch holes to partial annulus
h_etch.regions = cell(1,1);
h_etch.r = h_ihole.etch_r;
h_etch.rotation_angle = h_rot.theta;
h_etch.rotation_point = h_rot.p0;
h_etch.undercut = h_ihole.EH_UNDERCUT;
h_etch.layer = h_ihole.SOI_HOLE;
section.p0 = h_ihole.p0;
section.type = 'partial_annulus';
section.r = h_ihole.radii;
section.theta = [h_ihole.theta_2,h_ihole.theta_1]; 
h_etch.regions{1,1} = section;
C_holes = etch_hole(h_etch);

%% Add bottom right segment to half annulus
h_ihole.theta_1 = 3*pi/2;
h_ihole.theta_2 = 2*pi;

num_pts = round(h_ihole.radii(1)/2);               % How many points to be in arc
points = zeros((num_pts+1),2);    % Matrix to store annulus pts
mult = 1; 
j=1;
for i=1:num_pts+1
    points(i + (j-1)*(num_pts+1),1) = h_ihole.p0(1) + (h_ihole.radii(j)+.05)*cos(h_ihole.theta_1+mult*(i-1)*abs(h_ihole.theta_2-h_ihole.theta_1)/num_pts);
    points(i + (j-1)*(num_pts+1),2) = h_ihole.p0(2) + (h_ihole.radii(j)+.05)*sin(h_ihole.theta_1+mult*(i-1)*abs(h_ihole.theta_2-h_ihole.theta_1)/num_pts);
end
rotor_points = points;
%Add other points to fill out rectagular section
rotor_points = [rotor_points;rotor_points(end,:)+[h_ihole.radii(2)-h_ihole.radii(1) 0];];
rotor_points = [rotor_points;rotor_points(end,:)+[0 -h_ihole.radii(1)];];
remember_pt = rotor_points(end,:);

%Rotate points
h_rot.pts = rotor_points;
rotor_points = rotate_pts(h_rot);


be = gds_element('boundary', 'xy',rotor_points,'layer',h_ihole.SOI);
str_name = sprintf('BR_curve_[%d,%d],[%d]',round(h_ihole.p0(1)),round(h_ihole.p0(2)),round(h_ihole.radii(1)));
br_curved = gds_structure(str_name,be);

%Add etch holes to that segment
h_etch.regions = cell(1,1);
section2.p0 = remember_pt;  %Bottom left point of rectangle
section2.type = 'lcurve';
section2.w = h_ihole.radii(2);
section2.l = h_ihole.radii(1);
section2.rcurve = points;

h_etch.regions{1,1} = section2;
rc_etch_holes = etch_hole(h_etch);

%% Add bottom left segment to half annulus
h_ihole.theta_1 = pi;
h_ihole.theta_2 = 3*pi/2;

num_pts = round(h_ihole.radii(1)/2);               % How many points to be in arc
points = zeros((num_pts+1),2);    % Matrix to store annulus pts
mult = 1; 
j=1;
for i=1:num_pts+1
    points(i + (j-1)*(num_pts+1),1) = h_ihole.p0(1) + (h_ihole.radii(j)+.05)*cos(h_ihole.theta_1+mult*(i-1)*abs(h_ihole.theta_2-h_ihole.theta_1)/num_pts);
    points(i + (j-1)*(num_pts+1),2) = h_ihole.p0(2) + (h_ihole.radii(j)+.05)*sin(h_ihole.theta_1+mult*(i-1)*abs(h_ihole.theta_2-h_ihole.theta_1)/num_pts);
end
rotor_points = [points];
%Add other points to fill out rectagular section
rotor_points = [rotor_points;rotor_points(end,:)+[-h_ihole.radii(2) 0];];
rotor_points = [rotor_points;rotor_points(end,:)+[0 h_ihole.radii(1)];];
remember_pt2 = rotor_points(end-1,:);

%Rotate points
h_rot.pts = rotor_points;
rotor_points = rotate_pts(h_rot);

be = gds_element('boundary', 'xy',rotor_points,'layer',h_ihole.SOI);
str_name = sprintf('BL_curve_[%d,%d],[%d]',round(h_ihole.p0(1)),round(h_ihole.p0(2)),round(h_ihole.radii(1)));
bl_curved = gds_structure(str_name,be);

%Add etch holes to that segment
h_etch.regions = cell(1,1);
section2.p0 = remember_pt2;  %Bottom left point of rectangle
section2.type = 'rcurve';
section2.w = h_ihole.radii(2);
section2.l = h_ihole.radii(1);
section2.rcurve = points;

h_etch.regions{1,1} = section2;
lc_etch_holes = etch_hole(h_etch);
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