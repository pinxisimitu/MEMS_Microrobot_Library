% This script creates an example layour file (.GDS) containing some of the
% basic structures from the MEMS layout library. 

% 03/22/2018
% Joey Greenspun

clc
close all
clear 

tic                             % Timing Command, starts timer to see how long the sript runs for
gdsii_units(1e-6,1e-9)          % Sets units for (user, database) of the GDS file

% Define the GDS layer names 
default_layer_properties;

% Define text label properties
h_label.layer = SOI;
h_label.height = 10;             

% Etching/Release parameters 
default_etch_properties;

% Don't do any VREP simulation
global VREP_ignore;
VREP_ignore = 1;

%% Create a rectangle

h_rect.x = 0;                   % X coordinate of bottom left point
h_rect.y = 0;                   % Y coordinate of bottom left point
h_rect.w = 20;                  % Dimension in x direction
h_rect.l = 100;                 % Dimension in y direction
h_rect.layer = SOI;             % Set layer of rectangle
rect_ex = rect(h_rect);         % Function to create the rectangle GDS structure
                                % NOTE: Make sure GDS structures have unique names
                                
h_rect.x = 70;                  % Shift X coordinate to the right
h_rect.w = 50;                  % Change width of rectangle
h_rect.rounded = 5;             % Radius of curvature of edges (default is 0)
rect_rounded = rect(h_rect);    % Function to create a second rectangle GDS structure


h_rect.y = 120;                 % Shift X coordinate to the right
h_rect.w = 30;                  % Change width of rectangle
h_rect.etch = 1;
h_rect.ETCH_R = ETCH_R;
h_rect.etch_layer = SOIHOLE;
h_rect.UNDERCUT = UNDERCUT;
h_rect.rounded = 5;             % Radius of curvature of edges (default is 0)
rect_rounded3 = rect(h_rect);   % Function to create a second rectangle GDS structure

h_rect.y = 250;                 % Shift X coordinate to the right
h_rect.w = 30;                  % Change width of rectangle
h_rect.p0 = [h_rect.x h_rect.y];
h_rect.theta = 45*pi/180;
h_rect.rounded = 5;             % Radius of curvature of edges (default is 0)
rect_rounded4 = rect(h_rect);   % Function to create a second rectangle GDS structure
h_rect.theta = 0;
h_rect.rounded = 0;


h_label.text = 'rectangle';     
h_label.p0 = [0,-20];           
rect_lab = add_label(h_label);  % Add label underneath rectangle
h_rect.etch = 0;


%% Create an arbitrary path 

p0 = [150 0];
p1 = [160 0];
p2 = [160 10];
p3 = [170 10];
p4 = [170 0];
p5 = [190 10];                          % Define points in the path [x,y]
h_path.pts = [p0;p1;p2;p3;p4;p5];       % Make sure each row is one point [x,y]

h_path.w = 3;                           % Width of path
h_path.layer = SOI;                     % Set layer of path
path_ex = m_path(h_path);               % Function to create the path GDS structure

h_path.pts = h_path.pts + [0 30];       % Shift points up
h_path.smooth = 1;                      % Set flag to smooth out path [default is 0]
smoothed_path_ex = m_path(h_path);      % Create a second smoothed path

h_label.text = 'path';
h_label.p0 = [150,-20];
path_lab = add_label(h_label);          % Add label undneath path

%% Create a circle
h_cir.r = 30;                           % Circle radius
h_cir.layer = SOI;                      % Circle layer
h_cir.x = 330;                          % x position of circle center
h_cir.y = 30;                           % y position of circe center
circle_ex = circle(h_cir);              % Create GDS structure

h_label.text = 'circle';
h_label.p0 = [310,-20];
circ_lab = add_label(h_label);          % Add label undneath path

% Put etch holes in that circle           
h_etch.r = 2;                           % Set radius of etch holes
h_etch.undercut = 6;                    % Set oxide undercut (this is an estimate)
section.p0 = [h_cir.x, h_cir.y];        % Define center of etch hole array
section.type = 'annulus';               % Set etch flag to annulus
section.r = [0 h_cir.r];                % Set inner and outer radius of annulus to fill with etch holes
h_etch.regions = cell(1,1);             % Allocate a cell array to hold a region to etched  
h_etch.regions{1,1} = section;          % Assign the structure to the h_etch.regions field 
circle_holes = etch_hole(h_etch);       % Generate GDS structure containing all etch holes

% Create an array or circles
side_len = 10;                          % Number of rows/columns
x_pts = zeros(side_len,side_len);       % Preallocate x points vector
y_pts = zeros(side_len,side_len);       % Preallocate y points vector
gap = 10;                               % Pitch of circles

for i = 1:10                            % For loop to find centers of circles
    for j = 1:10
        x_pts(i,j) = (i-1)*gap;
        y_pts(i,j) = (j-1)*gap;
    end
end

xy_pts = [reshape(x_pts,numel(x_pts),1) reshape(y_pts,numel(y_pts),1)]; % This line reshapes the matricies into vectors and appends Y to X
h_cir.r = 2;                            % Set circle radius to 2
xy_pts = xy_pts + [280 70];             % Shift array above other circle
h_cir.xy = xy_pts;                      % h_cir.xy holds center points of all circles to be drawn
circle_array_ex = circle(h_cir);        % Create circle array GDS structure 

%% Create an arbitrary shape

p0 = [500 0];
p1 = [550 0];
p2 = [550 50];
p3 = [600 100];
p4 = [600 150];
p5 = [550 200];
p6 = [500 200]; 
p7 = [450 150];
p8 = [450 100];
p9 = [500 50];

% Define points in the path [x,y]
h_shape.pts = [p0;p1;p2;p3;p4;p5;p6;p7;p8;p9];       % Make sure each row is one point [x,y]
h_shape.layer = SOI;                      % Set layer of shape
shape_ex = m_shape(h_shape);              % Function to create the shape GDS structure

h_label.text = 'shape';
h_label.p0 = [500,-20];
shape_lab = add_label(h_label);           % Add label undneath path

%% Create serpentine spring

h_ss.p1 = [700 0];                        % First anchor point of the SS
h_ss.p2 = [700 300];                      % Second anchor point of the SS
h_ss.n = 11;                              % Number of meanders 
h_ss.w = 4;                               % Width of each beam 
h_ss.dpp = 95;                            % Peak to peak distance of meanders (beam length)
h_ss.layer = SOI;                         % Put SS on SOI layer 
ss_ex = s_spring(h_ss);                   % Function to create SS GDS structure

h_label.text = 'serpentine spring';
h_label.p0 = [625,-20];
ss_lab = add_label(h_label);              % Add label undneath SS


%% Create a pin joint

h_joint.p0 = [1000 200];                  % Center of circle for pin joint
h_joint.r = [50 75];                      % [pin radius, outer radius]
h_joint.opening_theta = 100;              % Joint openeing in degrees

h_joint.inner_arm_theta = 0;              % Angle (in degrees) of the inner arm
h_joint.inner_arm_w = 50;                 % Inner Arm Width
h_joint.inner_arm_l = 200;                % Inner Arm Length

h_joint.outer_arm_w = 50;                 % Outer arm width
h_joint.outer_arm_l = 200;                % Outer arm length
h_joint.outer_arm_theta = 90;             % Outer arm angle

h_joint.layer = SOI;
joint_ex = joint(h_joint);

h_label.text = 'pin joint';
h_label.p0 = [1000,-20];
joint_lab = add_label(h_label);           % Add label undneath SS

%% Make triangle with etch holes in it

tri_w = 500;
tri_height = 500;

p0 = [1400 0];
p1 = [1400+tri_w 0];
p2 = [1400 tri_height];

% Define points in the path [x,y]
h_shape.pts = [p0;p1;p2];                % Make sure each row is one point [x,y]
h_shape.layer = SOI;                     % Set layer of shape
triangle = m_shape(h_shape);             % Function to create the shape GDS structure

h_label.text = 'triangle etch holes';
h_label.p0 = [1400,-20];
t = add_label(h_label);                  % Add label undneath path

% %Add etch holes to the triangle
h_etch.regions = cell(1,1);              % Allocate regions 
section2.p0 = p0;                        % Bottom left point of rectangle
section2.type = 'rcurve';                % Set etch flag to right curve
section2.w = tri_w;                      % Set width of etch area
section2.l = tri_height;                 % Set height of etch area
section2.rcurve = [linspace(p0(1),p1(1),100)' linspace(p2(2),p0(2),100)'];  % Gernerate curve data from poins
 
h_etch.regions{1,1} = section2;
tri_etchholes = etch_hole(h_etch);

%% Make an array of structures using different methods

% Array rectangles using eval statements
num_rects = 10;
gap = 6;

h_rect.x = 2000;                % X coordinate of bottom left point
h_rect.y = 0;                   % Y coordinate of bottom left point
h_rect.w = 5;                   % Dimension in x direction
h_rect.l = 100;                 % Dimension in y direction
h_rect.layer = SOI;             % Set layer of rectangle
h_rect.ynum = 1;                % Reset number of copies to 1
h_rect.xnum = 1;                % Reset number of copies to 1

for i=1:num_rects
    h_rect.x = h_rect.x + h_rect.w + gap;
    eval_string = sprintf('rect_%d = rect(h_rect);',i); % make string that will create GDS structure
    eval(eval_string);                                  % Evaluate the string
end

% Array rectangles using cell arrays
rect_cell = cell(num_rects,1);         % Initialize cell array
gap = 6;                               % Gap between each rectangle 

h_rect.x = 2000;                       % X coordinate of bottom left point
h_rect.y = 200;                        % Y coordinate of bottom left point
h_rect.w = 5;                          % Dimension in x direction
h_rect.l = 100;                        % Dimension in y direction
h_rect.layer = SOI;                    % Set layer of rectangle

for i=1:num_rects
    h_rect.x = h_rect.x + h_rect.w + gap;
    rect_cell{i} = rect(h_rect);
end

h_label.text = 'rectangle arrays';
h_label.p0 = [2000,-20];
rect_array_lab = add_label(h_label);              % Add label undneath the motor


%% Inchworm motors

h_motor.pos = [3500 700];               % Shuttle position
h_motor.layer = [SOI, SOIHOLE];         % Layers for [SOI and SOIHOLE]
h_motor.N = 30;                         % Number of comb finger gaps per motor
h_motor.travel = 500;                   % Total travel of the shuttle
h_motor.angle = 0;                      % Angle of shuttle in degrees
[m1 moto_pts]= motor(h_motor);       % First output is the GDS structure, second contains points of interest for motor

h_motor.pos = [3000 2700];              % Shuttle position
h_motor.layer = [SOI, SOIHOLE];         % Layers for [SOI, SOIHOLE]
h_motor.N = 30;                         % Number of comb finger gaps per motor
h_motor.travel = 500;                   % Total travel of the shuttle
h_motor.angle = 45;                     % Angle of shuttle in degrees
h_motor.ground_serpentine = 1;          % Adds the serpentine srpings along GND contacts
[m2 moto_pts]= motor(h_motor);       % First output is the GDS structure, second contains points of interest for motor

h_motor.pos = [3700 5200];              % Shuttle position
h_motor.layer = [SOI, SOIHOLE];         % Layers for [SOI, SOIHOLE]
h_motor.N = 30;                         % Number of comb finger gaps per motor
h_motor.travel = 1500;                  % Total travel of the shuttle
h_motor.num_inch_sets = 2;              % Number of sets of GCAs to generate
h_motor.angle = 45;                     % Angle of shuttle in degrees
h_motor.ground_serpentine = 1;          % Adds the serpentine srpings along GND contacts
[m3 moto_pts]= motor(h_motor);       % First output is the GDS structure, second contains points of interest for motor


h_label.text = 'inchworm motor';
h_label.p0 = [2500,-20];
inch_lab = add_label(h_label);          % Add label undneath the motor

%% Collect all structures and save them to a GDS file

root = 'Tutorial';

% Nothing in this block should be changed except for the above "root" variable.
% This serves as the root name for saved GDS files. When you run this script
% consecutive times, it will append an incrementing number to this root.


% DO NOT CHANGE ANYTHING BELOW THIS LINE!!!
% This part of the code gathers all of the GDS structures and packages
% them into a .GDS file.

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
            fprintf('Empty Cell! Something went wrong in creating: %s!!\n',a(i).name)
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

% Create a top level structure
ts = gds_structure('TOP');
ts = add_ref(ts, b);  % add sref elements to top level

% Create a library to hold the structure
glib = gds_library('lib',b,ts);

% Finally write the library to a file, looking first for other files names
% the same thing and incrementing the appended digit if need be. 
ending = '.gds';
str = sprintf('%s*%s',root,ending);
a = dir([str]);
suffix = [];
if ~isempty(a)
    k = [];
    for ii=1:length(a)
        k = [k str2num(strtok(strtok(a(ii).name,'.'),root))];
    end
    if isempty(k)
        suffix = 1;
    else
        suffix = max(k)+1;
    end
end

write_gds_library(glib, [root num2str(suffix) ending]);

toc %Timing Command, stops the timer and prints how long the sript ran for