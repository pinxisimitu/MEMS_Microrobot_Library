function [bigX bigY EH_pts] = jdrawrack(P,N,PA,x,y,phi,position,n_teeth)
% 
% This MATLAB function draws a spur gear with involute teeth.
% P = diametral pitch
% N = number of teeth
% PA = pressure angle in radians
% x, y are linear offset of center
% phi is rotational offset in radians
% postition is +/-1. +1 the teeth point right, -1 the teeth point left

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  
%  Copyright 2009, 2010 Matt Moses   
%  mmoses152 at earthlink dot net
% 
%  This file is part of Multigear. 
%
%  Multigear is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  Multigear is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
%
%  You should have received a copy of the GNU General Public License
%  along with Multigear.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following websites were helpful in providing definitions of 
% terminology and formulas for gear geometry:
%
% The Involute Curve, Drafting a Gear in CAD and Applications
% by Nick Carter 
% http://www.cartertools.com/involute.html
%
% Introduction to Mechanisms
% Chapter 7 Gears
% by Yi Zhang, with Susan Finger and Stephannie Behrens 
% http://www.cs.cmu.edu/~rapidproto/mechanisms/chpt7.html
%
% Dynamic Mechanisms Tutorial 
% http://www.ul.ie/~nolk/gears.htm#GEARS

D = N/P; % pitch diameter
DB = D * cos(PA); % base circle diameter
a = 1/P; % addendum (also known as "module")
d = 1.25/P; % dedendum (this is set by a standard)
DO = D + 2*a; % outside diameter
DR = D - 2*d; % root diameter

% angular thickness of gear where profile intersects pitch circle
theta_mod = pi/N;

% angular coordinate at intersection of profile and pitch circle
theta_p = ccinv(D/2, DB/2);

% The involute profile is approximated by connected line segments
% nseg is the number of line segments in the involute profile 
% of one side of the tooth
nseg = 10;

% list of radius values at which we will calculate angular coordinates
r_pro1 = linspace(DB/2,DO/2,nseg);

% now clip r_pro1 in case the base circle is smaller than the root circle
r_temp = (DR/2)*(r_pro1 < DR/2) + r_pro1 .* (r_pro1 >= DR/2);
r_pro1 = r_temp;

% angular coordinates corresponding to radius values
t_pro1 = ccinv(r_pro1, DB/2);

% Rtooth and Ttooth are polar coordinates for one tooth of the gear.
% 'R' as in Radius, 'T' as in Theta
Rtooth = [DR/2, r_pro1, r_pro1(end:-1:1),DR/2];
Ttooth = [t_pro1(1), t_pro1, theta_mod + 2*theta_p - t_pro1(end:-1:1), ...
    theta_mod + 2*theta_p - t_pro1(1)];

%Find etch hole points for a single tooth 
%Convert to XY
X = Rtooth.*cos(Ttooth); 
Y = Rtooth.*sin(Ttooth);

etch_dist = 6;
EH = 4;

%First get height of tooth
height = max(X) - min(X);
base_width = max(Y) - min(Y);
theta = 16.9*pi/180;

rows = round(height/(etch_dist+EH));

vspacing = height/rows + EH/2;
x_start = min(X);
y_start = min(Y);

xpts = [];
ypts = [];
for i=1:rows
    %Find width of row
    curr_x = x_start + (i-1)*(vspacing);
    curr_w = base_width - 2*(vspacing*(i-1))*tan(theta); 
    num_eh = round((curr_w - etch_dist)/(etch_dist + EH));
    act_spacing = (curr_w - etch_dist)/(num_eh) - EH;
    curr_y = y_start + tan(theta)*vspacing*(i-1);
    xpts = [xpts curr_x*ones(1,num_eh)];
    w_step = curr_w/num_eh;
    ypts = [ypts curr_y - w_step/2 + cumsum(w_step*ones(1,num_eh))]; %/2
end

EH_pts_T = atan(ypts./xpts);
EH_pts_R = sqrt(xpts.^2 + ypts.^2);
bigT_eh = [];
bigR_eh = [];


% This is the angular thickness of the 'body' of one tooth,
% used to get phase alignment right.
angThick = Ttooth(end) - Ttooth(1);

% Initialize vectors for the whole gear
bigT = [];
bigR = [];

% This incrementally adds coordinates for the rest of the teeth to the
% vectors bigT and bigR

%Make one tooth
bigT = [bigT, Ttooth + pi/180*(90+90*position)];
bigR = [bigR, Rtooth];

bigT_eh = [bigT_eh, EH_pts_T + pi/180*(90+90*position)];
bigR_eh = [bigR_eh, EH_pts_R];


% convert polar coordinates to XY coordinates, and incorporate phase
% rotation and shifting of the center
bigX = [];
bigY = [];
bigX_eh = [];
bigY_eh = [];

for k=0:n_teeth-1
    bigX = [bigX [bigR] .* cos([bigT] - angThick/2 + phi) + x];
    bigY = [bigY [bigR] .* sin([bigT] - angThick/2 + phi) + y + k*pi/P];
    
    bigX_eh = [bigX_eh [bigR_eh] .* cos([bigT_eh] - angThick/2 + phi) + x];
    bigY_eh = [bigY_eh [bigR_eh] .* sin([bigT_eh] - angThick/2 + phi) + y + k*pi/P];

end

%Shift rack to where we want it
mbx = max(bigX);
mby = max(bigY);

bigX = bigX - mbx + x;
bigY = bigY - mby + y;

bigX_eh = bigX_eh - mbx + x;
bigY_eh = bigY_eh - mby + y;




EH_pts = [bigX_eh' bigY_eh'];

