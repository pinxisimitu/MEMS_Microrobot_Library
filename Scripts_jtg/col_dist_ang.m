function out = col_dist_ang(points)
%This function returns the distance between two points and the angle between them in radians
%points = [p1_x p1_y
%          p2_x p2_y];

x1 = points(1,1);
x2 = points(2,1);
y1 = points(1,2);
y2 = points(2,2);

dx = x2-x1;
dy = y2-y1;

out = [sqrt(dx^2+dy^2) atan2(dy,dx)];


