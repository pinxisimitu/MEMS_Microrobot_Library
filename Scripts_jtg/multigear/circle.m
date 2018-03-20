function handle = circle(x,y,r)
% handle = circle(x,y,r)
%
% This MATLAB function draws a circle of radius r with center at x,y.

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

n = 50; % number of line segments the circle will be made up of
t = linspace(0,2*pi,n); % angular coordinates
xx = x + r*cos(t); % x coordinates
yy = y + r*sin(t); % y coordinates
handle = line(xx,yy); % draw it