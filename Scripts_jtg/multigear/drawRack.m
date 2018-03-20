function handle = drawRack(P,N,PA,x,y)
% handle = drawRack(P, N, PA, x, y)
% 
% This MATLAB function draws a straight rack on the current plot.
% P = diametral pitch
% N = number of teeth
% PA = pressure angle in radians
% x, y = linear offset distances

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

a = 1/P; % addendum (also known as "module")
d = 1.25/P; % dedendum (this is set by a standard)

% find the tangent of pressure angle once
tanPA = tan(PA);

% length of bottom and top horizontal segments of tooth profile
botL = pi/P/2 - 2*d*tanPA;
topL = pi/P/2 - 2*a*tanPA;

% This is a vector of X coordinates for one tooth
toothX = x + [0, botL/2, botL/2 + (d+a)*tanPA, botL/2 + (d+a)*tanPA + ...
    topL, botL/2 + 2*(d+a)*tanPA + topL, botL + 2*(d+a)*tanPA + topL];

% This is a vector of Y coordinates for one tooth
toothY = y + [-d, -d, a, a, -d, -d]; 

% initialize big vector for the many-toothed rack
bigX = [toothX];
bigY = [toothY];

% add more teeth to the whole rack
for i = 1:N-1
    bigX = [bigX, toothX + i*pi/P];
    bigY = [bigY, toothY];
end

% add the rectangular envelope
bigX = [bigX, bigX(end), bigX(1), bigX(1)];
bigY = [bigY, bigY(end) - 4*a, bigY(end) - 4*a, bigY(1)];

% Draw the coordinates on the active plot
handle = line(bigX, bigY);