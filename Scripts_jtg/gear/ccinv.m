function theta = ccinv(r,rb)
% theta = ccinv(r, rb)
% 
% This MATLAB function returns the angular coordinates for points on 
% the "clockwise" involute of a circle, 
% given radius of the base circle rb, 
% and radial coordinates of the point r.

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

phi = (r.^2 ./ rb .^2 - 1).^0.5;
theta = phi - atan(phi);