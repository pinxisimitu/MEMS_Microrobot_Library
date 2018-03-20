% trigearRack.m
%
% This is a MATLAB script that calculates center locations for two gears 
% that are in simultaneous mesh with a third gear and a straight rack.

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

clear 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These are the parameters that can be freely changed to get different
% geometries.
%
% DP = diametral pitch. A gear of 1" pitch diameter has this many teeth.
% There are two instances of gear1, each in contact with the rack.
% Gear 2 is the single top gear.
%
DP = 12;
n1 = 16; % number of teeth on gear 1 
n2 = 8; % number of teeth on gear 2
nRack = 20; % number of teeth on rack
shiftRack = 5; % number of teeth rack is shifted - this is for aesthetics
rholes = 0.1; % radius of gear bore holes in drawing - for aesthetics
%
% Several non-unique solutions may be found.  Or maybe NO solutions will be
% found.  You can change this paramter to select different candidate sets
% of coordinates.  Typically it should be 0 or 1, but can be negative or
% larger positive...
solutionIndexShift = 0; 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1 = n1/DP/2; % pitch radius of gear 1
r2 = n2/DP/2; % pitch radius of gear 2
gamma = 2 * pi / n2; % angular spacing of teeth on gear 2

% radius to clearance circles around outside diameter of each gear
rod1 = (n1/DP + 2/DP)*1.05; % gear 1 clearance 
rod2 = (n2/DP + 2/DP)*1.05; % gear 2 clearance

% theta1 is the space of angles that we will search.
% More points means more accuracy, since we are just searching via brute
% force instead of using a more sophisticated zero-finder.
theta1 = linspace(0,pi/2,10000);

% see accompanying figure for definition of these angles
theta2 = 2 * (r1+r2)/r1 * sin(theta1) + theta1;
alpha = -theta1 - r1/r2 * (theta1 + theta2);

% This angle has to be an integer multiple of gamma
angDiff = theta1 - alpha;

% This finds integer multiples in a naive and inefficient way
% Good thetas end up in the vector thetaGoods
n = 0;
for k = 0:floor(angDiff(end)/gamma)    
    for i = 1:length(theta1)-1
        if (((angDiff(i) - k*gamma) < 0) && ((angDiff(i+1) - k*gamma) > 0))
            n = n + 1;
            thetaGoods(n) = (theta1(i) + theta1(i+1)) / 2;
        end
    end
end

% % Uncomment this block to plot the angular constraint function
% figure;
% hold on;
% for k = 0:floor(angDiff(end)/gamma)
%     plot([0 theta1(end)],[k*gamma k*gamma], 'k');
% end
% plot(theta1, angDiff);
% xlabel('\theta_1 in radians');
% ylabel('(\theta_1 - \alpha) and (k \gamma) in radians');

% vector of center-to-center distances corresponding to theta1
c2c = 2*(r1+r2)*sin(theta1);

% interpolate good center-to-center distances
c2cGoods = interp1(theta1,c2c,thetaGoods);

% Now scan the vector of good c2c distances to find the first c2c distance
% that actually allows the gears to clear each other
pind = 0;
for (i = 1:length(c2cGoods))
    if (c2cGoods(i)>rod1)
        pind = i;
        break
    end
end

% increment the solution index to plot solutions with larger or smaller
% clearance
pind = pind + solutionIndexShift;

% Here is center to center distance of the two gears in conctact with the
% rack
display(sprintf ...
    ('center-center distance of gears in contact with rack: %f', ...
    c2cGoods(pind)));

% And here is the location of the top gear relative to the first lower gear
xdis = (r1+r2) * sin(thetaGoods(pind));
ydis = (r1+r2) * cos(thetaGoods(pind));
display(sprintf('x distance of top gear: %f', xdis));
display(sprintf('y distance of top gear: %f', ydis));

% And finally this part draws a picture so we can see what is going on and 
% maybe print out some patterns.
figure;

% % Uncomment this block to draw the pitch diameters on the gears
% circle(r1*thetaGoods(pind),r1,r1);
% hold on
% circle(2*(r1+r2)*sin(thetaGoods(pind)) + r1*thetaGoods(pind),r1,r1);
% circle(r1*thetaGoods(pind)+xdis,r1+ydis,r2);

% Angular offset of gear 3
gear3phi = pi/2 -(1/r1) * ...
    (2*(r1+r2)*sin(thetaGoods(pind)) + r1*thetaGoods(pind));

% draw gears
gear1 = drawInvolute(DP, n1, 20*pi/180, ...
    r1*thetaGoods(pind), r1, -thetaGoods(pind)+pi/2,[0 360]);
hold on
gear2 = drawInvolute(DP, n2, 20*pi/180, ...
    r1*thetaGoods(pind) + xdis, r1+ydis, -thetaGoods(pind) - pi/2 + pi/n2,[0 360]);
gear3 = drawInvolute(DP, n1, 20*pi/180, ...
    2*(r1+r2)*sin(thetaGoods(pind)) + r1*thetaGoods(pind), r1, gear3phi,[0 360]);

% draw center holes
hole1 = circle(r1*thetaGoods(pind),r1,rholes);
hole2 = circle(2*(r1+r2)* ...
    sin(thetaGoods(pind)) + r1*thetaGoods(pind),r1,rholes);
hole3 = circle(r1*thetaGoods(pind)+xdis,r1+ydis,rholes);

% draw rack
rack1 = drawRack(DP, nRack, 20*pi/180, ...
    mod(pi*r1, pi/DP) - shiftRack * pi/DP, 0);

% set axes equal so things don't look funny
axis equal