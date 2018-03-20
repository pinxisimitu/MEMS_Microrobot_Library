%Testing drawInvolte script
close all

r1 = 1000;   %Radius
DP = .02;
angle = [270 360];
phi = 20*pi/180;

n1 = 2*r1*DP; % number of teeth on gear 1 

gear1 = drawInvolute(DP, n1, phi,0,0,0,angle);
d = 1.25/DP; % dedendum (this is set by a standard)

r1 = jdrawrack(DP, n1, 20*pi/180,r1+d,-340,0,1,10);
%r1 = jdrawrack(DP, n1, 20*pi/180,r1+d,130,0,1,10);
