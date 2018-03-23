function out = midpt(p1,p2,d)
%Function that finds a point on the line between two points
% p1 = point 1
% p2 = point 2
% d = how close to p2 we are. 0 gives midpoint at p1, 1 gives midpoint at
% p2. d=.5 gives the true midpoint

out = [(1-d)*p1(1)+d*p2(1),(1-d)*p1(2)+d*p2(2)];