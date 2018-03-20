function out = s_spring(h_ss)
% Function to make a serpentine spring between two points
% p1 = first anchor of spring  [x1,y1]
% p2 = second anchor of spring [x2,y2]
% n  = number of periods in spring
% w  = width of spring
% dpp = peak to peak distance of springs (perpendicular to line p1-p2)
% layer = GDS layer to place spring into

if ~isfield(h_ss,'rp')
    h_ss.rp = 0;
end


% Flags the set the endpoints of the serpentine springs
% If set to 1, the serpentine spring will terminate at the center of the
% meaneder, otherwise it will terminate at the bottom of the meander. 

if ~isfield(h_ss,'endpts') 
    h_ss.endpts = [.5 .5];
end

p1 = h_ss.p1;
p2 = h_ss.p2;
n = h_ss.n;
w = h_ss.w;
dpp = h_ss.dpp;
layer = h_ss.layer;

%Find the number of linear sections
sec = 2*(n+1);

xstep = (p2(1)-p1(1))/sec;
ystep = (p2(2)-p1(2))/sec;
theta = atan(ystep/xstep);

% Initialize first three and last three points
points = zeros(2*sec,2);
points(1,:)= p1;
points(2,:)= [p1(1)+xstep,p1(2)+ystep];
points(3,:)= [p1(1)+xstep-dpp/2*sin(theta),p1(2)+ystep+dpp/2*cos(theta)];

points(2*sec-2,:)=  [p2(1)-xstep+dpp/2*sin(theta),p2(2)-ystep-dpp/2*cos(theta)];
points(2*sec-1,:)= [p2(1)-xstep,p2(2)-ystep];
points(2*sec,:)= p2;

%Loop to create points on serpentine spring
for i = 4:2*sec-3
    if mod(i,2)==0
        points(i,:) = [points(i-1,1)+xstep,points(i-1,2)+ystep];
    else
        if mod(i,4)==1
           points(i,:) = [points(i-1,1)+dpp*sin(theta),points(i-1,2)-dpp*cos(theta)];
        else
           points(i,:) = [points(i-1,1)-dpp*sin(theta),points(i-1,2)+dpp*cos(theta)];
        end
    end
end

% Modify first two and last two points if need be
for i=3:length(points)
    points(i,:) = points(i,:) - (h_ss.endpts(1)-0.5)*[dpp*sin(theta),-dpp*cos(theta)];
end

for i=length(points)-1:length(points)
    points(i,:) = points(i,:) + (h_ss.endpts(2)-0.5)*[dpp*sin(theta),-dpp*cos(theta)];
end


be = gds_element('path', 'xy',points,'width', w,'layer',layer);
str_name = sprintf('SS_[%d,%d],[%d,%d]',round(p1(1)),round(p1(2)),round(p2(1)),round(p2(2)));

if(h_ss.rp == 1)
    out = points;
else
    out = gds_structure(str_name,be);
end