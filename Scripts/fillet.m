function out = fillet(h_fillet)
% Function to create fillets on corners
% p0 = center point of fillet
% p1 = left edge of fillet
% p2 = right edge of fillet (optional)
% theta = angle between p0->p1 and p0->p2
% d2 = length of p0->p2 
% d = factor between ~0.2 and 1 that determines how rounded the fillet is.
% layer = GDS layer to place circle into
% downsample = the number by which to downsample the fillets

if isfield(h_fillet,'p2')
    d1 = sqrt((h_fillet.p0(1) - h_fillet.p1(1))^2 + (h_fillet.p0(2) - h_fillet.p1(2))^2);
    d2 = sqrt((h_fillet.p0(1) - h_fillet.p2(1))^2 + (h_fillet.p0(2) - h_fillet.p2(2))^2);
    d3 = sqrt((h_fillet.p2(1) - h_fillet.p1(1))^2 + (h_fillet.p2(2) - h_fillet.p1(2))^2);
    h_fillet.theta = acos((d1^2+d2^2-d3^3)/(2*d1*d2));
    h_fillet.d2 = d2; 
end

if ~isfield(h_fillet,'p3mp')
    h_fillet.p3mp = .5;
end

if ~isfield(h_fillet,'rp')
    h_fillet.rp = 0;
end

if ~isfield(h_fillet,'downsample')
    h_fillet.downsample = 10;
end


%Find central point for spline fit

p3 = midpt(h_fillet.p1,h_fillet.p2,h_fillet.p3mp);
p4 = midpt(h_fillet.p0,p3,h_fillet.d);
init_points = [h_fillet.p1(1) p4(1) h_fillet.p2(1);h_fillet.p1(2) p4(2) h_fillet.p2(2)];

interm_pts = fnplt(cscvn(init_points));

%Downsample the number of points.
final_pts = [interm_pts(1,1) downsample(interm_pts(1,2:end-1),h_fillet.downsample) interm_pts(1,end);
             interm_pts(2,1) downsample(interm_pts(2,2:end-1),h_fillet.downsample) interm_pts(2,end)];

points = [final_pts';h_fillet.p0];
ca = gds_element('boundary', 'xy',points,'layer',h_fillet.layer);


str_name = sprintf('Fillet_[%d,%d]',round(h_fillet.p0(1)),round(h_fillet.p0(2)));
if h_fillet.rp ==1
    out = final_pts;
else
    out = gds_structure(str_name,ca);
end
