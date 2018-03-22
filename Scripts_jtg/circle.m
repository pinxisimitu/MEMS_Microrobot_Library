function out = circle(h_cir)
% Function to make a circle
% x = x coordinate of circle
% y = y coordinate of circle
% r  = radius of circle
% n = number of points
% layer = GDS layer to place circle into

if ~isfield(h_cir,'compact')
    h_cir.compact = 0;
end

if isfield(h_cir,'xy')
    h_cir.compact = 1;
end

if ~isfield(h_cir,'rp')
    h_cir.rp = 0;
end

if ~isfield(h_cir,'n')
    h_cir.n = round(pi*2*h_cir.r); % puts a point every um of circumference
end

if  h_cir.compact==0
    points = zeros(h_cir.n,2);
    
    for i=1:h_cir.n
        points(i,1) = h_cir.x + h_cir.r*cos(i*2*pi/h_cir.n);
        points(i,2) = h_cir.y + h_cir.r*sin(i*2*pi/h_cir.n);
    end
    
    be = gds_element('boundary', 'xy',points,'layer',h_cir.layer);
    str_name = sprintf('C_[%d,%d],%d%d',round(h_cir.x),round(h_cir.y),round(h_cir.r),h_cir.layer);
    out = gds_structure(str_name,be);
else
    be = cell(1,length(h_cir.xy));
    for j=1:length(h_cir.xy)
        h_cir.x = h_cir.xy(j,1);
        h_cir.y = h_cir.xy(j,2);
        
        points = zeros(h_cir.n,2);
        for i=1:h_cir.n
            points(i,1) = h_cir.x + h_cir.r*cos(i*2*pi/h_cir.n);
            points(i,2) = h_cir.y + h_cir.r*sin(i*2*pi/h_cir.n);
        end
        
        be{j} = gds_element('boundary', 'xy',points,'layer',h_cir.layer);
    end
    str_name = sprintf('C_[%d,%d],%d%d',round(h_cir.xy(1,1)),round(h_cir.xy(1,2)),round(h_cir.r),h_cir.layer);
    out = gds_structure(str_name,be);
end

if h_cir.rp==1
    out = points; 
end

