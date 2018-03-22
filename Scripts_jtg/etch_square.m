function out = etch_square(h_square)
% Function to make a square for etch holes
% x = x coordinate of circle
% y = y coordinate of circle
% r  = radius of circle
% n = number of points
% layer = GDS layer to place circle into

if ~isfield(h_square,'compact')
    h_square.compact = 0;
end

if isfield(h_square,'xy')
    h_square.compact = 1;
end

if ~isfield(h_square,'rp')
    h_square.rp = 0;
end


if ~isfield(h_square,'rotation_theta')
    h_square.rotation_theta = 0;
end

if ~isfield(h_square,'rotation_center')
    h_square.rotation_center = [0 0];
end


if ~isfield(h_square,'n')
    h_square.n = round(pi*2*h_square.r); % puts a point every um of circumference
end

if  h_square.compact==0
    points = zeros(h_square.n,2);
    
    for i=1:h_square.n
        points(i,1) = h_square.x + h_square.r*cos(i*2*pi/h_square.n);
        points(i,2) = h_square.y + h_square.r*sin(i*2*pi/h_square.n);
    end
    
    be = gds_element('boundary', 'xy',points,'layer',h_square.layer);
    str_name = sprintf('C_[%d,%d],%d%d',round(h_square.x),round(h_square.y),round(h_square.r),h_square.layer);
    out = gds_structure(str_name,be);
else
    be = cell(1,length(h_square.xy));
    for j=1:length(h_square.xy)
        h_square.x = h_square.xy(j,1);
        h_square.y = h_square.xy(j,2);
        
        points = zeros(4,2);
        points(1,:) = [h_square.x-h_square.r h_square.y-h_square.r];
        points(2,:) = [h_square.x+h_square.r h_square.y-h_square.r];
        points(3,:) = [h_square.x+h_square.r h_square.y+h_square.r];
        points(4,:) = [h_square.x-h_square.r h_square.y+h_square.r];
        
        % Rotate the square etch holes
        rot.pts = points;
        rot.theta = h_square.rotation_theta;
        rot.p0 = [h_square.x h_square.y];
        points =rotate_pts(rot);
        
        be{j} = gds_element('boundary', 'xy',points,'layer',h_square.layer);
    end
    str_name = sprintf('Sq_[%d,%d],%d%d',round(h_square.xy(1,1)),round(h_square.xy(1,2)),round(h_square.r),h_square.layer);
    out = gds_structure(str_name,be);
end

if h_square.rp==1
    out = points;
end

