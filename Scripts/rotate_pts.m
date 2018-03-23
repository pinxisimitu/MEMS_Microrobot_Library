function out=rotate_pts(h_rot)
    %pts = column vector with two columns [x y;x2 y2;...]
    %out = column vector of rotated pts
    %theta = angle to rotate about in radians
    %p0 = point about which to rotate all other points
    theta = h_rot.theta; 
    out = zeros(size(h_rot.pts));
    for j=1:size(out,1)
        out(j,1) = cos(h_rot.theta) * (h_rot.pts(j,1) - h_rot.p0(1)) - sin(theta) * (h_rot.pts(j,2) - h_rot.p0(2)) + h_rot.p0(1);
        out(j,2) = sin(h_rot.theta) * (h_rot.pts(j,1) - h_rot.p0(1)) + cos(theta) * (h_rot.pts(j,2) - h_rot.p0(2)) + h_rot.p0(2); 
    end
end