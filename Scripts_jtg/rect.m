function out = rect(h_rect)
% Function to create and possibly array a rectangle
% x = x coordinate bottom left
% y = y coordinate bottom left
% w = width (in x direction)
% l = length (in y direction)
% layer = GDS layer to place circle into
% xnum = number of repeats in x direction
% xspace = gap between x repeats
% ynum = number of repeats in y direction
% yspace = gap between y repeats
% theta = angle by which to rotate rectangle (in radians)
% p0 = point about which to rotate
% rp = return points, if rp==1, returns column vector of verticies
% rounded = if non zero, set radius of curvature of corner to rounded 

%% Default Values
global VREP_ignore

if ~isfield(h_rect,'xnum')
    h_rect.xnum = 1;
end

if ~isfield(h_rect,'ynum')
    h_rect.ynum = 1;
end

if ~isfield(h_rect,'xspace')
    h_rect.xspace = 0;
end

if ~isfield(h_rect,'layer')
    h_rect.layer = 6;
end

if ~isfield(h_rect,'yspace')
    h_rect.yspace = 0;
end

if ~isfield(h_rect,'theta')             % Rotation angle
    h_rect.theta = 0;
end

if ~isfield(h_rect,'p0')                % Bottom left coordinate of rectangle
    h_rect.p0 = [0 0];
end

if ~isfield(h_rect,'rp')                % Return points instead of GDS structure
    h_rect.rp = 0;
end

if ~isfield(h_rect,'rounded')           % Radius of curvature on points of rectangles        
    h_rect.rounded = 0;
end

if ~isfield(h_rect,'etch')              % If set to 1, will add etch holes to rectangle
    h_rect.etch = 0;
end

if ~isfield(h_rect,'etch_layer')        % Layer to put etch holes in   
    default_layer_properties;
    h_rect.etch_layer = SOIHOLE;
end

if ~isfield(h_rect,'circle_etch')       % Circular etch holes or square?
    default_etch_properties;
    h_rect.circle_etch = CIRCULAR_ETCH;
end

if ~isfield(h_rect,'UNDERCUT')          % Rough Undercut for release
    h_rect.UNDERCUT = 6;
end


if ~isfield(h_rect,'ETCH_R')            % Radius of etch holes
     h_rect.ETCH_R = 2;
end

if ~isfield(h_rect,'VREP_group')
    h_rect.VREP_group = 'dummy_group';
end

%% 

points = [h_rect.x,h_rect.y;h_rect.x+h_rect.w,h_rect.y;h_rect.x+h_rect.w,h_rect.y+h_rect.l;h_rect.x,h_rect.y+h_rect.l;];
pts_orig = points;
if h_rect.w<0
    if h_rect.rounded>0
        p1 = points(1,:) + [-h_rect.rounded h_rect.rounded];
        p2 = points(2,:) + [h_rect.rounded h_rect.rounded];
        p3 = points(3,:) + [h_rect.rounded -h_rect.rounded];
        p4 = points(4,:) + [-h_rect.rounded -h_rect.rounded];
        
        pts = [p1;p2;p3;p4];
        
        %Draw quarter circles from each new point
        step = round(2*h_rect.rounded);
        angle = 0;
        points2 = [];
        for i=1:4
            theta = linspace(angle,angle-pi/2,step);
            points2 = [points2;[h_rect.rounded*cos(theta) + pts(i,1)]', [h_rect.rounded*sin(theta) + pts(i,2)]'];
            angle = angle - pi/2;
        end
        points = points2;
    end
else
    if h_rect.rounded>0
        p1 = points(1,:) + [h_rect.rounded h_rect.rounded];
        p2 = points(2,:) + [-h_rect.rounded h_rect.rounded];
        p3 = points(3,:) + [-h_rect.rounded -h_rect.rounded];
        p4 = points(4,:) + [h_rect.rounded -h_rect.rounded];
        
        pts = [p1;p2;p3;p4];
        %Draw quarter circles from each new point
        step = round(2*h_rect.rounded);
        angle = pi;
        points2 = [];
        for i=1:4
            theta = linspace(angle,angle+pi/2,step);
            points2 = [points2;[h_rect.rounded*cos(theta) + pts(i,1)]', [h_rect.rounded*sin(theta) + pts(i,2)]'];
            angle = angle + pi/2;
        end
        points = points2;
    end
end


ca = {};
count = 1;
for i = 1:h_rect.xnum
    for j = 1:h_rect.ynum
        points2 = bsxfun(@plus, points, [(i-1)*(h_rect.xspace+h_rect.w), (j-1)*(h_rect.yspace+h_rect.l)]);
        
        %Rotate the points here
        rot.pts = points2;
        rot.theta = h_rect.theta;
        rot.p0 = h_rect.p0;
        points2=rotate_pts(rot);

        ca{count} = gds_element('boundary', 'xy',points2,'layer',h_rect.layer);
        count = count + 1;
       
        if ~VREP_ignore
            % For VREP, find the center of the rectangle and place coboid there.
            % First sort points by Y value
            if h_rect.rounded>0
                points2 = bsxfun(@plus, pts_orig, [(i-1)*(h_rect.xspace+h_rect.w), (j-1)*(h_rect.yspace+h_rect.l)]);
                
                %Rotate the points here
                rot.pts = points2;
                rot.theta = h_rect.theta;
                rot.p0 = h_rect.p0;
                points2=rotate_pts(rot);
            end
            
            points3 = sortrows(points2,2);
            
            
            % Find width (distance between points3(1,:) and points3(2,:)
            VREP_w = sqrt((points3(1,2)-points3(2,2))^2+(points3(1,1)-points3(2,1))^2);
            
            % Find length (distance between points3(2,:) and either
            % points3(3,:) or points(4,:)...whichever is shorter
            VREP_la = sqrt((points3(3,2)-points3(2,2))^2+(points3(3,1)-points3(2,1))^2);
            VREP_lb = sqrt((points3(2,2)-points3(4,2))^2+(points3(2,1)-points3(4,1))^2);
            if VREP_la < VREP_lb
                VREP_l = VREP_la;
            else
                VREP_l = VREP_lb;
            end
            
            % Find the angle between points3(1,:) and points3(2,:)
            angle = atan2(points3(2,2)-points3(1,2),points3(2,1)-points3(1,1));
            
            % Find the center of the rectangle
            cent = midpt(points2(1,:),points2(3,:),.5);
            
            %generate name (cant use negative signs)
            if ~isfield(h_rect,'VREP_name')
                name = sprintf('R%d_%d%d',round(cent(1)),round(cent(2)),h_rect.layer(1));
                name = strrep(name,'-','n');        % Replace any negative signs with 'n'
                vgroup = h_rect.VREP_group;
            else
                name = h_rect.VREP_name;
                vgroup = 'none';
            end
            
            if h_rect.layer(1) ~= 8 && h_rect.rp == 0
                if ~isempty(strfind(h_rect.VREP_group,'_a')) || ~isempty(strfind(name,'_a'))
                    VREP_cuboid(cent,[VREP_w VREP_l],name,angle,1,vgroup);
                else
                    % This will put everything in contact with the floor, but
                    % need to do it to make joints work...
                    %VREP_cuboid(cent,[VREP_w VREP_l],name,angle,0,vgroup);
                    VREP_cuboid(cent,[VREP_w VREP_l],name,angle,1,vgroup);
                end
            end
        end
    end
end

nname = midpt(points(1,:),points(2,:),.8);
str_name = sprintf('R_[%d,%d]_%d',round(nname(1)),round(nname(2)),h_rect.layer);

ca = gds_structure(str_name,ca);

if h_rect.etch
    % Add etch holes if required
    if h_rect.l < 0
        fix_neg = abs(h_rect.l);
    else
        fix_neg = 0;
    end
    
    if h_rect.w < 0
        fix_negw = abs(h_rect.w);
    else
        fix_negw = 0;
    end
    
    for i = 1:h_rect.xnum
        for j = 1:h_rect.ynum
            h_etch.r = h_rect.ETCH_R;
            h_etch.layer = h_rect.etch_layer;
            h_etch.undercut = h_rect.UNDERCUT;
            section.p0 = [h_rect.x-fix_negw+(i-1)*(h_rect.xspace+h_rect.w) h_rect.y-fix_neg+(j-1)*(h_rect.yspace+h_rect.l)];
            section.type = 'rect';
            section.w = abs(h_rect.w);
            section.l = abs(h_rect.l);
            h_etch.circle_etch = h_rect.circle_etch;
            h_etch.rotation_theta = h_rect.theta;
            h_etch.rotation_center = h_rect.p0;
            
            h_etch.regions{1,1} = section;
            EH_1 = etch_hole(h_etch);
            
            %ca = join_gds_structures(ca,EH_1);
            ca = join_gds_structures(ca,EH_1{1});
        end
    end
end


if(h_rect.rp == 1)
    out = points2;
else
    %out = gds_structure(str_name,ca);
    out = ca;
end