function out = etch_hole(h_etch)
% Function to create etch holes in a structure
% p0 = bottom left coordinate of the rectangle to add etch holes to
% w = width (in x direction) of rectangle
% l = length (in y direction) of rectangle
% r = etch hole radius
% undercut = the largest distance from etch hole to etch hole
% overhang = [0 0 0 0] - Corresponds to the left, bottom, right or top side
%    of the rectangle. Should this side have etch holes that hang off edge?
% layer = GDS layer to place circle into


% A rectangular object with regularly spaced etch holes
%Assuming ystep and xstep are the same
if ~isfield(h_etch,'noetch')
    h_etch.noetch = 0;
end

if ~isfield(h_etch,'undercut')
    h_etch.undercut = 6;
end

if ~isfield(h_etch,'r')
    h_etch.r = 2;
end

if ~isfield(h_etch,'layer')
    h_etch.layer = 2;
end

if ~isfield(h_etch,'rotation_theta')
    h_etch.rotation_theta = 0;
end

if ~isfield(h_etch,'rotation_center')
    h_etch.rotation_center = [0 0];
end

if ~isfield(h_etch,'circle_etch')   % If set to 1, make circular etch holes
    h_etch.circle_etch = 0;
end


if h_etch.noetch == 0
    
    nregions = length(h_etch.regions);
    
    for k=1:nregions
        
        sec = h_etch.regions{k};
        if isfield(sec,'p0')
            h_etch.p0 = sec.p0;
        end
        if isfield(sec,'w')
            h_etch.w = sec.w;
        end
        if isfield(sec,'l')
            h_etch.l = sec.l;
        end
        
        switch sec.type
            case 'rect' %The shape to put etch holes in is rectangular
%                 step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
%                 
%                 ny = round((h_etch.l - 2*(h_etch.undercut + h_etch.r))/step);
%                 nx = round((h_etch.w - 2*(h_etch.undercut + h_etch.r))/step);
%                 
%                 
%                 ystep = (h_etch.l - 2*(h_etch.undercut + h_etch.r))/ny;
%                 xstep = (h_etch.w - 2*(h_etch.undercut + h_etch.r))/nx;
%                 
%                 if nx == 0
%                     nx = 1;
%                     ny = ny + 1;
%                     xstep = 0;
%                     xy_pts = zeros(nx*ny,2);
%                 else
%                     nx = nx+1;
%                     ny = ny+1;
%                     xy_pts = zeros(nx*ny,2);
%                     
%                     for i=1:nx
%                         for j=1:ny
%                             xy_pts((i-1)*ny+j,1) = h_etch.p0(1) + h_etch.undercut + h_etch.r + xstep*(i-1);
%                             xy_pts((i-1)*ny+j,2) = h_etch.p0(2) + h_etch.undercut + h_etch.r + ystep*(j-1);
%                         end
%                     end
%                 end


                % Find x coordinates
                xcoords = etch_points(h_etch,h_etch.w);
                ycoords = etch_points(h_etch,h_etch.l);
                
                count = 1;
                for i=1:numel(xcoords)
                    for j=1:numel(ycoords)
                        xy_pts(count,1) = h_etch.p0(1) + xcoords(i);
                        xy_pts(count,2) = h_etch.p0(2) + ycoords(j);
                        count = count + 1;
                    end
                end
                if h_etch.circle_etch % makes circular etch holes
                    h_cir.r = h_etch.r;
                    h_cir.n = 25;
                    h_cir.layer = h_etch.layer;
                    
                    % Rotate points
                    rot.pts = xy_pts;
                    rot.theta = h_etch.rotation_theta;
                    rot.p0 = h_etch.rotation_center;
                    h_cir.xy =rotate_pts(rot);
                    
                    h_cir.compact = 1;
                    b_circ = circle(h_cir);
                    str = sprintf('k%d = b_circ;',k);
                    eval(str);
                    h_cir.compact = 0;
                    clear xy_pts b_circ
                else % makes square etch holes
                    h_square.r = h_etch.r;
                    h_square.layer = h_etch.layer;
                    
                    % Rotate points
                    rot.pts = xy_pts;
                    rot.theta = h_etch.rotation_theta;
                    rot.p0 = h_etch.rotation_center;
                    h_square.xy =rotate_pts(rot);
                    h_square.rotation_theta = h_etch.rotation_theta;
                    h_square.rotation_center = h_etch.rotation_center;
                    b_circ = etch_square(h_square);
                    str = sprintf('k%d = b_circ;',k);
                    eval(str);
                    h_cir.compact = 0;
                    clear xy_pts b_circ
                end
            case 'rcurve' %Shape has a curve on the right hand side
                
                
                %For rotation
                if ~isfield(h_etch,'rotation_angle')
                    h_etch.rotation_angle = 0;
                end
                
                if ~isfield(h_etch,'rotation_point')
                    h_etch.rotation_point = [0 0];
                end
                
                
                
                % Updated to be compact
                step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                h_cir.r = h_etch.r;
                
                ny = round((h_etch.l - 2*(h_etch.undercut + h_etch.r))/step);
                ystep = (h_etch.l - 2*(h_etch.undercut + h_etch.r))/ny;
                ny = ny+1;
                
                xy_pts = zeros(ny*ny,2);
                count = 0;
                for i=1:ny
                    %Y coordinate of the etch hole
                    h_cir.y = h_etch.p0(2) + h_etch.undercut + h_etch.r + ystep*(i-1);
                    
                    %At each y slice need to find width
                    yvals = sec.rcurve(:,2);
                    
                    [val1 index1] = min(abs(yvals-h_cir.y));
                    [val2 index2] = min(abs(yvals-h_cir.y-h_cir.r));
                    [val3 index3] = min(abs(yvals-h_cir.y+h_cir.r));
                    
                    wid = min([sec.rcurve(index1,1) sec.rcurve(index2,1) sec.rcurve(index3,1)])-h_etch.p0(1);
                    
                    if(abs(sec.rcurve(index2,1) - sec.rcurve(index3,1))>30)
                        wid = wid - 10;
                    end
                    
                    if wid<30
                        nx = round((wid - (h_etch.undercut + h_etch.r))/step);
                        xstep = (wid - 2*(.5*h_etch.undercut + h_etch.r))/nx;
                        
                        %nx = nx+1;
                        
                        for j=1:nx
                            count = count+1;
                            xy_pts(count,1) = h_etch.p0(1) + j*wid/(nx+1);
                            xy_pts(count,2) = h_cir.y;
                        end
                    else
                        nx = round((wid - 2*(h_etch.undercut + h_etch.r))/step);
                        xstep = (wid - 2*(h_etch.undercut + h_etch.r))/nx;
                        %if nx>1
                         nx = nx+1; %testing with >1, usually always done regardless of nx value
                        %end
                        
                        for j=1:nx
                            count = count+1;
                            xy_pts(count,1) = h_etch.p0(1) + h_etch.undercut + h_etch.r + xstep*(j-1);
                            xy_pts(count,2) = h_cir.y;
                        end
                    end
                end
                %remove extra points
                xy_pts(count+1:end,:) = [];
                
                h_cir.r = h_etch.r;
                h_cir.n = 25;
                h_cir.layer = h_etch.layer;
                h_cir.compact = 1;
                %h_cir.xy = xy_pts;
                h_rot.p0 = h_etch.rotation_point;
                h_rot.theta = h_etch.rotation_angle;
                h_rot.pts = xy_pts;
                h_cir.xy =rotate_pts(h_rot);
                
                str = sprintf('k%d = circle(h_cir);',k);
                eval(str);
                h_cir.compact = 0;
                clear xy_pts
            case 'lcurve' %Shape has a curve on the left hand side
                % Updated to be compact
                
                
                %For rotation
                if ~isfield(h_etch,'rotation_angle')
                    h_etch.rotation_angle = 0;
                end
                
                if ~isfield(h_etch,'rotation_point')
                    h_etch.rotation_point = [0 0];
                end
                
                step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                h_cir.r = h_etch.r;
                
                ny = round((h_etch.l - 2*(h_etch.undercut + h_etch.r))/step);
                ystep = (h_etch.l - 2*(h_etch.undercut + h_etch.r))/ny;
                ny = ny+1;
                
                xy_pts = zeros(ny*ny,2);
                count = 0;
                for i=1:ny
                    %Y coordinate of the etch hole
                    h_cir.y = h_etch.p0(2) + h_etch.undercut + h_etch.r + ystep*(i-1);
                    
                    %At each y slice need to find width
                    yvals = sec.rcurve(:,2);
                    
                    [val1 index1] = min(abs(yvals-h_cir.y));
                    [val2 index2] = min(abs(yvals-h_cir.y-h_cir.r));
                    [val3 index3] = min(abs(yvals-h_cir.y+h_cir.r));
                    
                    wid = min([h_etch.p0(1)-sec.rcurve(index1,1) h_etch.p0(1)-sec.rcurve(index2,1) h_etch.p0(1)-sec.rcurve(index3,1)]);
                    %wid = min([sec.rcurve(index1,1) sec.rcurve(index2,1) sec.rcurve(index3,1)])-h_etch.p0(1);
                    
                    if(abs(sec.rcurve(index2,1) - sec.rcurve(index3,1))>30)
                        wid = wid - 30;
                    end
                    
                    if wid<30
                        nx = round((wid - (h_etch.undercut + h_etch.r))/step);
                        xstep = (wid - 2*(.5*h_etch.undercut + h_etch.r))/nx;
                        
                        for j=1:nx
                            count = count+1;
                            xy_pts(count,1) = h_etch.p0(1) - j*wid/(nx+1);
                            xy_pts(count,2) = h_cir.y;
                        end
                    else
                        nx = round((wid - 2*(h_etch.undercut + h_etch.r))/step);
                        xstep = (wid - 2*(h_etch.undercut + h_etch.r))/nx;
                        nx = nx+1;
                        
                        for j=1:nx
                            count = count+1;
                            
                            xy_pts(count,1) = h_etch.p0(1) - h_etch.undercut - h_etch.r - xstep*(j-1);
                            
                            %xy_pts(count,1) = h_etch.p0(1) - j*wid/(nx+1);
                            xy_pts(count,2) = h_cir.y;
                        end
                    end
                end
                %remove extra points
                xy_pts(count+1:end,:) = [];
                
                h_cir.r = h_etch.r;
                h_cir.n = 25;
                h_cir.layer = h_etch.layer;
                h_cir.compact = 1;
                %h_cir.xy = xy_pts;
                h_rot.p0 = h_etch.rotation_point;
                h_rot.theta = h_etch.rotation_angle;
                h_rot.pts = xy_pts;
                h_cir.xy =rotate_pts(h_rot);
                
                str = sprintf('k%d = circle(h_cir);',k);
                eval(str);
                h_cir.compact = 0;
                clear xy_pts
                
            case 'semicircle_t' %Shape is a circle, and needs etch holes
                step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                %ny = round((sec.r - 2*(h_etch.undercut + h_etch.r))/step);
                ny = round((sec.r - 2*(.5*h_etch.undercut + h_etch.r))/step);
                if ny==0
                    ystep = (sec.r - 2*(.5*h_etch.undercut + h_etch.r));
                else
                    ystep = (sec.r - 2*(.5*h_etch.undercut + h_etch.r))/ny;
                end
                
                ny = ny+1;
                
                %nx = round((h_etch.w - 2*(h_etch.undercut + h_etch.r))/step);
                %xstep = (h_etch.w - 2*(h_etch.undercut + h_etch.r))/nx;
                %nx = nx+1;
                
                b_circ={}; %This is not a great way to initialize this...
                count = 0;
                for i=1:ny
                    %Y coordinate of the etch hole
                    %h_cir.y = h_etch.p0(2) + h_etch.undercut + h_etch.r + ystep*(i-1);
                    h_cir.y = h_etch.p0(2) + h_etch.r + ystep*(i-1);
                    h_cir.r = h_etch.r;
                    h_cir.n = 25;
                    h_cir.layer = h_etch.layer;
                    
                    %At each y slice need to find width
                    wid = 2*sqrt((sec.r)^2-(h_cir.r + h_cir.y-h_etch.p0(2))^2);
                    
                    nx = round((wid - 2*(h_etch.undercut + h_etch.r))/step);
                    xstep = (wid - 2*(h_etch.undercut + h_etch.r))/nx;
                    nx = nx+1;
                    
                    for j=1:nx
                        h_cir.x = h_etch.p0(1) - sec.r + h_etch.undercut + h_etch.r + xstep*(j-1) + (2*sec.r - wid)/2;
                        count = count+1;
                        b_circ{count} = circle(h_cir);
                        %fprintf('x: %.1f y: %.1f\n',h_cir.x,h_cir.y)
                    end
                end
                str = sprintf('k%d = b_circ;',k);
                eval(str);
                clear b_circ
            case 'annulus' %Shape is an annulus, and needs etch holes
                % Updated to be compact
                % Initial radius of ring
                init_rad = sec.r(2) - h_etch.undercut - h_etch.r;
                
                %Figure out space between each etch hole in a ring
                cent_to_cent = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                nrings = round((sec.r(2)-sec.r(1) - 2*(h_etch.undercut + h_etch.r))/cent_to_cent);
                radial_step = (sec.r(2)-sec.r(1) - 2*(h_etch.undercut + h_etch.r))/nrings;
                nrings = nrings + 1;
                xy_pts=zeros(100*nrings,2); %find better way to initialize this!
                count = 0;
                for j=1:nrings
                    %Find circumference one undercut away from outer radius edge
                    cur_rad = init_rad - radial_step*(j-1);
                    
                    circ = 2*pi*(cur_rad); %circumference of ring
                    
                    %Find number of etch holes per ring
                    n = round(circ/cent_to_cent);
                    for i=1:n
                        %Y coordinate of the etch hole
                        count = count + 1;
                        xy_pts(count,1) =  sec.p0(1) + cur_rad*cos(i*2*pi/n);
                        xy_pts(count,2) =  sec.p0(2) + cur_rad*sin(i*2*pi/n);
                    end
                end
                xy_pts(count+1:end,:)=[];
                h_cir.r = h_etch.r;
                h_cir.n = 25;
                h_cir.layer = h_etch.layer;
                h_cir.xy = xy_pts;
                h_cir.compact = 1;
                str = sprintf('k%d = circle(h_cir);',k);
                eval(str);
                h_cir.compact = 0;
            case 'partial_annulus' %Shape is a partial annulus, such as the C of a joint
                % Updated to be compact
                
                %For rotation
                if ~isfield(h_etch,'rotation_angle')
                    h_etch.rotation_angle = 0;
                end
                
                if ~isfield(h_etch,'rotation_point')
                    h_etch.rotation_point = [0 0];
                end
                
                
                % Initial radius of ring
                init_rad = sec.r(2) - h_etch.undercut - h_etch.r;
                
                %Figure out space between each etch hole in a ring
                cent_to_cent = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                nrings = round((sec.r(2)-sec.r(1) - 2*(h_etch.undercut + h_etch.r))/cent_to_cent);
                radial_step = (sec.r(2)-sec.r(1) - 2*(h_etch.undercut + h_etch.r))/nrings;
                nrings = nrings + 1;
                xy_pts=zeros(100*nrings,2); %find better way to initialize this!
                count = 0;
                for j=1:nrings
                    %Find circumference one undercut away from outer radius edge
                    cur_rad = init_rad - radial_step*(j-1);
                    
                    circ = 2*pi*(cur_rad); %circumference of ring (should be arc length of partial ring)
                    circ = sec.theta(2)*(cur_rad); %Arc length of ring at current radius
                    
                    %Find number of etch holes per ring
                    n = round(circ/cent_to_cent);
                    for i=1:n
                        %Y coordinate of the etch hole
                        count = count + 1;
                        xy_pts(count,1) =  sec.p0(1) + cur_rad*cos(sec.theta(1)+(i-1)*sec.theta(2)/n + h_etch.undercut/(2*cur_rad));
                        xy_pts(count,2) =  sec.p0(2) + cur_rad*sin(sec.theta(1)+(i-1)*sec.theta(2)/n + h_etch.undercut/(2*cur_rad));
                    end
                end
                xy_pts(count+1:end,:)=[];
                h_cir.r = h_etch.r;
                h_cir.n = 25;
                h_cir.layer = h_etch.layer;
                %rotate these points
                h_rot.p0 = h_etch.rotation_point;
                h_rot.theta = h_etch.rotation_angle;
                h_rot.pts = xy_pts;
                h_cir.xy =rotate_pts(h_rot);
                h_cir.compact = 1;
                str = sprintf('k%d = circle(h_cir);',k);
                eval(str);
                h_cir.compact = 0;
            case 'arc' %Shape is an arc (such as the hammer bird's beak)
                %Walk along arc until find a spot where one etch hole can go
                circ_to_arc = 2; %minimum distance from edge of etch hole to arced edge of beak
                dist = 0;
                count = 1;
                while dist<circ_to_arc
                    dist = abs(sqrt((sec.arc(count,1)-sec.p0(1))^2 + (sec.arc(count,2)-sec.p0(2))^2)-sec.r(2)) - h_etch.r;
                    count = count + 1;
                end
                
                %starting and ending angles
                theta1 = atan2((sec.arc(count,2)-sec.p0(2)),(sec.arc(count,1)-sec.p0(1)));
                if theta1<0
                    theta1=theta1+2*pi;
                end
                theta2 = pi - sec.theta2;
                
                %Additional angle to help get etch holes as close to fillet as
                %possible
                theta3 = sec.bbfillet/sec.r(2);
                
                %Find number of steps in total angular span
                step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                theta_step = step/sec.r(2);
                
                temp_theta = abs(theta2 - theta1) + theta3;
                
                
                %while(temp_theta>2*pi)
                %    temp_theta = temp_theta-2*pi;
                %end
                
                if temp_theta>3*pi/2
                    temp_theta = temp_theta - 2*pi;
                end
                
                n_steps = abs(round(temp_theta/theta_step));
                %theta_step = abs(theta2 - theta1)/n_steps;
                theta_step = abs(temp_theta)/n_steps;
                
                
                %Create array that has angle of each (x,y) coordinate of arc
                sec.angles = atan2((sec.arc(:,2)-sec.p0(2)),(sec.arc(:,1)-sec.p0(1)));
                for jj=1:length(sec.angles)
                    if sec.angles(jj)<0
                        sec.angles(jj) =  sec.angles(jj) + 2*pi;
                    end
                end
                sec.radius = sqrt((sec.arc(:,2)-sec.p0(2)).^2 + (sec.arc(:,1)-sec.p0(1)).^2);
                
                find_pts = [];
                
                sec.angles = unwrap(sec.angles);
                
                %This is terrible. Needed to do it to make etch holes in
                %compact ZIF hammer chiplet work in a pinch. Fix this ASAP.
                
                if n_steps<5
                    n_steps = n_steps + 1;
                end
                
                for j=1:n_steps
                    %for j=1:3
                    %Find closest angle to theta1
                    theta = theta1 +  sec.orientation*(j-1)*theta_step;
                    
                    [minval index] = min(abs(sec.angles-theta));
                    
                    %Determine how many etch holes to put in reach radial line
                    t_dist = sec.radius(index) -  sec.r(2) + h_etch.undercut - h_etch.r;
                    
                    n_radial = round(t_dist/(h_etch.undercut + 2*h_etch.r));
                    
                    pitch = t_dist/(n_radial + 1);
                    
                    
                    
                    b_circ=cell(1,n_radial); %This is not a great way to initialize this...
                    
                    for i=1:n_radial
                        %h_cir.y = sec.p0(2) + (sec.r(2) - h_etch.undercut + h_etch.r + (i)*pitch)*sin(sec.angles(index));
                        %h_cir.x = sec.p0(1) + (sec.r(2) - h_etch.undercut + h_etch.r + (i)*pitch)*cos(sec.angles(index));
                        h_cir.y = sec.p0(2) + (sec.r(2) - h_etch.undercut - h_etch.r + i*(pitch + h_etch.r))*sin(sec.angles(index));
                        h_cir.x = sec.p0(1) + (sec.r(2) - h_etch.undercut - h_etch.r + i*(pitch + h_etch.r))*cos(sec.angles(index));
                        find_pts = [find_pts;h_cir.x h_cir.y];
                        h_cir.r = h_etch.r;
                        h_cir.n = 25;
                        h_cir.layer = h_etch.layer;
                        b_circ{i} = circle(h_cir);
                    end
                    
                    str = sprintf('k%d = b_circ;',j);
                    eval(str);
                    clear b_circ
                end
            case 'tcurve' %Shape has a curve on the top side
                % Updated to be compact
                
                if ~isfield(h_etch,'rotation_angle')
                    h_etch.rotation_angle = 0;
                end
                
                if ~isfield(h_etch,'rotation_point')
                    h_etch.rotation_point = [0 0];
                end
                
                step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                nx = round((h_etch.w - 2*(h_etch.undercut + h_etch.r))/step);
                xstep = (h_etch.w - 2*(h_etch.undercut + h_etch.r))/nx;
                nx = nx+1;
                
                h_rot.p0 = h_etch.rotation_point;
                h_rot.theta = -(h_etch.rotation_angle)*pi/180;
                
                xy_pts = zeros(100*nx,2); %Find better way to initialize!!
                count = 0;
                for i=1:nx
                    %X coordinate of the etch hole
                    h_cir.x = h_etch.p0(1) + h_etch.undercut + h_etch.r + xstep*(i-1);
                    
                    
                    %At each x slice need to find length
                    xvals = sec.tcurve(:,1);
                    
                    [val1 index1] = min(abs(xvals-h_cir.x));
                    [val2 index2] = min(abs(xvals-h_cir.x-h_etch.r));
                    [val3 index3] = min(abs(xvals-h_cir.x+h_etch.r));
                    
                    %len = min([sec.tcurve(index1,2) sec.tcurve(index2,2) sec.tcurve(index3,2)])-h_etch.p0(2) + h_etch.undercut;
                    len = min(abs([sec.tcurve(index1,2)-h_etch.p0(2) sec.tcurve(index2,2)-h_etch.p0(2) sec.tcurve(index3,2)-h_etch.p0(2)])) + h_etch.undercut;
                    
                    ny = round((len - 2*(h_etch.undercut + h_etch.r))/step);
                    ystep = (len - 2*(h_etch.undercut + h_etch.r))/ny;
                    ny = ny + 1;
                    h_cir.x_init = h_cir.x;
                    for j=1:ny
                        h_cir.y = h_etch.p0(2) + h_etch.undercut + h_etch.r + ystep*(j-1);
                        h_rot.pts = [h_cir.x_init h_cir.y];
                        rotated =rotate_pts(h_rot);
                        count = count + 1;
                        xy_pts(count,1) = rotated(1);
                        xy_pts(count,2) = rotated(2);
                    end
                end
                xy_pts(count+1:end,:) = [];
                
                h_cir.r = h_etch.r;
                h_cir.n = 25;
                h_cir.layer = h_etch.layer;
                h_cir.xy = xy_pts;
                h_cir.compact = 1;
                str = sprintf('k%d = circle(h_cir);',k);
                eval(str);
                h_cir.compact = 0;
            case 'semicircle_l' %Shape is a circle, and needs etch holes on the left
                
                if ~isfield(h_etch,'rotation_angle')
                    h_etch.rotation_angle = 0;
                end
                
                if ~isfield(h_etch,'rotation_point')
                    h_etch.rotation_point = [0 0];
                end
                
                h_rot.p0 = h_etch.rotation_point;
                h_rot.theta = -(h_etch.rotation_angle)*pi/180;
                
                
                step = (sqrt(2)*h_etch.undercut + sqrt(2)*h_etch.r);
                
                
                %Find total number of x slices (used to be ny)
                nx = round((sec.r - 2*(.5*h_etch.undercut + h_etch.r))/step);
                if nx==0
                    xstep = (sec.r - 2*(.5*h_etch.undercut + h_etch.r));
                else
                    xstep = (sec.r - 2*(.5*h_etch.undercut + h_etch.r))/nx;
                end
                nx = nx+1;
                
                
                b_circ={}; %This is not a great way to initialize this...
                count = 0;
                for i=1:nx
                    %X coordinate of the etch hole
                    h_cir.x = h_etch.p0(1) - h_etch.r - xstep*(i-1);
                    h_cir.x_init = h_cir.x;
                    %h_cir.x = h_etch.p0(1) - xstep*(i-1);
                    h_cir.r = h_etch.r;
                    h_cir.n = 25;
                    h_cir.layer = h_etch.layer;
                    
                    %At each x slice need to find the length
                    %len = 2*sqrt((sec.r)^2-(h_cir.r + h_cir.x-h_etch.p0(1))^2);
                    len = 2*sqrt((sec.r)^2-(h_cir.x-h_etch.p0(1) - h_cir.r)^2);
                    
                    ny = round((len - 2*(h_etch.undercut + h_etch.r))/step);
                    
                    ystep = (len - 2*(h_etch.undercut + h_etch.r))/ny;
                    ny = ny+1;
                    
                    for j=1:ny
                        h_cir.y = h_etch.p0(2) - sec.r + h_etch.undercut + h_etch.r + ystep*(j-1) + (2*sec.r - len)/2;
                        count = count+1;
                        
                        h_rot.pts = [h_cir.x_init h_cir.y];
                        rotated =rotate_pts(h_rot);
                        h_cir.x = rotated(1);
                        h_cir.y = rotated(2);
                        b_circ{count} = circle(h_cir);
                        %fprintf('x: %.1f y: %.1f\n',h_cir.x,h_cir.y)
                    end
                end
                str = sprintf('k%d = b_circ;',k);
                eval(str);
                clear b_circ
            case 'raw_points'  %Feed the function raw points and have it put etch holes on those coordinates
                h_cir.r = h_etch.r;
                h_cir.n = 25;
                h_cir.layer = h_etch.layer;
                h_cir.xy = h_etch.pts;
                h_cir.compact = 1;
                b_circ = circle(h_cir);
                str = sprintf('k%d = b_circ;',k);
                eval(str);
                h_cir.compact = 0;
                clear xy_pts b_circ
                
        end
    end
    
    %% Collect and save all GDS structures to an output variable
    
    %Find all gds structures and save into cell array b
    a=whos();
    b={};
    c = 0;
    for i=1:length(a)
        if(strcmp(a(i).class,'gds_structure'))
            c = c+1;
            str = sprintf('b{c} = %s;',a(i).name);
            eval(str);
        elseif(strcmp(a(i).class,'cell') && a(i).bytes>0)
            str = sprintf('strcmp(class(%s{1}),''gds_structure'');',a(i).name);
            if(eval(str))
                str = sprintf('temp = %s;',a(i).name);
                eval(str)
                for i=1:length(temp)
                    c = c+1;
                    b{c} = temp{i};
                end
            end
        end
    end
    
    % Outputs a cell array of
    out = b;
    
else
    out = 0;
end