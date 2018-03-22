function out = cavity(h_cavity)
% Function to create a backside cavity
% p0 = center of
% r = [r1 r2] inner and outer radii


if ~isfield(h_cavity,'chiplet')
    h_cavity.chiplet = 0;
end

if ~isfield(h_cavity,'zif')
    h_cavity.zif = 0;
end

if ~isfield(h_cavity,'noham')
    h_cavity.noham = 0;
end

if ~isfield(h_cavity,'notether')
    h_cavity.notether = 0;
end


h_cavity.n = 100;
h_cavity.strut_w = 5;



%Create rotor-like structure
points = zeros(2*(h_cavity.n+1),2);
for j=1:2
    if j==2
        mult = -1;
    else
        mult = 1;
    end
    for i=1:h_cavity.n+1
        %Adding 50nm to each radius to fix rounding errors
        points(i + (j-1)*(h_cavity.n+1),1) = h_cavity.p0(1) + (h_cavity.r(j)+.05)*cos(mult*i*2*pi/h_cavity.n);
        points(i + (j-1)*(h_cavity.n+1),2) = h_cavity.p0(2) + (h_cavity.r(j)+.05)*sin(mult*i*2*pi/h_cavity.n);
    end
end
cavity_pts = [points(end,:);points];
%outer_rotor_ring = cavity_pts(h_latch.n+1:end,:);

%be = gds_element('boundary', 'xy',cavity_pts,'layer',2);
%str_name = sprintf('Rotor_[%d,%d],[%d]',round(h_cavity.p0(1)),round(h_cavity.p0(2)),round(h_cavity.r(1)));
%rot = gds_structure(str_name,be);


%Dummy fill backside etch 
be = gds_element('boundary', 'xy',cavity_pts,'layer',8);
str_name = sprintf('R_df_[%d,%d],[%d]',round(h_cavity.p0(1)),round(h_cavity.p0(2)),round(h_cavity.r(1)));
rot = gds_structure(str_name,be);


%Dummy fill backside etch 
%h_cir.x = h_cavity.p0(1);
%h_cir.y = h_cavity.p0(2);
%h_cir.r = h_cavity.r(2);
%h_cir.layer = 8;
%h_cir.n = 100;
%df_backside_cav = circle(h_cir);

%Add backside hole
h_cir.x = h_cavity.p0(1);
h_cir.y = h_cavity.p0(2);
h_cir.r = h_cavity.r(2);
h_cir.layer = 5;
h_cir.n = 100;
backside_cav = circle(h_cir);


%Adds struts to the plate

if h_cavity.nstruts>0
    ang_step = 2*pi/h_cavity.nstruts;
    
    theta_0=0;
    snap = .3;
    holder = cell(1,h_cavity.nstruts);
    for i=1:h_cavity.nstruts
        theta = theta_0 + i*ang_step;
        p1 = h_cavity.p0 + (h_cavity.r(1)-snap)*[cos(theta) sin(theta)];
        p2 = h_cavity.p0 + (h_cavity.r(2)+snap)*[cos(theta) sin(theta)];
        pts = [p1;p2];
        
        pathh = gds_element('path', 'xy',pts,'width', h_cavity.strut_w,'layer',6);
        str_name = sprintf('p_[%d,%d]',round(pts(1,1)),round(pts(1,2)));
        holder{i} = gds_structure(str_name,pathh);
    end
end

%% To create a cavity chiplet

if h_cavity.chiplet
    %Total size of cavity
    if h_cavity.zif
        t_device = 40;
        t_handle = 550;
        t_box = 2;
        t_uncertainty = 15; %unknown from the handle thickness
        t_backside_recess = 25;
        trench_w = 400;
        
        switch h_cavity.zif
            case 1
                % from hammer -> chiplet_w = 2900;
                chiplet_w = 6100;
                chiplet_l = 3000;
                peg_dist = 3200;
                
                if h_cavity.noham
                    chiplet_w = 4700;
                end
            case 2
                % from hammer -> chiplet_w = 3900;
                chiplet_w = 7100;
                chiplet_l = 4000;
                peg_dist = 4200;
            otherwise
                joey = 1;
        end

        
        peg_w = 350; 
        peg_l = t_device + t_handle + t_box + t_uncertainty;
        

        
        p1=  [h_cavity.p0(1) - chiplet_w/2 - trench_w/2 h_cavity.p0(2)+chiplet_l/2+trench_w/2];
        p2 = [h_cavity.p0(1) + chiplet_w/2 + trench_w/2 h_cavity.p0(2)+chiplet_l/2+trench_w/2];
        p3 = [h_cavity.p0(1) + chiplet_w/2 + trench_w/2 h_cavity.p0(2)-chiplet_l/2-trench_w/2];
        p4 = [h_cavity.p0(1) - chiplet_w/2 - trench_w/2 h_cavity.p0(2)-chiplet_l/2-trench_w/2];
        p5 = p1;
        
        be = gds_element('path', 'xy',[p1;p2;p3;p4;p5],'width', trench_w,'layer',5);
        str_name = sprintf('trnchcav_[%d,%d]',round(p1(1)),round(p1(2)));
        cv2t_cav = gds_structure(str_name,be);
        
        
        %Dummy fill
        be = gds_element('path', 'xy',[p1;p2;p3;p4;p5],'width', trench_w,'layer',8);
        str_name = sprintf('trncav_[%d,%d]%d',round(p1(1)),round(p1(2)),8);
        f_df_3_cav = gds_structure(str_name,be);
        
        %Add tethers
        tether_w = 200;
        pt = midpt(p3,p4,.5);
        h_rect.x = pt(1) - tether_w/2;
        h_rect.y = pt(2) - trench_w/2;
        h_rect.w = tether_w;
        h_rect.l = trench_w;
        h_rect.layer = 4;
        tether_strut = rect(h_rect);
        
        %Add peg holes
        tether_w = 200;
        h_rect.x = h_cavity.p0(1) - peg_dist/2 - peg_w;
        h_rect.y = h_cavity.p0(2) - t_device/2;
        h_rect.w = peg_w;
        h_rect.l = peg_l;
        h_rect.layer = 5;
        peg_hole_1 = rect(h_rect);
        h_rect.l = peg_l + t_backside_recess;
        h_rect.layer = 8;
        peg_hole_1_not_dummy = rect(h_rect);
        
        %Add left hammer
        if ~h_cavity.noham
            h_hammer.p0 = [h_rect.x+peg_w/2 h_rect.y];
            h_hammer.num_springs = 10;
            h_hammer.stages = 1;
            h_hammer.spring_l = 2000;
            h_hammer.p0_to_head = 1;
            
            h_hammer.displacement = 80;     %Microns of travel for the hammer
            h_hammer.hhead_r = h_hammer.displacement - 5;
            
            
            l_hammer=make_hammer(h_hammer);
        end
        
        h_rect.x = h_cavity.p0(1) + peg_dist/2;
        h_rect.y = h_cavity.p0(2) - t_device/2;
        h_rect.w = peg_w;
        h_rect.l = peg_l + t_backside_recess;   %Recess the device layer Si back a bit
        h_rect.layer = 5;
        peg_hole_2 = rect(h_rect);
        h_rect.layer = 8;
        h_rect.l = peg_l + t_backside_recess;
        peg_hole_2_not_dummy = rect(h_rect);
        
        
        
        %Add right hammer
        if ~h_cavity.noham
            h_hammer.p0 = [h_rect.x+peg_w/2 h_rect.y];
            h_hammer.num_springs = 10;
            h_hammer.stages = 1;
            h_hammer.p0_to_head = 1;
            r_hammer=make_hammer(h_hammer);
        end
        

        
        
        
    else
        chiplet_w = h_cavity.r(2)*2 + 1000;
        chiplet_l = h_cavity.r(2)*2 + 1000;
        trench_w = 400;
        
        p1=  [h_cavity.p0(1) - chiplet_w/2 - trench_w h_cavity.p0(2)+chiplet_l/2+trench_w/2];
        p2 = [h_cavity.p0(1) + chiplet_w/2 + trench_w/2 h_cavity.p0(2)+chiplet_l/2+trench_w/2];
        p3 = [h_cavity.p0(1) + chiplet_w/2 + trench_w/2 h_cavity.p0(2)-chiplet_l/2-trench_w/2];
        p4 = [h_cavity.p0(1) - chiplet_w/2 - trench_w/2 h_cavity.p0(2)-chiplet_l/2-trench_w/2];
        p5 = p1 + [trench_w/2 0];
        
        be = gds_element('path', 'xy',[p1;p2;p3;p4;p5],'width', trench_w,'layer',5);
        str_name = sprintf('trnchcav_[%d,%d]',round(p1(1)),round(p1(2)));
        cv2t_cav = gds_structure(str_name,be);
                
        %Dummy fill
        be = gds_element('path', 'xy',[p1;p2;p3;p4;p5],'width', trench_w,'layer',8);
        str_name = sprintf('trncav_[%d,%d]%d',round(p1(1)),round(p1(2)),8);
        f_df_3_cav = gds_structure(str_name,be);
        
        %Add fillet rounding to corners
        fillet_len = 100;
        h_fillet.d = .5;
        h_fillet.layer = 5;
        h_fillet.p0 = p1 + [trench_w -trench_w/2];
        h_fillet.p1 = p1 + [trench_w -trench_w/2] + [fillet_len 0];
        h_fillet.p2 = p1 + [trench_w -trench_w/2] + [0 -fillet_len];
        fillet_1 = fillet(h_fillet);
        h_fillet.layer = 8;
        fillet_1a = fillet(h_fillet);
        
        
        h_fillet.d = .5;
        h_fillet.layer = 5;
        h_fillet.p0 = p2 + [-trench_w/2 -trench_w/2];
        h_fillet.p1 = p2 + [-trench_w/2 -trench_w/2] + [-fillet_len 0];
        h_fillet.p2 = p2 + [-trench_w/2 -trench_w/2] + [0 -fillet_len];
        fillet_2 = fillet(h_fillet);
        h_fillet.layer = 8;
        fillet_2a = fillet(h_fillet);
        
        h_fillet.d = .5;
        h_fillet.layer = 5;
        h_fillet.p0 = p3 + [-trench_w/2 trench_w/2];
        h_fillet.p1 = p3 + [-trench_w/2 trench_w/2] + [-fillet_len 0];
        h_fillet.p2 = p3 + [-trench_w/2 trench_w/2] + [0 fillet_len];
        fillet_3 = fillet(h_fillet);
        h_fillet.layer = 8;
        fillet_3a = fillet(h_fillet);
        
        h_fillet.d = .5;
        h_fillet.layer = 5;
        h_fillet.p0 = p4 + [trench_w/2 trench_w/2];
        h_fillet.p1 = p4 + [trench_w/2 trench_w/2] + [fillet_len 0];
        h_fillet.p2 = p4 + [trench_w/2 trench_w/2] + [0 fillet_len];
        fillet_4 = fillet(h_fillet);
        h_fillet.layer = 8;
        fillet_4a = fillet(h_fillet);
        
        if h_cavity.notether==0
            %Add tethers
            tether_w = 200;
            pt = midpt(p3,p4,.5);
            h_rect.x = pt(1) - tether_w/2;
            h_rect.y = pt(2) - trench_w/2;
            h_rect.w = tether_w;
            h_rect.l = trench_w;
            h_rect.layer = 4;
            tether_strut = rect(h_rect);
        end
    end
end



%% Grab all the GDS structures and arrays of structures

%Find all gds structures
a=whos();
b={};
c = 0;
for i=1:length(a)
    if(strcmp(a(i).class,'gds_structure'))
        c = c+1;
        str = sprintf('b{c} = %s;',a(i).name);
        eval(str);
    elseif(strcmp(a(i).class,'cell'))
        str = sprintf('temp = %s;',a(i).name);
        eval(str);
        if(isempty(temp))
            fprintf('Empty Cell! Something went wrong with %s!!\n',a(i).name)
            break;
        end
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