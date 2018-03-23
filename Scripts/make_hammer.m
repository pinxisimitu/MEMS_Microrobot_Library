function out=make_hammer(h_hammer)
% Function to create a full hammer with a latch etc
% p0 = bottom left coordinate of the hammer
% w = width (in x direction)
% l = length (in y direction)
% layer = GDS layer to place circle into
% xnum = number of repeats in x direction
% xspace = gap between x repeats
% ynum = number of repeats in y direction
% yspace = gap between y repeats
% theta = angle by which to rotate rectangle
% p0 = point about which to rotate
% rp = return points, if rp==1, returns column vector of verticies

if ~isfield(h_hammer,'stages')
    h_hammer.stages = 1;
end

if ~isfield(h_hammer,'piezo')
    h_hammer.piezo = 0;
end

if ~isfield(h_hammer,'r_stages')
    h_hammer.r_stages = 100;
end

if ~isfield(h_hammer,'no_bottom')
    h_hammer.no_bottom = 0;
end


if ~isfield(h_hammer,'no_bottom')
    h_hammer.no_bottom = 0;
end

if ~isfield(h_hammer,'no_head')
    h_hammer.no_head = 0;
end

if ~isfield(h_hammer,'inchworm')
    h_hammer.inchworm = 0;
end

if ~isfield(h_hammer,'h_joint')
    h_hammer.h_joint = 0;
end

if ~isfield(h_hammer,'mech_latch')
    h_hammer.mech_latch = 1;
end

if ~isfield(h_hammer,'alumina')
    h_hammer.alumina = 0;
end

if ~isfield(h_hammer,'displacement')
    h_hammer.displacement = 40;
end

if ~isfield(h_hammer,'theta_arm')
    h_hammer.theta_arm = -90;
end


if ~isfield(h_hammer,'spring_gap')
    h_hammer.spring_gap = 10;
end


if ~isfield(h_hammer,'noetch')
    h_hammer.noetch = 0;
end

if ~isfield(h_hammer,'closed')
    h_hammer.closed = 0;
end

if ~isfield(h_hammer,'chiplet')
    h_hammer.chiplet = 0;
end

if ~isfield(h_hammer,'large_dx')
    h_hammer.large_dx = 0;
end

if ~isfield(h_hammer,'dage')
    h_hammer.dage = 0;
end

if ~isfield(h_hammer,'cavity_chiplet')
    h_hammer.cavity_chiplet = 0;
end


if ~isfield(h_hammer,'spring_l')
    h_hammer.spring_l = 1000;
end

if isfield(h_hammer,'ham_w')
    h_ham.ham_w = h_hammer.ham_w;
end

if isfield(h_hammer,'anchor_w')
    h_ham.anchor_w = h_hammer.anchor_w;
else
    h_hammer.anchor_w = 50;
end

if ~isfield(h_hammer,'latch_actuation_angle')
    h_hammer.latch_actuation_angle = 40;
end

if ~isfield(h_hammer,'latch_initial_angle')
    h_hammer.latch_initial_angle = 40;
end

if ~isfield(h_hammer,'chiplet_w')
    h_hammer.chiplet_w = 2000;
end

if ~isfield(h_hammer,'rotor_gap')
    h_hammer.rotor_gap = 3;
end


if ~isfield(h_hammer,'ss_contact_w')
    h_hammer.ss_contact_w = 7;
end

if ~isfield(h_hammer,'ss_meander_spline')
    h_hammer.ss_meander_spline = 0;
end

if ~isfield(h_hammer,'no_backstops')
    h_hammer.no_backstops = 0;
end

if ~isfield(h_hammer,'compact_latch')
    h_hammer.compact_latch = 0;
end

if ~isfield(h_hammer,'hhead_r')
    h_hammer.hhead_r = 40;
end

if ~isfield(h_hammer,'p0_to_head')
    h_hammer.p0_to_head = 0;
end

if ~isfield(h_hammer,'zif')
    h_hammer.zif = 0;
end


if ~isfield(h_hammer,'h_latch')
        h_hammer.h_latch.ss.n = 6;
        h_hammer.h_latch.ss.dpp = 70;
        h_hammer.h_latch.ss.dist_from_rotor = 200;
end

if ~isfield(h_hammer,'rp')
    h_hammer.rp = 0;
end

if ~isfield(h_hammer,'h_latch_r')
    h_hammer.h_latch_r = [50 100];
end


%Creates hammer
h_ham.x = h_hammer.p0(1);                             % x coordinate of hammer anchor
h_ham.y = h_hammer.p0(2);                             % y coordinate of hammer anchor
h_ham.num_springs = h_hammer.num_springs;             % Number of springs on each side of hammer
h_ham.spring_l = h_hammer.spring_l/2;                 % Half Length of the springs
h_ham.spring_w = 6;                                   % Width of the springs
h_ham.chiplet = h_hammer.chiplet;
h_ham.spring_gap = h_hammer.spring_gap;
h_ham.stages = h_hammer.stages;
h_ham.displacement = h_hammer.displacement;
h_ham.large_dx = h_hammer.large_dx;
h_ham.dage = h_hammer.dage;
h_ham.no_bottom = h_hammer.no_bottom;
h_ham.anchor_w=h_hammer.anchor_w;
h_ham.hhead_r = h_hammer.hhead_r;
h_ham.zif = h_hammer.zif;  
h_ham.ham_w = 400;
h_ham.no_bottom = h_hammer.no_bottom;
h_ham.piezo = h_hammer.piezo;

%Radius of hammer head

h_ham.noetch = h_hammer.noetch;
h_ham.no_head = h_hammer.no_head;
[g_hammer ref_points] = ham(h_ham);
latch_point = ref_points(1,:);
backside_etch = ref_points(2,:) - [0 3*h_hammer.displacement];
backside_etch_anchor = ref_points(3,:);
fill_body_contact = ref_points(4,:) - [50 h_hammer.displacement];
fill_ham_anchor_t = ref_points(5,:);
fill_ham_anchor_b = ref_points(6,:);
fill_ham_head_cent = ref_points(7,:);            %The center of the hammer head

if h_hammer.p0_to_head
    clear g_hammer
    h_ham.x = h_hammer.p0(1) - (fill_ham_head_cent(1)-h_hammer.p0(1));                             % x coordinate of hammer anchor
    h_ham.y = h_hammer.p0(2) - (fill_ham_head_cent(2)-h_hammer.p0(2));
    [g_hammer ref_points] = ham(h_ham);
    latch_point = ref_points(1,:);
    backside_etch = ref_points(2,:) - [0 3*h_hammer.displacement];
    backside_etch_anchor = ref_points(3,:);
    fill_body_contact = ref_points(4,:) - [50 h_hammer.displacement];
    fill_ham_anchor_t = ref_points(5,:);
    fill_ham_anchor_b = ref_points(6,:);
    fill_ham_head_cent = ref_points(7,:);            %The center of the hammer head
end


%Creates a latch
h_latch.orientation = -1;

h_latch.r = h_hammer.h_latch_r;

h_latch.ss.n = h_hammer.h_latch.ss.n;
h_latch.ss.dpp = h_hammer.h_latch.ss.dpp;
h_latch.ss.dist_from_rotor = h_hammer.h_latch.ss.dist_from_rotor;
h_latch.pin_gap = h_hammer.rotor_gap; 
h_latch.h_joint = h_hammer.h_joint;
h_latch.inchworm = h_hammer.inchworm; 
updated_p0 = latch_point + [h_latch.r(2)+2 -h_hammer.displacement];
h_latch.actuation_angle =  h_hammer.latch_actuation_angle;               % Total angle the
h_latch.init_angle = h_hammer.latch_initial_angle;                    % Angle above the negative x axis
h_latch.ss_contact_w = h_hammer.ss_contact_w;
h_latch.ss_meander_spline = h_hammer.ss_meander_spline;


%Variables used yo determine how big the chiplet will be
layout_right = -25000;
layout_bottom = 25000;


if h_hammer.stages >0
    for k = 1:h_hammer.stages
        if k <= h_hammer.r_stages
            h_latch.r = h_hammer.h_latch_r;
        else
            h_latch.r = [50 100];
        end

        
        
        h_latch.bbeak_fillet_length = 5;
        h_latch.p0 = updated_p0;
        
        %Removes backstops on latches that don't have voltage on them
        if k<h_hammer.stages
            h_latch.no_backstops = h_hammer.no_backstops;
        else
            h_latch.no_backstops = 0;
        end
        h_latch.n = 100;
        h_latch.arml =900;
        h_latch.armw = 2*h_latch.r(1);
        
        h_latch.hhead_r = h_ham.hhead_r;
        h_latch.closed = h_hammer.closed;
        h_latch.alumina = h_hammer.alumina;
        
        h_latch.orientation = h_latch.orientation*-1;
        h_latch.compact_latch = h_hammer.compact_latch;
        if h_latch.compact_latch
            switch k
                case 1
                    h_latch.theta = h_latch.init_angle;
                    h_latch.theta_arm = -121;
                    h_latch.manual_theta = 0;
                    h_latch.compact_ang1 = 197;
                    h_latch.flip_spring = 0;
                case 2
                    h_latch.theta = (180-h_latch.actuation_angle + h_latch.theta_arm)-h_latch.actuation_angle;
                    h_latch.theta_arm = -102;
                    h_latch.manual_theta = 1.2*pi;
                    h_latch.compact_ang1 = 197;
                    h_latch.flip_spring = 1;
                case 3
                    h_latch.theta = 180 + h_latch.theta_arm + h_latch.actuation_angle + h_latch.actuation_angle;
                    h_latch.theta_arm = 40;
                    h_latch.manual_theta = 1.3*pi;
                    h_latch.compact_ang1 = 207;
                    h_latch.flip_spring = 1;
                case 4
                    h_latch.theta = (180-h_latch.actuation_angle + h_latch.theta_arm)-h_latch.actuation_angle;
                    h_latch.theta_arm = 20;
                    h_latch.compact_ang1 = 35;
                    h_latch.manual_theta = 1.1*pi;
                    h_latch.flip_spring = 0;
                otherwise
                    %Run the standard non-optimized latches
                    if mod(k,2)==1
                        h_latch.theta = h_latch.init_angle;
                        h_latch.theta_arm = h_hammer.theta_arm;
                    else
                        h_latch.theta = (90-h_latch.init_angle)-h_latch.actuation_angle;
                        h_latch.theta_arm = 180 - h_latch.actuation_angle;
                    end
            end
                    
        else
            %Standard non-optimized latches
            if mod(k,2)==1
                h_latch.theta = h_latch.init_angle;
                h_latch.theta_arm = h_hammer.theta_arm;
            else
                h_latch.theta = (90-h_latch.init_angle)-h_latch.actuation_angle;
                h_latch.theta_arm = 180 - h_latch.actuation_angle;
            end
        end
        h_latch.chiplet = h_hammer.chiplet;
        h_latch.layer = 6;
        h_latch.zif = h_hammer.zif;
        h_latch.noetch = h_hammer.noetch;
        h_latch.stages = h_hammer.stages;
        h_latch.cstage = k;
        h_latch.mech_latch = h_hammer.mech_latch;
        ev = sprintf('[g_latch_%d latch_ref_pts]= ham_latch(h_latch);',k);
        eval(ev);
        
        %Find the right and bottom most points
        ev = sprintf('[l_min r_max t_max b_min] = find_bb(g_latch_%d);',k);
        eval(ev);
        
        layout_right = max([layout_right r_max]);
        layout_bottom = min([layout_bottom b_min]);
        
        if h_hammer.zif && ~h_hammer.compact_latch
            layout_right = layout_right + 215;
        end
        
               
       
        h_rot.pts = latch_ref_pts(4,:);
        h_rot.p0 = h_latch.p0;
        h_rot.theta = h_latch.orientation*h_latch.actuation_angle*pi/180;
        updated_p0 = rotate_pts(h_rot);
        %updated_p0 = latch_ref_pts(4,:);
    end
    
    backside_trench_len = latch_ref_pts(1,2);
    fill_gap_stop_1 = latch_ref_pts(2,:);
    fill_gap_stop_2 = latch_ref_pts(3,:);
    
end

%Create backside square around it all and all the fill
if (h_hammer.chiplet ==1)
    trench_w = 400; %Width of the trench around the hammer chiplets
    h_rect.x = h_ham.x;
    h_rect.y = h_ham.y;
    h_rect.w = sqrt((h_ham.x-backside_etch(1))^2 + (h_ham.y-backside_etch(2))^2) + 2*h_latch.r(2) + trench_w + 240 +60*(h_hammer.num_springs/5);
    h_rect.l = backside_trench_len + 150 + trench_w;
    h_rect.theta = -atan((h_ham.y-backside_etch(2))/(-h_ham.x+backside_etch(1)));
    tether_ang1 = h_rect.theta;
    h_rect.phi = pi/2 + h_rect.theta;
    tether_ang2 = h_rect.phi;
    dist = sqrt( (backside_etch_anchor(1) - h_ham.x)^2 + (backside_etch_anchor(2) - h_ham.y)^2);
    elongate = (dist)*cos(pi/2+h_rect.theta)+trench_w/2;
    h_rect.x = h_ham.x - elongate*cos(h_rect.theta) - trench_w/2*cos(h_rect.phi);
    h_rect.y = h_ham.y - elongate*sin(h_rect.theta) - trench_w/2*sin(h_rect.phi);
    
    
    h_rect.p0 = [h_rect.x h_rect.y];
    h_rect.layer = 1;
    h_rect.rp = 1;
    backside_cut_points = rect(h_rect);
    backside_cut_points = [backside_cut_points;backside_cut_points(1,:)];
    h_rect.rp = 0;
    elongate = dist*cos(pi/2+h_rect.theta);
    h_rect.x = h_ham.x - elongate*cos(h_rect.theta);
    h_rect.y = h_ham.y - elongate*sin(h_rect.theta);
    h_rect.p0 = [h_rect.x h_rect.y];
    h_rect.w = sqrt((h_ham.x-backside_etch(1))^2 + (h_ham.y-backside_etch(2))^2) + 2*h_latch.r(2) + 300;
    h_rect.l = backside_trench_len + 150;
    
    %backside_inside = rect(h_rect);
    h_rect.rp = 1;
    backside_inside_points = rect(h_rect);
    
    be = gds_element('path', 'xy',backside_cut_points,'width', trench_w,'layer',5);
    str_name = sprintf('trn_[%d,%d]%d',round(h_rect.x),round(h_rect.y),h_rect.layer);
    backside_cut = gds_structure(str_name,be);
    
    %Dummy fill    
    be = gds_element('path', 'xy',backside_cut_points,'width', trench_w,'layer',8);
    str_name = sprintf('trn_[%d,%d]%d',round(h_rect.x),round(h_rect.y),8);
    f_df = gds_structure(str_name,be);
    
    be = gds_element('boundary', 'xy',backside_cut_points,'layer',8);
    str_name = sprintf('df_[%d,%d]%d',round(h_rect.x),round(h_rect.y),8);
    f_df_2 = gds_structure(str_name,be);
    
    
    
    
    %Add fill to the chiplets
    
    %triangle at left side of chip
    tri = [h_ham.x h_ham.y;backside_inside_points(1,:);backside_etch_anchor];
    
    %gtext=gds_structure('label');
    %gtext(end+1)=gdsii_boundarytext(label,[0,0],textsize);
    
    %Larger triangle on bottom left
    p1 = midpt(backside_inside_points(1,:),backside_inside_points(2,:),(backside_inside_points(1,2)-fill_body_contact(2))/(backside_inside_points(1,2)-backside_inside_points(2,2)));
    %p2 = fill_body_contact;
    p2 = fill_body_contact;
    p3 = midpt(backside_inside_points(1,:),backside_inside_points(2,:),(fill_body_contact(1)-backside_inside_points(1,1))/(backside_inside_points(2,1)-backside_inside_points(1,1)));
    ele = gds_element('boundary', 'xy',[p1;p2;p3],'layer',h_latch.layer);
    str_name = sprintf('tri2_[%d,%d]',round(p1(1)),round(p1(2)));
    fill_tri2 = gds_structure(str_name,ele);
    
    %Triangle on the top
    fillamt = .8;
    p1 = midpt(backside_inside_points(4,:),backside_inside_points(1,:),fillamt);
    p2 = backside_inside_points(4,:);
    p3 = midpt(backside_inside_points(4,:),backside_inside_points(3,:),fillamt);
    
    ele = gds_element('boundary', 'xy',[p1;p2;p3],'layer',h_latch.layer);
    str_name = sprintf('tri3_[%d,%d]',round(p1(1)),round(p1(2)));
    fill_tri3 = gds_structure(str_name,ele);
    
    % Add 40 um rectangle on the bottom edge
    dist = 35;
    p1 = backside_inside_points(1,:);
    p2 = backside_inside_points(4,:);
    p3 = midpt(backside_inside_points(4,:),backside_inside_points(3,:),dist/h_rect.w);
    p4 = midpt(backside_inside_points(1,:),backside_inside_points(2,:),dist/h_rect.w);
    ele = gds_element('boundary', 'xy',[p1;p2;p3;p4],'layer',h_latch.layer);
    str_name = sprintf('r_[%d,%d]',round(p1(1)),round(p1(2)));
    fill_rect_2 = gds_structure(str_name,ele);
    
    % Add 30 um rectangle on the left edge
    dist = 40;
    p1 = backside_inside_points(1,:);
    p2 = backside_inside_points(2,:);
    p3 = midpt(backside_inside_points(2,:),backside_inside_points(3,:),dist/h_rect.l);
    p4 = midpt(backside_inside_points(1,:),backside_inside_points(4,:),dist/h_rect.l);
    ele = gds_element('boundary', 'xy',[p1;p2;p3;p4],'layer',h_latch.layer);
    str_name = sprintf('rect2_[%d,%d]',round(p1(1)),round(p1(2)));
    fill_rect = gds_structure(str_name,ele);
    
    %Add fill to right side hammer anchor
    p1 = fill_ham_anchor_t;
    p2 = fill_ham_anchor_b;
    p3 = fill_gap_stop_1;
    p4 = fill_gap_stop_2;
    ele = gds_element('boundary', 'xy',[p1;p2;p3;p4],'layer',h_latch.layer);
    str_name = sprintf('fill_[%d,%d]',round(p1(1)),round(p1(2)));
    fill_rect2 = gds_structure(str_name,ele);
    
    %Add tethers
    snap = 20;
    tether_w = 200;
    pt_t1 = midpt(backside_inside_points(1,:),backside_inside_points(4,:),.5); %temp point
    pt_t2 = midpt(backside_inside_points(2,:),backside_inside_points(3,:),.5); %temp point 2
    pt = midpt(pt_t2,pt_t1,.99); %temp point
    
    h_rect.x = pt(1) - tether_w/2;
    %h_rect.y = pt(2) - trench_w/2;
    h_rect.y = pt(2);
    h_rect.w = tether_w;
    h_rect.theta = tether_ang2;
    h_rect.p0 = pt;
    h_rect.l = trench_w+snap;
    h_rect.layer = 4;
    h_rect.rp = 0;
    tether_strut = rect(h_rect);
    
elseif (h_hammer.chiplet == 2) %Standard not rotated chiplets

    if h_hammer.large_dx
        fill_gap = 30;
    else
        fill_gap = 0;
    end
    
    %The middle of the structure, where the hammer head is
    central_x = fill_ham_head_cent;
    
    %chiplet_w = h_hammer.chiplet_w;
    
    if h_hammer.zif 
        switch h_hammer.zif 
            case 1
                chiplet_w = 3900;
            case 2
                chiplet_w = 4900;
            case 3
                chiplet_w = 4000;
            otherwise
                joey = 1;
        end
   
    else
        chiplet_w = 2*(layout_right - central_x(1) + 25); %Makes chiplet width symmetric around hammer head
    end

    trench_w = 400;
    %fill_y_bottom_lever_arm = latch_point(2) - h_latch.arml - 2*h_latch.r(2) - 25;
    fill_y_bottom_lever_arm = layout_bottom - 25;
    
    
    if h_hammer.compact_latch
        p1 = [central_x(1) - chiplet_w/2 - trench_w/2 central_x(2)+trench_w/2];
        p2 = [central_x(1) + chiplet_w/2 + trench_w/2 central_x(2)+trench_w/2];
        p3 = [central_x(1) + chiplet_w/2 + trench_w/2 fill_y_bottom_lever_arm - trench_w/2];
        p4 = [central_x(1) - chiplet_w/2 - trench_w/2 fill_y_bottom_lever_arm - trench_w/2];
        p5 = p1;
    else
        p1 = [central_x(1) - chiplet_w/2 - trench_w/2 central_x(2)+trench_w/2];
        p2 = [layout_right + 150 central_x(2)+trench_w/2];
        p3 = [layout_right + 150 fill_y_bottom_lever_arm - trench_w/2];
        p4 = [central_x(1) - chiplet_w/2 - trench_w/2 fill_y_bottom_lever_arm - trench_w/2];
        p5 = p1;
    end
    
    be = gds_element('path', 'xy',[p1;p2;p3;p4;p5],'width', trench_w,'layer',5);
    str_name = sprintf('trnch_[%d,%d]',round(p1(1)),round(p1(2)));
    cv2t = gds_structure(str_name,be);
    
            
    %Dummy fill  
    be = gds_element('path', 'xy',[p1;p2;p3;p4;p5],'width', trench_w,'layer',8);
    str_name = sprintf('trn_[%d,%d]%d',round(p1(1)),round(p1(2)),8);
    f_df_3 = gds_structure(str_name,be);
    
    be = gds_element('boundary', 'xy',[p1;p2;p3;p4;p5],'layer',8);
    str_name = sprintf('df_[%d,%d]%d',round(p1(1)),round(p1(2)),8);
    %f_df_3_bound = gds_structure(str_name,be);
    
    %Code to return pegs to proper place
    p1 = [central_x(1) - chiplet_w/2 - trench_w/2 central_x(2)+trench_w/2];
    p2 = [central_x(1) + chiplet_w/2 + trench_w/2 central_x(2)+trench_w/2];
    p3 = [central_x(1) + chiplet_w/2 + trench_w/2 fill_y_bottom_lever_arm - trench_w/2];
    p4 = [central_x(1) - chiplet_w/2 - trench_w/2 fill_y_bottom_lever_arm - trench_w/2];
    p5 = p1;
    
    %Add fill onthe left side
    h_rect.x = central_x(1) - chiplet_w/2;
    h_rect.y = fill_y_bottom_lever_arm;
    h_rect.w = h_ham.x - (central_x(1) - chiplet_w/2) - fill_gap;
    h_rect.l = central_x(2) - fill_y_bottom_lever_arm;
    h_rect.layer = 6;
    %fill_ham2 = rect(h_rect);
        
    %Add fill onthe left side
    h_rect.x = central_x(1) - chiplet_w/2;
    h_rect.y = fill_y_bottom_lever_arm;
    h_rect.w = h_ham.x - (central_x(1) - chiplet_w/2) + h_ham.spring_l;
    h_rect.l = h_ham.y - fill_y_bottom_lever_arm - h_hammer.displacement;
    h_rect.layer = 6;
    %fill_ham2_2 = rect(h_rect);
    
    %Add fill closest to hammer latch
    h_rect.x = central_x(1) - chiplet_w/2;
    h_rect.y = fill_y_bottom_lever_arm;
    h_rect.w = central_x(1) - (central_x(1) - chiplet_w/2) - 40;
    h_rect.l = latch_point(2) - 250 - fill_y_bottom_lever_arm - h_hammer.displacement;
    h_rect.layer = 6;
    %fill_ham2_3 = rect(h_rect);
    
    %Add fll on the right side
    fill_down_length = 550;
    h_rect.x = fill_gap + fill_ham_anchor_b(1);
    h_rect.y = fill_ham_anchor_b(2) - fill_down_length;
    h_rect.w = (central_x(1) + chiplet_w/2)-fill_ham_anchor_b(1)-fill_gap;
    h_rect.l = fill_down_length + (central_x(2)-fill_ham_anchor_b(2));
    h_rect.layer = 6;
    %fill_ham4 = rect(h_rect);
    
    %Add tethers
    tether_w = 200;
    pt = midpt(p3,p4,.5);
    h_rect.x = pt(1) - tether_w/2;
    h_rect.y = pt(2) - trench_w/2;
    h_rect.w = tether_w;
    h_rect.l = trench_w;
    h_rect.layer = 4;
    tether_strut = rect(h_rect);

    if h_hammer.zif==1 % Will not make pegs for zif==3 (jumper)
        peg_w = 350;
        peg_l = 300;
        
        %Left Peg
        h_rect.x = p1(1) + trench_w/2;
        h_rect.y = p1(2) - trench_w/2;
        h_rect.w = peg_w;
        h_rect.l = peg_l;
        h_rect.layer = 4;
        left_peg = rect(h_rect);
        h_rect.layer = 6;
        left_peg_SOI = rect(h_rect);

        %Add cutout for hammer head to slot into
        slot_w = 110;
        slot_h = 43;
        h_rect.x = p1(1) + trench_w/2 + peg_w/2 - slot_w/2;
        h_rect.y = p1(2) - trench_w/2 - 1.5;
        h_rect.w = slot_w;
        h_rect.l = slot_h;
        h_rect.layer = 2;
        left_peg_slot = rect(h_rect);
        

        %Right Peg
        h_rect.x = p2(1) - trench_w/2 - peg_w;
        h_rect.y = p2(2) - trench_w/2;
        h_rect.w = peg_w;
        h_rect.l = peg_l;
        h_rect.layer = 4;
        right_peg = rect(h_rect);
        h_rect.layer = 6;
        right_peg_SOI = rect(h_rect);

        %Add cutout for hammer head to slot into
        h_rect.x = p2(1) - trench_w/2 - peg_w/2 - slot_w/2;
        h_rect.y = p2(2) - trench_w/2 - 1.5;
        h_rect.w = slot_w;
        h_rect.l = slot_h;
        h_rect.layer = 2;
        right_peg_slot = rect(h_rect);
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
    elseif(strcmp(a(i).class,'cell'))
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

out = b;

if h_hammer.rp
    out = [];
    out.body_tl = fill_ham_head_cent - [h_ham.ham_w/2 0];
    out.body_w = h_ham.ham_w;
end