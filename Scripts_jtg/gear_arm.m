function [out out_pts] = gear_arm(h_arm)
% Function to make a gear arm with etch holes
if ~isfield(h_arm,'noetch')
    h_arm.noetch = 0;
end

if ~isfield(h_arm,'rack_offset')
    h_arm.rack_offset = 159;
end



%Points will hold particulat points of interest for layout
out_pts = [];


% Create angled arc with etch holes
points = zeros(2*(h_arm.n+1),2);
theta_gap = h_arm.opening_theta*pi/180;
theta_init = h_arm.inner_arm_theta*pi/180 + theta_gap/2;

fillet_pts = zeros(2,2);
fillet_p1s = zeros(2,2);
fillet_p2s = zeros(2,2);

for j=1:2
    if j==2
        mult = -1;
        theta_init = h_arm.inner_arm_theta*pi/180 - theta_gap/2;
    else
        mult = 1;
    end
    for i=1:h_arm.n+1
        points(i + (j-1)*(h_arm.n+1),1) = h_arm.p0(1) + (h_arm.r(j)+.05)*cos(theta_init+mult*(i-1)*(2*pi-theta_gap)/h_arm.n);
        points(i + (j-1)*(h_arm.n+1),2) = h_arm.p0(2) + (h_arm.r(j)+.05)*sin(theta_init+mult*(i-1)*(2*pi-theta_gap)/h_arm.n);
        if i==1 && j==1
            fillet_pts(1,1) = points(i + (j-1)*(h_arm.n+1),1);
            fillet_pts(1,2) = points(i + (j-1)*(h_arm.n+1),2);
        elseif i==h_arm.n+1 && j==1
            fillet_pts(2,1) = points(i + (j-1)*(h_arm.n+1),1);
            fillet_pts(2,2) = points(i + (j-1)*(h_arm.n+1),2);
        elseif i==2 && j==1
            fillet_p1s(1,1) = points(i + (j-1)*(h_arm.n+1),1);
            fillet_p1s(1,2) = points(i + (j-1)*(h_arm.n+1),2);
        elseif i==h_arm.n+1 && j==2
            fillet_p2s(1,1) = points(i + (j-1)*(h_arm.n+1),1);
            fillet_p2s(1,2) = points(i + (j-1)*(h_arm.n+1),2);
        elseif i==h_arm.n && j==1
            fillet_p1s(2,1) = points(i + (j-1)*(h_arm.n+1),1);
            fillet_p1s(2,2) = points(i + (j-1)*(h_arm.n+1),2);
        elseif i==1 && j==2
            fillet_p2s(2,1) = points(i + (j-1)*(h_arm.n+1),1);
            fillet_p2s(2,2) = points(i + (j-1)*(h_arm.n+1),2);
        end
    end
end
rotor_points = [points(end,:);points];
outer_rotor_ring = rotor_points(h_arm.n+1:end,:);

be = gds_element('boundary', 'xy',rotor_points,'layer',h_arm.layer);
str_name = sprintf('Joint_C_[%d,%d],[%d]',round(h_arm.p0(1)),round(h_arm.p0(2)),round(h_arm.r(1)));
c_joint = gds_structure(str_name,be);

%Add fillets at the fillet points
h_fillet.d = .8;
h_fillet.layer = 6;
h_fillet.radial_dist = pi/16;
h_fillet.axial_dist = 10;
h_fillet.p0 =  fillet_pts(1,:);

% Find the correct p1
slope = (fillet_p1s(1,2)- h_fillet.p0(2))/(fillet_p1s(1,1)- h_fillet.p0(1));
ext = abs(25/slope);
%if slope>0
    h_fillet.p1 = fillet_p1s(1,:) - [ext ext*slope];
% else
%     h_fillet.p1 = fillet_p1s(1,:) - [ext ext*slope];
% end

h_fillet.p2 = midpt(fillet_p2s(1,:),h_fillet.p0,.97);

% %Rotate P0,P1, and P2
% h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
% h_rot.theta = 0;
% h_rot.p0 = h_fillet.p0;
% out=rotate_pts(h_rot);
% 
% h_fillet.p0 = out(1,:);
% h_fillet.p1 = out(2,:);
% h_fillet.p2 = out(3,:);

first_fillet_rotor = fillet(h_fillet);

%Add etch hole to fillet
pt = midpt(midpt(h_fillet.p2,h_fillet.p1,.25),h_fillet.p0,.5);

h_cir.r = 2;
h_cir.n = 15;
h_cir.layer = 2;
h_cir.x = pt(1);
h_cir.y = pt(2);
fillet_EH_first = circle(h_cir);


% Add other fillet
%Add fillets at the fillet points
h_fillet.d = .8;
h_fillet.layer = 6;
h_fillet.radial_dist = pi/16;
h_fillet.axial_dist = 10;
h_fillet.p0 =  fillet_pts(2,:);

% Find the correct p1
slope = (fillet_p1s(2,2)- h_fillet.p0(2))/(fillet_p1s(2,1)- h_fillet.p0(1));
ext = abs(25/slope);
%if slope>0
    h_fillet.p1 = fillet_p1s(2,:) + [ext ext*slope];
% else
%     h_fillet.p1 = fillet_p1s(1,:) - [ext ext*slope];
% end

h_fillet.p2 = midpt(fillet_p2s(2,:),h_fillet.p0,.97);

% %Rotate P0,P1, and P2
% h_rot.pts = [h_fillet.p0;h_fillet.p1;h_fillet.p2];
% h_rot.theta = 0;
% h_rot.p0 = h_fillet.p0;
% out=rotate_pts(h_rot);
% 
% h_fillet.p0 = out(1,:);
% h_fillet.p1 = out(2,:);
% h_fillet.p2 = out(3,:);

second_fillet_rotor = fillet(h_fillet);

%Add etch hole to fillet
pt = midpt(midpt(h_fillet.p2,h_fillet.p1,.25),h_fillet.p0,.5);

h_cir.r = 2;
h_cir.n = 15;
h_cir.layer = 2;
h_cir.x = pt(1);
h_cir.y = pt(2);
fillet_EH_second = circle(h_cir);






%Add dummy fill
tooth_height = 42 + 10;              % Does not set tooth height, used to create dummy fill
% Create angled arc with etch holes
points = zeros(2*(h_arm.n+1),2);
theta_extra = 1*(pi/180);         %Add this much to both sides of dummy fill
theta_gap2 = h_arm.opening_theta*pi/180;
theta_gap2 = 2*pi - 2*(2*pi - theta_gap2+theta_extra);


%angular_span = 2*pi - 2*(2*pi - theta_gap2);

%theta_init2 = (180)*pi/180 - h_arm.orientation*angular_span/2 + theta_gap2/2;
theta_init2 = 3*pi/2 + theta_gap2/2 - h_arm.orientation*pi/2;

for j=1:2
    if j==2
        mult = -1;
        %theta_init2 = (180)*pi/180 - h_arm.orientation*angular_span/2 - theta_gap2/2;
        theta_init2 = 3*pi/2 - theta_gap2/2 - h_arm.orientation*pi/2;
        h_arm.r(j) = h_arm.r(j) + tooth_height;
    else
        mult = 1;
    end
    for i=1:h_arm.n+1
        points(i + (j-1)*(h_arm.n+1),1) = h_arm.p0(1) + (h_arm.r(j)+.05)*cos(theta_init2+mult*(i-1)*(2*pi-theta_gap2)/h_arm.n);
        points(i + (j-1)*(h_arm.n+1),2) = h_arm.p0(2) + (h_arm.r(j)+.05)*sin(theta_init2+mult*(i-1)*(2*pi-theta_gap2)/h_arm.n);
    end
end
rotor_points = [points(end,:);points];
outer_rotor_ring = rotor_points(h_arm.n+1:end,:);

be = gds_element('boundary', 'xy',rotor_points,'layer',8);
str_name = sprintf('DF_d_[%d,%d],[%d]',round(h_arm.p0(1)),round(h_arm.p0(2)),round(h_arm.r(1)));
c_joint_df = gds_structure(str_name,be);

h_arm.r(2) = h_arm.r(2) - tooth_height;

%Add gears to the arm
DP = .053;                       %Teeth per 1 unit of diameter
tooth_bottom = 1.25/DP;     	%Distance to bottom of tooth from middle
tooth_top = 1/DP;
r1 = h_arm.r(2)+tooth_bottom;   %Inner radius of gear teeth

angle = 180/pi*[theta_gap+theta_init theta_init+2*pi];
%angle = 180/pi*[0 2*pi];
phi = 20*pi/180;


n1 = 2*r1*DP; % number of teeth on gear 1

[X Y EHP] = drawInvolute(DP, n1, phi,h_arm.p0(1),h_arm.p0(2),3.7*pi/500,angle);

if h_arm.noetch == 0
    %Draw etch holes at EHP
    h_cir.r = 2;
    h_cir.n = 15;
    h_cir.layer = 2;
    h_cir.xy = EHP;
    h_cir.compact = 1;
    EHP_arm = circle(h_cir);
end

be = gds_element('boundary', 'xy',[X', Y'],'layer',h_arm.layer);
str_name = sprintf('gearz_[%d,%d]',round(h_arm.p0(1)),round(h_arm.p0(2)));
gear_teeth = gds_structure(str_name,be);

%Draw rack
rack_offset = h_arm.rack_offset;
rack_teeth = h_arm.rack_teeth;
if h_arm.orientation > 0
    [X Y EHP] = jdrawrack(DP, n1, phi,h_arm.p0(1)+h_arm.orientation*r1+tooth_bottom,h_arm.p0(2)-rack_offset,0,h_arm.orientation,rack_teeth);
else
    [X Y EHP] = jdrawrack(DP, n1, phi,h_arm.p0(1)+h_arm.orientation*r1+tooth_top,h_arm.p0(2)-rack_offset,0,h_arm.orientation,rack_teeth);
end


if h_arm.noetch == 0
    %Add etch holes to rack involutes (already done in the drawrack function
    h_cir.r = 2;
    h_cir.n = 15;
    h_cir.layer = 2;
    h_cir.xy = EHP;
    h_cir.compact = 1;
    EHP_rackz = circle(h_cir);
end
%Add shuttle to rack

%Sort the points to make adding shuttle easier
[aa bb] = sort(Y);
Y = Y(bb);
X = X(bb);


X = [X(1)+h_arm.orientation*h_arm.shuttle_w X X(end)+h_arm.orientation*h_arm.shuttle_w];
Y = [Y(1) Y Y(end)];

out_pts = [out_pts; X(2), Y(2)];

% Add etch holes to the bottom rectangle of the shuttle
rack_p0 = [X(2) Y(2)];

closer_to_edge = 0;                     %Pushing etch holes closer to edges.
h_etch.regions = cell(1,1);
h_etch.r = 2;
if h_arm.orientation ==1
    section.p0 = rack_p0 - [0 closer_to_edge];
else
    section.p0 = rack_p0 - [h_arm.shuttle_w closer_to_edge];
end
section.type = 'rect';
section.w = h_arm.shuttle_w;
section.l = max(Y)-min(Y);
h_etch.regions{1,1} = section;
h_etch.circle_etch = 1;


if h_arm.noetch == 0
    rack_ehhh = etch_hole(h_etch);
end

be = gds_element('boundary', 'xy',[X', Y'],'layer',h_arm.layer);
str_name = sprintf('rack_[%d,%d]',round(h_arm.p0(1)),round(h_arm.p0(2)));
rack_gd = gds_structure(str_name,be);

%Add dummy and pull tabs for Rack
pull_tab_len = 800;

if h_arm.orientation < 0
    h_rect.x = X(1)-2;
    h_rect.y = Y(1)-100;
    h_rect.w = h_arm.shuttle_w + 2 + tooth_height;
    h_rect.l = 2*(max(Y) - min(Y))+300;
    h_rect.layer = 8;
    df_shuttles = rect(h_rect);
    
    if h_arm.manual
        %Lengthen shuttle upwards
        h_rect.x = X(1);
        h_rect.y = Y(end);
        h_rect.w = h_arm.shuttle_w;
        h_rect.l = pull_tab_len;
        h_rect.layer = 6;
        pull_tab = rect(h_rect);
        
        %Add etch holes
        closer_to_edge = 6;                     %Pushing etch holes closer to edges.
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        section.p0 = [h_rect.x h_rect.y - closer_to_edge];
        section.type = 'rect';
        section.w = h_arm.shuttle_w;
        section.l = pull_tab_len;
        h_etch.regions{1,1} = section;
        
        if h_arm.noetch == 0
            rack_ehhh_manl = etch_hole(h_etch);
        end
        %Add probe tip hold at the top
        %Add square hold for probe tip.
        s = 150;
        ppt = [h_rect.x - (s-h_arm.shuttle_w)/2 h_rect.y+h_rect.l-2*closer_to_edge];
        
        h_rect.x = ppt(1);
        h_rect.y = ppt(2);
        h_rect.w = s;
        h_rect.l = s;
        h_rect.layer = 6;
        holder2 = rect(h_rect);
        
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        h_etch.undercut = 6;
        section.p0 = [h_rect.x h_rect.y];
        section.type = 'rect';
        section.w = s;
        section.l = s;
        h_etch.regions{1,1} = section;
        if h_arm.noetch == 0
            eh_holder2 = etch_hole(h_etch);
        end
        %Add in more not-SOI
        h_rect.x = h_rect.x + 2*(h_etch.undercut+h_etch.r);
        h_rect.y = h_rect.y + 2*(h_etch.undercut+h_etch.r);
        h_rect.w = h_rect.w - 4*(h_etch.undercut+h_etch.r);
        h_rect.l = h_rect.l - 4*(h_etch.undercut+h_etch.r);
        h_rect.layer = 2;
        holder_not_SOI2 = rect(h_rect);
        
        %Add serpentine spring to keep probetip pull arm in place
        h_ss.dpp = 400;
        h_ss.p1 = [ppt(1)-50-h_ss.dpp/2 ppt(2)-50];         % First point that the SS will span from
        h_ss.p2 = h_ss.p1 + [0 170];          % Second point that the SS will span to
        h_ss.n = 10;                          % Number of meanders
        h_ss.w = 4;                           % Width of beams
        % Peak to peak distance of meanders
        h_ss.layer = 6;                       % Layer
        serpentine_hhh = s_spring(h_ss);
        
        
        %Add more not dummy for the pully system
        h_rect.x = ppt(1)-100 -  h_ss.dpp;
        h_rect.y = ppt(2)-50;
        h_rect.w = s+150+h_ss.dpp;
        h_rect.l = s+pull_tab_len;
        h_rect.layer = 8;
        holder2_ddff = rect(h_rect);
        
        %Add cross T to attach ss to pull tab
        h_rect.x = h_ss.p2(1) - h_ss.w/2;
        h_rect.y = h_ss.p2(2);
        h_rect.w = h_ss.dpp/2+50+ h_ss.w/2;
        h_rect.l = 30;
        h_rect.layer = 6;
        dpp_t_cross = rect(h_rect);
        
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        h_etch.undercut = 6;
        section.p0 = [h_rect.x h_rect.y];
        section.type = 'rect';
        section.w = h_rect.w+h_etch.undercut;
        section.l = h_rect.l;
        h_etch.regions{1,1} = section;
        if h_arm.noetch == 0
            cross_eh = etch_hole(h_etch);
        end        
    end
else
    h_rect.x = X(1)-h_arm.shuttle_w - tooth_height;
    h_rect.y = Y(1)-100;
    h_rect.w = h_arm.shuttle_w + 2+ tooth_height;
    h_rect.l = 2*(max(Y) - min(Y)) + 300;
    h_rect.layer = 8;
    df_shuttles = rect(h_rect);
    
    
    if h_arm.manual
        %Lengthen shuttle upwards
        h_rect.x = X(1)-h_arm.shuttle_w;
        h_rect.y = Y(end);
        h_rect.w = h_arm.shuttle_w;
        h_rect.l = pull_tab_len;
        h_rect.layer = 6;
        pull_tab = rect(h_rect);
        
        %Add etch holes
        closer_to_edge = 6;                     %Pushing etch holes closer to edges.
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        section.p0 = [h_rect.x h_rect.y - closer_to_edge];
        section.type = 'rect';
        section.w = h_arm.shuttle_w;
        section.l = pull_tab_len;
        h_etch.regions{1,1} = section;
        if h_arm.noetch == 0
            rack_ehhh_manl = etch_hole(h_etch);
        end
        %Add probe tip hold at the top
        %Add square hold for probe tip.
        
        s = 150;
        ppt = [h_rect.x - (s-h_arm.shuttle_w)/2 h_rect.y+h_rect.l-2*closer_to_edge];
        
        h_rect.x = ppt(1);
        h_rect.y = ppt(2);
        h_rect.w = s;
        h_rect.l = s;
        h_rect.layer = 6;
        holder2 = rect(h_rect);
        
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        h_etch.undercut = 6;
        section.p0 = [h_rect.x h_rect.y];
        section.type = 'rect';
        section.w = s;
        section.l = s;
        h_etch.regions{1,1} = section;
        if h_arm.noetch == 0
            eh_holder2 = etch_hole(h_etch);
        end
        %Add in more not-SOI
        
        h_rect.x = h_rect.x + 2*(h_etch.undercut+h_etch.r);
        h_rect.y = h_rect.y + 2*(h_etch.undercut+h_etch.r);
        h_rect.w = h_rect.w - 4*(h_etch.undercut+h_etch.r);
        h_rect.l = h_rect.l - 4*(h_etch.undercut+h_etch.r);
        h_rect.layer = 2;
        holder_not_SOI2 = rect(h_rect);
        
        %Add serpentine spring to keep probetip pull arm in place
        h_ss.dpp = 400;
        h_ss.p1 = [ppt(1)+s+50+h_ss.dpp/2 ppt(2)-50];         % First point that the SS will span from
        h_ss.p2 = h_ss.p1 + [0 170];          % Second point that the SS will span to
        h_ss.n = 10;                          % Number of meanders
        h_ss.w = 4;                           % Width of beams
        % Peak to peak distance of meanders
        h_ss.layer = 6;                       % Layer
        serpentine_hhh = s_spring(h_ss);
        
        
        %Add more not dummy for the pully system
        h_rect.x = ppt(1)-50;
        h_rect.y = ppt(2)-50;
        h_rect.w = s+150+h_ss.dpp;
        h_rect.l = s+pull_tab_len;
        h_rect.layer = 8;
        holder2_ddff = rect(h_rect);
        
        %Add cross T to attach ss to pull tab
        h_rect.x = h_ss.p2(1) -h_ss.dpp/2-50;
        h_rect.y = h_ss.p2(2);
        h_rect.w = h_ss.dpp/2+50+ h_ss.w/2;
        h_rect.l = 30;
        h_rect.layer = 6;
        dpp_t_cross = rect(h_rect);
        
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        h_etch.undercut = 6;
        section.p0 = [h_rect.x h_rect.y];
        section.type = 'rect';
        section.w = h_rect.w+h_etch.undercut;
        section.l = h_rect.l;
        h_etch.regions{1,1} = section;
        if h_arm.noetch == 0
            cross_eh = etch_hole(h_etch);
        end
        
    end
    
end






%Make etch holes in angled part
h_etch.regions = cell(1,1);
h_etch.r = 2;
h_etch.undercut = 6;
section.p0 = [h_arm.p0(1), h_arm.p0(2)];
section.type = 'partial_annulus';
section.r = h_arm.r;
section.theta = [h_arm.inner_arm_theta*pi/180 + theta_gap/2 + (h_etch.undercut+h_etch.r)/h_arm.r(2), 2*pi-theta_gap-(h_etch.undercut+h_etch.r)/h_arm.r(2)];
h_etch.regions{1,1} = section;
if h_arm.noetch == 0
    C_holes = etch_hole(h_etch);
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

% Outputs a cell array of GDS Structures
out = b;