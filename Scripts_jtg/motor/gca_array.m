function [gca_struct gca_pts] = gca_array(h_gca) % gap stop parameters
% Function to create Gap Closing Actuators

global joints fid VREP_ignore

VREP_vertical_joint_dist = 4;

% h_gca.top = is the actuator above the shuttle (1=yes, assuming shuttle horiz)
%% Default values

if ~isfield(h_gca,'top')
    h_rect.gca = 1;
end

if ~isfield(h_gca,'NOTDUMMY')
    default_layer_properties;
    h_rect.NOTDUMMY = NOTDUMMY;
end

if ~isfield(h_gca,'VREP_group')
    h_gca.VREP_group = 'dummy_group';
end



%%
gca_pts = [];


gstopL=110;             % length of T piece that holds bumps near angled arm
gstopW=10;              % width of T piece that holds bumps near angled arm
gstopBumpW=2;           % width of bumps on T piece
gstopBumpL=5;           % length of bumps on T piecs
gstopSpace=gstopW+gstopBumpW+h_gca.gstopGap;    % Space from top of T to top of anchor
midBarOL=50;            % Extension past end of GCA for support springs to contact midbar
wbar=20;                % Width of central rotor bar
lbar=h_gca.nFingers/2*(2*h_gca.width+h_gca.gap1+h_gca.gap2)+midBarOL+gstopSpace+h_gca.anchorW; % Length of central rotor bar
dim = [wbar, lbar];

% Generate central rotor (moveable part) bar (jtg)
if h_gca.top
    h_rect.l = lbar;
else
    h_rect.l = -lbar;
end

h_rect.VREP_group = h_gca.VREP_group;

h_rect.x = h_gca.pos(1)-wbar/2;
h_rect.y = h_gca.pos(2);
h_rect.w = wbar;
h_rect.p0 = h_gca.rotation_center;
h_rect.theta = h_gca.rotation_theta;
h_rect.layer = h_gca.layers(1);
h_rect.etch = 1;
h_rect.ETCH_R = 4;
h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
h_rect.etch_layer = h_gca.SOI_HOLE;
h_rect.UNDERCUT = h_gca.UNDERCUT;
h_rect.circle_etch = 0;
rotor_bar = rect(h_rect);
h_rect.etch = 0;

% Add backstops for the rotors
bump_stop = 2;
anchorW_2 = 50;

% This can be deleted after all DC's funcs are changed
rot.theta = h_gca.rotation_theta;
rot.p0 = h_gca.rotation_center;

%% Add backstop for motor

h_rect.VREP_group =   [h_gca.VREP_group  '_a'];

if h_gca.ang == 90
    % Add square anchor for backstop
    p0 = h_gca.pos + [-anchorW_2/2 lbar+2*bump_stop];
    
    h_rect.x = p0(1);
    h_rect.y = p0(2);
    h_rect.w = anchorW_2;
    h_rect.l = anchorW_2;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    anch_1 = rect(h_rect);
    
    % Add small bump on the anchor
    bl = 6;     %bump h_gca.length
    p1 = p0 - [-(anchorW_2-bl)/2 bump_stop];
    
    h_rect.x = p1(1);
    h_rect.y = p1(2);
    h_rect.w = bl;
    h_rect.l = bump_stop;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    anch_1_bump = rect(h_rect);
    
else
    % Add square anchor for backstop
    p0 = h_gca.pos + [-anchorW_2/2 -lbar-2*bump_stop-anchorW_2];
    
    h_rect.x = p0(1);
    h_rect.y = p0(2);
    h_rect.w = anchorW_2;
    h_rect.l = anchorW_2;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    anch_1 = rect(h_rect);
    
    % Add small bump on anchor
    bl = 6;     %bump h_gca.length
    p1 = p0 + [(anchorW_2-bl)/2 anchorW_2];
    
    h_rect.x = p1(1);
    h_rect.y = p1(2);
    h_rect.w = bl;
    h_rect.l = bump_stop;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    anch_1_bump = rect(h_rect);
end

%%

% make support springs for GCA rotor
% orientation determines the direction the springs will go
% orientation==0 -> springs extend out to the left
% orientation==1 -> springs extend out to the right
% h_gca.ang determines whether the array extends out upward or downward from
% its origin

%Set rotation for Comb Arrays
h_ca.rotation_theta = h_gca.rotation_theta;
h_ca.rotation_center = h_gca.rotation_center;
h_ca.length = h_gca.length;
h_ca.length_support = h_gca.lengthSupport;
h_ca.layer = h_gca.layers(1);
h_ca.N = h_gca.nFingers;

df_spring_gap = 10;         %Gap to place around all of the springs

% Set respondable masks 

if ~VREP_ignore
    group = [h_gca.VREP_group  '_r'];
    VREP_SetParam(group,3019,2^4+2^5+2^6+2^7,0);
    
    group = [h_gca.VREP_group  '_a'];
    VREP_SetParam(group,3019,2^0+2^1+2^2+2^3,0);
end

if h_gca.ang==90    % Top GCAs (above the shuttle)
    
    % Check the spring orientation, this determines if the support springs
    % for the rotor extend to the right of left
    if h_gca.springOrient == 0
        % Bottom support spring (close to pawl) for top GCA
        posTemp = h_gca.pos+[-dim(1)/2,h_gca.space];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.springL;
        h_rect.l = h_gca.springW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        BSS_t = rect(h_rect);
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = -h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_BSS_t = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp(1,:) + [wbar h_gca.springW/2]];
        
        % Bottom support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.anchorW;
        h_rect.l = h_gca.anchorW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        g1 = h_rect.VREP_group;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        BSS_t_anchor = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp + [0 h_gca.anchorW/2]];
        
        % Top support spring, top GCA
        posTemp=h_gca.pos+[-dim(1)/2,dim(2)-h_gca.space-h_gca.springW];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.springL;
        h_rect.l = h_gca.springW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        g2 = h_rect.VREP_group;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_t = rect(h_rect);
        
        h_rect.rp = 1;
        pts = rect(h_rect);
        h_rect.rp = 0;
        
        if ~VREP_ignore
            % Add joints at the support spring
            g1 = [h_gca.VREP_group  '_a'];
            g2 = [h_gca.VREP_group  '_d'];
            g3 = [h_gca.VREP_group  '_r'];
            name = [h_gca.VREP_group  '_j1'];
            joints = VREP_add_joint(joints,g3,g2,'prismatic',pts(2,:),h_rect.theta+pi/2,[0 VREP_vertical_joint_dist],name);
            
            name = [h_gca.VREP_group  '_j2'];
            
            joints = VREP_add_joint(joints,g2,g1,'prismatic',pts(2,:),h_rect.theta,[0 VREP_vertical_joint_dist],name);
            
            %Add dummy block that isn't respondable
            VREP_Properties;
            dummy_name = [h_gca.VREP_group '_d'];
            dummy_mass = .01;
            fprintf(fid,'%s = simCreatePureShape(0,8,{%g,%g,%g},%g,nil);\n',dummy_name,1e-2,1e-2,1e-2,dummy_mass);
            %fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
            fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
            fprintf(fid,'simSetObjectName(%s,''%s'')\n\n',dummy_name,dummy_name);
            VREP_SetParam(dummy_name,3019,2^8-1,0);
        end
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = -h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_TSS_t = rect(h_rect);
        
        % top support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,dim(2)-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.anchorW;
        h_rect.l = h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_t_anchor = rect(h_rect);
    else
        % Bottom support spring (close to pawl) for top GCA
        posTemp = h_gca.pos+[dim(1)/2,h_gca.space];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.springL;
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        h_rect.l = h_gca.springW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        BSS_t = rect(h_rect);
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_BSS_t = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp(1,:) + [wbar h_gca.springW/2]];
        
        % Bottom support spring anchor
        posTemp = h_gca.pos+[dim(1)/2+h_gca.springL+h_gca.anchorW,h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.anchorW;
        h_rect.l = h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        g1 = h_rect.VREP_group;
        h_rect.layer = h_gca.layers(1);
        BSS_t_anchor = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp + [0 h_gca.anchorW/2]];
        
        % Top support spring, top GCA
        posTemp=h_gca.pos+[dim(1)/2,dim(2)-h_gca.space-h_gca.springW];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        h_rect.p0 = h_gca.rotation_center;
        g2 = h_rect.VREP_group;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_t = rect(h_rect);
        
        h_rect.rp = 1;
        pts = rect(h_rect);
        h_rect.rp = 0;
        
        % Add joints at the support spring
        g1 = [h_gca.VREP_group  '_a'];
        g2 = [h_gca.VREP_group  '_d'];
        g3 = [h_gca.VREP_group  '_r'];
        name = [h_gca.VREP_group  '_j1'];
        joints = VREP_add_joint(joints,g3,g2,'prismatic',pts(2,:),h_rect.theta+pi/2,[0 VREP_vertical_joint_dist],name);
        
        name = [h_gca.VREP_group  '_j2'];

        joints = VREP_add_joint(joints,g2,g1,'prismatic',pts(2,:),h_rect.theta,[0 VREP_vertical_joint_dist],name);
        
        %Add dummy block that isn't respondable
        VREP_Properties;
        dummy_name = [h_gca.VREP_group '_d'];
        dummy_mass = .01;
        fprintf(fid,'%s = simCreatePureShape(0,8,{%g,%g,%g},%g,nil);\n',dummy_name,1e-2,1e-2,1e-2,dummy_mass);
        %fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
        fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
        fprintf(fid,'simSetObjectName(%s,''%s'')\n\n',dummy_name,dummy_name);
        VREP_SetParam(dummy_name,3019,2^8-1,0);
        
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_TSS_t = rect(h_rect);
        
        % top support spring anchor
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL+h_gca.anchorW,dim(2)-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.l = h_gca.anchorW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_t_anchor = rect(h_rect);
    end
    
    % Define gap stop structures

    % Left anchor for gap stops
    posGstop=h_gca.pos+[-dim(1)/2-h_gca.space-h_gca.anchorW,2*h_gca.space+h_gca.springW];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
    h_rect.w = h_gca.anchorW;
    h_rect.l = h_gca.anchorW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    LGS_t_anchor = rect(h_rect);
    
    % Right anchor for gap stops
    posGstop=h_gca.pos+[dim(1)/2+h_gca.space,2*h_gca.space+h_gca.springW];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
    h_rect.w = h_gca.anchorW;
    h_rect.l = h_gca.anchorW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    RGS_t_anchor = rect(h_rect);
    
    % define gapstop arm on GCA rotor
    posGstop=h_gca.pos+[-gstopL/2,2*h_gca.space+h_gca.springW+h_gca.gstopGap+gstopBumpW+h_gca.anchorW];
    
    % Horizontal bar that will attach gap stop bumps to rotor
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = gstopL;
    h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
    h_rect.l = gstopW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    GCA_horiz_bar_t = rect(h_rect);
    
    % Add left bumps onto horizontal bar
    posGstop=h_gca.pos+[-gstopL/2+(1-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,2*h_gca.space+h_gca.springW+h_gca.anchorW+h_gca.gstopGap];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = gstopBumpL;
    h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
    h_rect.l = gstopBumpW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    h_rect.xnum = 3;
    h_rect.xspace = 10;
    left_top_GCA_bumps = rect(h_rect);
    
    h_rect.xnum = 1;
    h_rect.xspace = 10;
    
    % Add right bumps onto horizontal bar
    posGstop= h_gca.pos+[dim(1)/2+(1-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,2*h_gca.space+h_gca.springW+h_gca.anchorW+h_gca.gstopGap];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = gstopBumpL;
    h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
    h_rect.l = gstopBumpW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    h_rect.xnum = 3;
    h_rect.xspace = 10;
    right_top_GCA_bumps = rect(h_rect);
    
    h_rect.xnum = 1;
    h_rect.xspace = 10;
    
    % Generate left comb array
    posComb=h_gca.pos+[-h_gca.length-2*h_gca.lengthSupport-dim(1)/2,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW];
    
    h_ca.pos = posComb;
    h_ca.width = h_gca.width;
    h_ca.length = h_gca.length;
    h_ca.length_support = h_gca.lengthSupport;
    h_ca.gap1 = h_gca.gap1;
    h_ca.gap2 = h_gca.gap2;
    h_ca.side = 0;
    h_ca.VREP_group = h_gca.VREP_group;
    h_ca.N = h_gca.nFingers;
    h_ca.layer =h_gca.layers(1);
    comb_arr1 = make_comb_array(h_ca);
    
    % Generate Right comb array
    posComb=h_gca.pos+[dim(1)/2+h_gca.lengthSupport,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW];
    h_ca.pos = posComb;
    h_ca.side =1;
    comb_arr2=make_comb_array(h_ca);
    
    % Define stator anchors
    posAnchor=h_gca.pos+[-(h_gca.length+2*h_gca.lengthSupport+dim(1)/2+h_gca.anchorW),2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space];
    posAnchor2=h_gca.pos+[dim(1)/2+h_gca.length+2*h_gca.lengthSupport,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space];
    
    h_rect.x = posAnchor(1);
    h_rect.y = posAnchor(2);
    h_rect.w = h_gca.anchorW;
    h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
    h_rect.l = h_gca.anchorL;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    h_rect.xnum = 2;
    h_rect.xspace = posAnchor2(1) - posAnchor(1) - h_gca.anchorW;
    stator_anch = rect(h_rect);
    h_rect.xnum = 1;
    
elseif h_gca.ang==270 % Bottom GCAs (below the shuttle)
    
    % Check the spring orientation, this determines if the support springs
    % for the rotor extend to the right of left
    if h_gca.springOrient == 0
        % Bottom support spring (close to pawl) for bottom GCA
        posTemp=h_gca.pos+[-dim(1)/2,-h_gca.space-h_gca.springW];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.springL;
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        h_rect.l = h_gca.springW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        BSS_b = rect(h_rect);
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = -h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_BSS_b = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp(1,:) + [wbar h_gca.springW/2]];
        
        % Anchor for Bottom support spring
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.l = h_gca.anchorW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        BSS_b_anchor = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp + [0 h_gca.anchorW/2]];
        
        % Top support spring, bottom GCA
        posTemp=h_gca.pos+[-dim(1)/2,-dim(2)+h_gca.space];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.springL;
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        h_rect.l = h_gca.springW;
        g2 = h_rect.VREP_group;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_b = rect(h_rect);
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = -h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_TSS_b = rect(h_rect);
        
        
        % top support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,-dim(2)+h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.l = h_gca.anchorW;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        g1 = h_rect.VREP_group;
        h_rect.layer = h_gca.layers(1);
        TSS_b_anchor = rect(h_rect);
        
        h_rect.rp = 1;
        pts = rect(h_rect);
        h_rect.rp = 0;
        
        if ~VREP_ignore
            % Add joints at the support spring
            g1 = [h_gca.VREP_group  '_a'];
            g2 = [h_gca.VREP_group  '_d'];
            g3 = [h_gca.VREP_group  '_r'];
            name = [h_gca.VREP_group  '_j1'];
            joints = VREP_add_joint(joints,g3,g2,'prismatic',pts(2,:),h_rect.theta+pi/2,[0 VREP_vertical_joint_dist],name);
            
            name = [h_gca.VREP_group  '_j2'];
            
            joints = VREP_add_joint(joints,g2,g1,'prismatic',pts(2,:),h_rect.theta,[0 VREP_vertical_joint_dist],name);
            
            %Add dummy block that isn't respondable
            VREP_Properties;
            dummy_name = [h_gca.VREP_group '_d'];
            dummy_mass = .01;
            fprintf(fid,'%s = simCreatePureShape(0,8,{%g,%g,%g},%g,nil);\n',dummy_name,1e-2,1e-2,1e-2,dummy_mass);
            %fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
            fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
            fprintf(fid,'simSetObjectName(%s,''%s'')\n\n',dummy_name,dummy_name);
            VREP_SetParam(dummy_name,3019,2^8-1,0);
        end
    else
        % Bottom support spring (close to pawl) for bottom GCA
        posTemp=h_gca.pos+[dim(1)/2,-h_gca.space-h_gca.springW];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        BSS_b = rect(h_rect);
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_BSS_b = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp(1,:) + [wbar h_gca.springW/2]];
        
        % Anchor for Bottom support spring
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL+h_gca.anchorW,-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.anchorW;
        h_rect.l = h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        BSS_b_anchor = rect(h_rect);
        
        gca_pts = [gca_pts; posTemp + [0 h_gca.anchorW/2]];
        
        % Top support spring, bottom GCA
        posTemp=h_gca.pos+[dim(1)/2,-dim(2)+h_gca.space];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW;
        g2 = h_rect.VREP_group;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_b = rect(h_rect);
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]];
        h_rect.x = dfxy1(1);
        h_rect.y = dfxy1(2);
        h_rect.w = h_gca.springL;
        h_rect.l = h_gca.springW+2*df_spring_gap;
        h_rect.p0 = h_gca.rotation_center;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_rect.NOTDUMMY;
        df_TSS_b = rect(h_rect);
        
        
        % top support spring anchor
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL+h_gca.anchorW,-dim(2)+h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        h_rect.x = posTemp(1);
        h_rect.y = posTemp(2);
        h_rect.w = -h_gca.anchorW;
        h_rect.l = h_gca.anchorW;
        h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
        h_rect.p0 = h_gca.rotation_center;
        g1 = h_rect.VREP_group;
        h_rect.theta = h_gca.rotation_theta;
        h_rect.layer = h_gca.layers(1);
        TSS_b_anchor = rect(h_rect);
        
        h_rect.rp = 1;
        pts = rect(h_rect);
        h_rect.rp = 0;
        
        % Add joints at the support spring
        g1 = [h_gca.VREP_group  '_a'];
        g2 = [h_gca.VREP_group  '_d'];
        g3 = [h_gca.VREP_group  '_r'];
        name = [h_gca.VREP_group  '_j1'];
        joints = VREP_add_joint(joints,g3,g2,'prismatic',pts(2,:),h_rect.theta+pi/2,[0 VREP_vertical_joint_dist],name);
        
        name = [h_gca.VREP_group  '_j2'];

        joints = VREP_add_joint(joints,g2,g1,'prismatic',pts(2,:),h_rect.theta,[0 VREP_vertical_joint_dist],name);
        
        %Add dummy block that isn't respondable
        VREP_Properties;
        dummy_name = [h_gca.VREP_group '_d'];
        dummy_mass = .01;
        fprintf(fid,'%s = simCreatePureShape(0,8,{%g,%g,%g},%g,nil);\n',dummy_name,1e-2,1e-2,1e-2,dummy_mass);
        %fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
        fprintf(fid,'simSetObjectPosition(%s,-1,{%g,%g,%g});\n',dummy_name,VREP_S*pts(2,1),VREP_S*pts(2,2),VREP_t);
        fprintf(fid,'simSetObjectName(%s,''%s'')\n\n',dummy_name,dummy_name);
        VREP_SetParam(dummy_name,3019,2^8-1,0);
        
    end
    
    % Define gap stop structures
    
    % Left anchor for gap stops
    posGstop=h_gca.pos+[-dim(1)/2-h_gca.space-h_gca.anchorW,-2*h_gca.space-h_gca.springW-h_gca.anchorW];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = h_gca.anchorW;
    h_rect.l = h_gca.anchorW;
    h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    LGS_b_anchor = rect(h_rect);
    
    % Right anchor for gap stops
    posGstop=h_gca.pos+[dim(1)/2+h_gca.space,-2*h_gca.space-h_gca.springW-h_gca.anchorW];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
    h_rect.w = h_gca.anchorW;
    h_rect.l = h_gca.anchorW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    RGS_b_anchor = rect(h_rect);
    
    % Horizontal bar that will attach gap stop bumps to rotor
    posGstop=h_gca.pos+[-gstopL/2,-(2*h_gca.space+h_gca.springW+h_gca.gstopGap+gstopBumpW+h_gca.anchorW+gstopW)];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = gstopL;
    h_rect.l = gstopW;
    h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    GCA_horiz_bar_b = rect(h_rect);
    
    % Add left bumps onto horizontal bar
    posGstop=h_gca.pos+[-gstopL/2+(1-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,-2*h_gca.space-h_gca.springW-h_gca.anchorW-h_gca.gstopGap-gstopBumpW];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = gstopBumpL;
    h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
    h_rect.l = gstopBumpW;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    h_rect.xnum = 3;
    h_rect.xspace = 10;
    left_top_GCA_bumps_b = rect(h_rect);
    
    h_rect.xnum = 1;
    h_rect.xspace = 10;
    
    % Add right bumps onto horizontal bar
    posGstop=h_gca.pos+[dim(1)/2+(1-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,-2*h_gca.space-h_gca.springW-h_gca.anchorW-h_gca.gstopGap-gstopBumpW];
    h_rect.x = posGstop(1);
    h_rect.y = posGstop(2);
    h_rect.w = gstopBumpL;
    h_rect.l = gstopBumpW;
    h_rect.VREP_group =   [h_gca.VREP_group  '_r'];
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    h_rect.xnum = 3;
    h_rect.xspace = 10;
    right_top_GCA_bumps_b = rect(h_rect);
    
    h_rect.xnum = 1;
    h_rect.xspace = 10;
    
    % Generate left comb array
    posComb=h_gca.pos+[-h_gca.length-2*h_gca.lengthSupport-dim(1)/2,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW)];
    h_ca.pos = posComb;
    h_ca.width = -h_gca.width;
    h_ca.length = h_gca.length;
    h_ca.length_support = h_gca.lengthSupport;
    h_ca.gap1 = -h_gca.gap1;
    h_ca.gap2 = -h_gca.gap2;
    h_ca.side = 0;
    h_ca.N = h_gca.nFingers;
    h_ca.layer = h_gca.layers(1);
    h_ca.VREP_group = h_gca.VREP_group;
    comb_arr1 = make_comb_array(h_ca);
    
    % Generate Right comb array
    posComb = h_gca.pos+[dim(1)/2+h_gca.lengthSupport,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW)];
    h_ca.pos = posComb;
    h_ca.side =1;
    comb_arr2=make_comb_array(h_ca);
    
    % Define stator anchors
    posAnchor=h_gca.pos+[-(h_gca.length+2*h_gca.lengthSupport+dim(1)/2+h_gca.anchorW),-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space)];
    posAnchor2=h_gca.pos+[dim(1)/2+h_gca.length+2*h_gca.lengthSupport,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space)];
    
    h_rect.x = posAnchor(1);
    h_rect.y = posAnchor(2);
    h_rect.w = h_gca.anchorW;
    h_rect.VREP_group =   [h_gca.VREP_group  '_a'];
    h_rect.l = -h_gca.anchorL;
    h_rect.p0 = h_gca.rotation_center;
    h_rect.theta = h_gca.rotation_theta;
    h_rect.layer = h_gca.layers(1);
    h_rect.xnum = 2;
    h_rect.xspace = posAnchor2(1) - posAnchor(1) - h_gca.anchorW;
    stator_anch = rect(h_rect);
    h_rect.xnum = 1;
end

%% Grab all the GDS structures and arrays of structures

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
gca_struct = b;

