function [gca_struct gca_pts] = make_GCA_array(h_gca) % gap stop parameters
% Function to create Gap Closing Actuators 


% h_gca.top = is the actuator above the shuttle (1=yes, assuming shuttle horiz)
%% Default values

if ~isfield(h_gca,'top')
    h_rect.gca = 1;
end


if ~isfield(h_gca,'NOTDUMMY')
    default_layer_properties;
    h_rect.NOTDUMMY = NOTDUMMY;
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

% Make central movable rotor
mid_bar = gds_structure(['mid_bar_' num2str(round(h_gca.pos(1))) '_'  num2str(round(h_gca.pos(2)))]);

% Generate central rotor (moveable part) bar (jtg)
h_rect.x = h_gca.pos(1)-wbar/2;
h_rect.y = h_gca.pos(2);
h_rect.w = wbar;
if h_gca.top
    h_rect.l = lbar;
else
    h_rect.l = -lbar;
end
h_rect.p0 = h_gca.rotation_center;
h_rect.theta = h_gca.rotation_theta;
h_rect.layer = h_gca.layers(1);
h_rect.etch = 1;
h_rect.ETCH_R = h_gca.ETCH_R;
h_rect.etch_layer = h_gca.SOI_HOLE;
h_rect.UNDERCUT = h_gca.UNDERCUT;
rotor_bar = rect(h_rect);

h_rect.etch = 0;
dim = [wbar, lbar];

% Add backstops for the rotors
bump_stop = 2;
anchorW_2 = 50;

%% Add backstop for motor
if h_gca.ang == 90
    p0 = h_gca.pos + [-anchorW_2/2 lbar+2*bump_stop];
    xy1=[p0;
        p0+[0,anchorW_2];
        p0+[anchorW_2,anchorW_2];
        p0+[anchorW_2,0];
        p0];
    
    rot.pts = xy1;
    rot.theta = h_gca.rotation_theta;
    rot.p0 = h_gca.rotation_center;
    xy1=rotate_pts(rot);
    
    mid_bar(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_gca.layers(1));
    
    % add small bump
    bl = 6;     %bump h_gca.length
    p1 = p0 - [-(anchorW_2-bl)/2 bump_stop];
    
    xy1=[p1;
        p1+[0,bump_stop];
        p1+[bl,bump_stop];
        p1+[bl,0];
        p1];
    
    rot.pts = xy1;
    rot.theta = h_gca.rotation_theta;
    rot.p0 = h_gca.rotation_center;
    xy1=rotate_pts(rot);
    
    mid_bar(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_gca.layers(1));
    
else
    p0 = h_gca.pos + [-anchorW_2/2 -lbar-2*bump_stop-anchorW_2];
    xy1=[p0;
        p0+[0,anchorW_2];
        p0+[anchorW_2,anchorW_2];
        p0+[anchorW_2,0];
        p0];
    
    rot.pts = xy1;
    rot.theta = h_gca.rotation_theta;
    rot.p0 = h_gca.rotation_center;
    xy1=rotate_pts(rot);
    
    
    mid_bar(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_gca.layers(1));
    
    
    %add small bump
    bl = 6;     %bump h_gca.length
    p1 = p0 + [(anchorW_2-bl)/2 anchorW_2];
    
    xy1=[p1;
        p1+[0,bump_stop];
        p1+[bl,bump_stop];
        p1+[bl,0];
        p1];
    
    rot.pts = xy1;
    rot.theta = h_gca.rotation_theta;
    rot.p0 = h_gca.rotation_center;
    xy1=rotate_pts(rot);
    
    
    mid_bar(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_gca.layers(1));
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

df_spring_gap = 10;         %Gap to place around allof the springs

if h_gca.springOrient==0
    if h_gca.ang==90
        
        % bottom support spring
        posTemp = h_gca.pos+[-dim(1)/2,h_gca.space];
        xy1 = [posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[-h_gca.springL,h_gca.springW] ; posTemp+[-h_gca.springL,0] ; posTemp];
        
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]; 
                 posTemp + [0,h_gca.springW+df_spring_gap]; 
                 posTemp + [-h_gca.springL,h_gca.springW+df_spring_gap];
                 posTemp + [-h_gca.springL,-df_spring_gap];
                 posTemp - [0 df_spring_gap]];
        
        rot.pts = dfxy1;
        dfxy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',dfxy2,'layer',h_rect.NOTDUMMY);
        
        
        
        gca_pts = [gca_pts; xy1(1,:) + [wbar h_gca.springW/2]];
        
        % bottom support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        gca_pts = [gca_pts; posTemp + [0 h_gca.anchorW/2]];
        
        
        % top support spring
        posTemp=h_gca.pos+[-dim(1)/2,dim(2)-h_gca.space-h_gca.springW];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[-h_gca.springL,h_gca.springW] ; posTemp+[-h_gca.springL,0] ; posTemp];
        
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap]; 
                 posTemp + [0,h_gca.springW+df_spring_gap]; 
                 posTemp + [-h_gca.springL,h_gca.springW+df_spring_gap];
                 posTemp + [-h_gca.springL,-df_spring_gap];
                 posTemp - [0 df_spring_gap]];
        
        rot.pts = dfxy1;
        dfxy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',dfxy2,'layer',h_rect.NOTDUMMY);
        
        
        
        % top support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,dim(2)-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % start defining gap stop structures
        % define anchors that the gapstop will bump against
        posGstop=h_gca.pos+[-dim(1)/2-h_gca.space-h_gca.anchorW,2*h_gca.space+h_gca.springW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posGstop=h_gca.pos+[dim(1)/2+h_gca.space,2*h_gca.space+h_gca.springW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % define gapstop arm on GCA rotor
        posGstop=h_gca.pos+[-gstopL/2,2*h_gca.space+h_gca.springW+h_gca.gstopGap+gstopBumpW+h_gca.anchorW];
        xy1=[posGstop ; posGstop+[0,gstopW] ; posGstop+[gstopL,gstopW] ; posGstop+[gstopL,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % make gapstop bumps on GCA rotor
        for i=1:3
            posGstop=h_gca.pos+[-gstopL/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,2*h_gca.space+h_gca.springW+h_gca.anchorW+h_gca.gstopGap];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        for i=1:3
            posGstop=h_gca.pos+[dim(1)/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,2*h_gca.space+h_gca.springW+h_gca.anchorW+h_gca.gstopGap];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        % make comb array
        % left comb array
        posComb=h_gca.pos+[-h_gca.length-2*h_gca.lengthSupport-dim(1)/2,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW];
        
        %Start here, replace all make_comb_array(...) with make_comb_array(h_ca)
        h_ca.pos = posComb;
        h_ca.width = h_gca.width;
        h_ca.length = h_gca.length;
        h_ca.length_support = h_gca.lengthSupport;
        h_ca.gap1 = h_gca.gap1;
        h_ca.gap2 = h_gca.gap2;
        h_ca.side = 0;
        h_ca.N = h_gca.nFingers;
        h_ca.layer =h_gca.layers(1);
        comb_arr1=make_comb_array(h_ca);
        
        
        mid_bar=join_gds_structures(mid_bar,comb_arr1);
        
        % right comb array
        posComb=h_gca.pos+[dim(1)/2+h_gca.lengthSupport,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW];
        
        
        h_ca.pos = posComb;
        h_ca.side =1;
        comb_arr2=make_comb_array(h_ca);
        mid_bar=join_gds_structures(mid_bar,comb_arr2);
        
        % define stator anchors
        posAnchor=h_gca.pos+[-(h_gca.length+2*h_gca.lengthSupport+dim(1)/2+h_gca.anchorW),2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space];
        xy1=[posAnchor ; posAnchor+[0,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posAnchor=h_gca.pos+[dim(1)/2+h_gca.length+2*h_gca.lengthSupport,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space];
        xy1=[posAnchor ; posAnchor+[0,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
                
    elseif h_gca.ang==270
        % bottom support spring
        posTemp=h_gca.pos+[-dim(1)/2,-h_gca.space-h_gca.springW];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[-h_gca.springL,h_gca.springW] ; posTemp+[-h_gca.springL,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);        
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap];
            posTemp + [0,h_gca.springW+df_spring_gap];
            posTemp + [-h_gca.springL,h_gca.springW+df_spring_gap];
            posTemp + [-h_gca.springL,-df_spring_gap];
            posTemp - [0 df_spring_gap]];
        
        rot.pts = dfxy1;
        dfxy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',dfxy2,'layer',h_rect.NOTDUMMY);
        
        
        
        gca_pts = [gca_pts; xy1(1,:) + [wbar h_gca.springW/2]];
        
        % bottom support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        gca_pts = [gca_pts; posTemp + [0 h_gca.anchorW/2]];
        
        
        % top support spring
        posTemp=h_gca.pos+[-dim(1)/2,-dim(2)+h_gca.space];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[-h_gca.springL,h_gca.springW] ; posTemp+[-h_gca.springL,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        %Dummy fill for support spring
        dfxy1 = [posTemp - [0 df_spring_gap];
            posTemp + [0,h_gca.springW+df_spring_gap];
            posTemp + [-h_gca.springL,h_gca.springW+df_spring_gap];
            posTemp + [-h_gca.springL,-df_spring_gap];
            posTemp - [0 df_spring_gap]];
        
        rot.pts = dfxy1;
        dfxy2=rotate_pts(rot);
        
        mid_bar(end+1)=gds_element('boundary', 'xy',dfxy2,'layer',h_rect.NOTDUMMY);
        
        
        
        
        % top support spring anchor
        posTemp=h_gca.pos+[-dim(1)/2-h_gca.springL-h_gca.anchorW,-dim(2)+h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % start defining gap stop structures
        % define anchors that the gapstop will bump against
        posGstop=h_gca.pos+[-dim(1)/2-h_gca.space-h_gca.anchorW,-2*h_gca.space-h_gca.springW-h_gca.anchorW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % define gapstop arm on GCA rotor
        posGstop=h_gca.pos+[dim(1)/2+h_gca.space,-2*h_gca.space-h_gca.springW-h_gca.anchorW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posGstop=h_gca.pos+[-gstopL/2,-(2*h_gca.space+h_gca.springW+h_gca.gstopGap+gstopBumpW+h_gca.anchorW+gstopW)];
        xy1=[posGstop ; posGstop+[0,gstopW] ; posGstop+[gstopL,gstopW] ; posGstop+[gstopL,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % make gapstop bumps on GCA rotor
        for i=1:3
            posGstop=h_gca.pos+[-gstopL/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,-2*h_gca.space-h_gca.springW-h_gca.anchorW-h_gca.gstopGap-gstopBumpW];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        for i=1:3
            posGstop=h_gca.pos+[dim(1)/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,-2*h_gca.space-h_gca.springW-h_gca.anchorW-h_gca.gstopGap-gstopBumpW];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        % make comb array
        % left comb array
        posComb=h_gca.pos+[-h_gca.length-2*h_gca.lengthSupport-dim(1)/2,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW)];
        
        h_ca.pos = posComb;
        h_ca.width = -h_gca.width;
        h_ca.gap1 = -h_gca.gap1;
        h_ca.gap2 = -h_gca.gap2;
        h_ca.side = 0;
        
        h_ca.length =h_gca.length;
        h_ca.length_support =h_gca.lengthSupport;
        h_ca.N = h_gca.nFingers;
        h_ca.layer =h_gca.layers(1);
        
        
        comb_arr1=make_comb_array(h_ca);
        
        %comb_arr1=make_comb_array(posComb,-h_gca.width,h_gca.length,h_gca.lengthSupport,-h_gca.gap1,-h_gca.gap2,0,h_gca.nFingers,'1',h_gca.layers(1));
        mid_bar=join_gds_structures(mid_bar,comb_arr1);
        
        % right comb array
        posComb = h_gca.pos+[dim(1)/2+h_gca.lengthSupport,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW)];
        
        h_ca.pos = posComb;
        h_ca.side = 1;
        comb_arr2 = make_comb_array(h_ca);
        %comb_arr2 = make_comb_array(posComb,-h_gca.width,h_gca.length,h_gca.lengthSupport,-h_gca.gap1,-h_gca.gap2,1,h_gca.nFingers,'2',h_gca.layers(1));
        mid_bar = join_gds_structures(mid_bar,comb_arr2);
        
        % define stator anchors
        posAnchor=h_gca.pos+[-(h_gca.length+2*h_gca.lengthSupport+dim(1)/2+h_gca.anchorW),-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space)];
        xy1=[posAnchor ; posAnchor+[0,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posAnchor=h_gca.pos+[dim(1)/2+h_gca.length+2*h_gca.lengthSupport,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space)];
        xy1=[posAnchor ; posAnchor+[0,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot); 
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
    end
elseif h_gca.springOrient==1
    if h_gca.ang==90
        % bottom support spring
        posTemp=h_gca.pos+[dim(1)/2,h_gca.space];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[h_gca.springL,h_gca.springW] ; posTemp+[h_gca.springL,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        gca_pts = [gca_pts; posTemp];
        
        % bottom support spring anchor
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL,h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        gca_pts = [gca_pts; posTemp];
        
        % top support spring
        posTemp=h_gca.pos+[dim(1)/2,dim(2)-h_gca.space-h_gca.springW];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[h_gca.springL,h_gca.springW] ; posTemp+[h_gca.springL,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % top support spring anchor
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL,dim(2)-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % start defining gap stop structures
        % define anchors that the gapstop will bump against
        posGstop=h_gca.pos+[-dim(1)/2-h_gca.space-h_gca.anchorW,2*h_gca.space+h_gca.springW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posGstop=h_gca.pos+[dim(1)/2+h_gca.space,2*h_gca.space+h_gca.springW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % define gapstop arm on GCA rotor
        posGstop=h_gca.pos+[-gstopL/2,2*h_gca.space+h_gca.springW+h_gca.gstopGap+gstopBumpW+h_gca.anchorW];
        xy1=[posGstop ; posGstop+[0,gstopW] ; posGstop+[gstopL,gstopW] ; posGstop+[gstopL,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % make gapstop bumps on GCA rotor
        for i=1:3
            posGstop=h_gca.pos+[-gstopL/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,2*h_gca.space+h_gca.springW+h_gca.anchorW+h_gca.gstopGap];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        for i=1:3
            posGstop=h_gca.pos+[dim(1)/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,2*h_gca.space+h_gca.springW+h_gca.anchorW+h_gca.gstopGap];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        % make comb array
        % left comb array
        posComb=h_gca.pos+[-h_gca.length-2*h_gca.lengthSupport-dim(1)/2,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW];
        
        h_ca.pos = posComb;
        h_ca.width = h_gca.width;
        h_ca.gap1 = h_gca.gap1;
        h_ca.gap2 = h_gca.gap2;
        h_ca.side = 0;
        h_ca.length =h_gca.length;
        
        h_ca.length =h_gca.length;
        h_ca.length_support =h_gca.lengthSupport;
        h_ca.N = h_gca.nFingers;
        h_ca.layer = h_gca.layers(1);
        
        
        comb_arr1=make_comb_array(h_ca);
        
        %comb_arr1=make_comb_array(posComb,h_gca.width,h_gca.length,h_gca.lengthSupport,h_gca.gap1,h_gca.gap2,0,h_gca.nFingers,'1',h_gca.layers(1));
        mid_bar=join_gds_structures(mid_bar,comb_arr1);
        
        % right comb array
        posComb=h_gca.pos+[dim(1)/2+h_gca.lengthSupport,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW];
        h_ca.pos = posComb;
        h_ca.side = 1;
        comb_arr2=make_comb_array(h_ca);
        
        %comb_arr2=make_comb_array(posComb,h_gca.width,h_gca.length,h_gca.lengthSupport,h_gca.gap1,h_gca.gap2,1,h_gca.nFingers,'2',h_gca.layers(1));
        mid_bar=join_gds_structures(mid_bar,comb_arr2);
        
        % define stator anchors
        posAnchor=h_gca.pos+[-(h_gca.length+2*h_gca.lengthSupport+dim(1)/2+h_gca.anchorW),2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space];
        xy1=[posAnchor ; posAnchor+[0,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posAnchor=h_gca.pos+[dim(1)/2+h_gca.length+2*h_gca.lengthSupport,2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space];
        xy1=[posAnchor ; posAnchor+[0,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
    elseif h_gca.ang == 270
        % bottom support spring
        posTemp=h_gca.pos+[dim(1)/2,-h_gca.space-h_gca.springW];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[h_gca.springL,h_gca.springW] ; posTemp+[h_gca.springL,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        gca_pts = [gca_pts; posTemp];
        
        % bottom support spring anchor
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL,-h_gca.space-h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        gca_pts = [gca_pts; posTemp];
        
        % top support spring
        posTemp=h_gca.pos+[dim(1)/2,-dim(2)+h_gca.space];
        xy1=[posTemp ; posTemp+[0,h_gca.springW] ; posTemp+[h_gca.springL,h_gca.springW] ; posTemp+[h_gca.springL,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % top support spring anchor
        posTemp=h_gca.pos+[dim(1)/2+h_gca.springL,-dim(2)+h_gca.space+h_gca.springW/2-h_gca.anchorW/2];
        xy1=[posTemp ; posTemp+[0,h_gca.anchorW] ; posTemp+[h_gca.anchorW,h_gca.anchorW] ; posTemp+[h_gca.anchorW,0] ; posTemp];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % start defining gap stop structures
        % define anchors that the gapstop will bump against
        posGstop=h_gca.pos+[-dim(1)/2-h_gca.space-h_gca.anchorW,-2*h_gca.space-h_gca.springW-h_gca.anchorW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posGstop=h_gca.pos+[dim(1)/2+h_gca.space,-2*h_gca.space-h_gca.springW-h_gca.anchorW];
        xy1=[posGstop ; posGstop+[0,h_gca.anchorW] ; posGstop+[h_gca.anchorW,h_gca.anchorW] ; posGstop+[h_gca.anchorW,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % define gapstop arm on GCA rotor
        posGstop=h_gca.pos+[-gstopL/2,-(2*h_gca.space+h_gca.springW+h_gca.gstopGap+gstopBumpW+h_gca.anchorW+gstopW)];
        xy1=[posGstop ; posGstop+[0,gstopW] ; posGstop+[gstopL,gstopW] ; posGstop+[gstopL,0] ; posGstop];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        % make gapstop bumps on GCA rotor
        for i=1:3
            posGstop=h_gca.pos+[-gstopL/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,-2*h_gca.space-h_gca.springW-h_gca.anchorW-h_gca.gstopGap-gstopBumpW];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        for i=1:3
            posGstop=h_gca.pos+[dim(1)/2+(i-1/2)*(gstopL/2-dim(1)/2)/3-gstopBumpL/2,-2*h_gca.space-h_gca.springW-h_gca.anchorW-h_gca.gstopGap-gstopBumpW];
            xy1=[posGstop ; posGstop+[0,gstopBumpW] ; posGstop+[gstopBumpL,gstopBumpW] ; posGstop+[gstopBumpL,0] ; posGstop];
            rot.pts = xy1;
            xy2=rotate_pts(rot);
            mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        end
        
        
        % make comb array
        % left comb array
        posComb=h_gca.pos+[-h_gca.length-2*h_gca.lengthSupport-dim(1)/2,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW)];
        h_ca.pos = posComb;
        h_ca.width = -h_gca.width;
        h_ca.gap1 = -h_gca.gap1;
        h_ca.gap2 = -h_gca.gap2;
        h_ca.side = 0;
        comb_arr1=make_comb_array(h_ca);
        
        
        %comb_arr1=make_comb_array(posComb,-h_gca.width,h_gca.length,h_gca.lengthSupport,-h_gca.gap1,-h_gca.gap2,0,h_gca.nFingers,'1',h_gca.layers(1));
        mid_bar=join_gds_structures(mid_bar,comb_arr1);
        
        % right comb array
        posComb=h_gca.pos+[dim(1)/2+h_gca.lengthSupport,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW)];
        h_ca.pos = posComb;
        h_ca.width = -h_gca.width;
        h_ca.gap1 = -h_gca.gap1;
        h_ca.gap2 = -h_gca.gap2;
        h_ca.side = 1;
        comb_arr2=make_comb_array(h_ca);
        
        %comb_arr2=make_comb_array(posComb,-h_gca.width,h_gca.length,h_gca.lengthSupport,-h_gca.gap1,-h_gca.gap2,1,h_gca.nFingers,'2',h_gca.layers(1));
        mid_bar=join_gds_structures(mid_bar,comb_arr2);
        
        % define stator anchors
        posAnchor=h_gca.pos+[-(h_gca.length+2*h_gca.lengthSupport+dim(1)/2+h_gca.anchorW),-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space)];
        xy1=[posAnchor ; posAnchor+[0,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
        
        posAnchor=h_gca.pos+[dim(1)/2+h_gca.length+2*h_gca.lengthSupport,-(2*h_gca.space+h_gca.gap2+h_gca.springW+h_gca.anchorW+h_gca.gstopGap+gstopBumpW+gstopW-h_gca.space)];
        xy1=[posAnchor ; posAnchor+[0,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,-h_gca.anchorL] ; posAnchor+[h_gca.anchorW,0] ; posAnchor];
        rot.pts = xy1;
        xy2=rotate_pts(rot);
        mid_bar(end+1)=gds_element('boundary', 'xy',xy2,'layer',h_gca.layers(1));
    end
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

