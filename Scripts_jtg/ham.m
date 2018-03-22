function [out p1] = ham(h_ham)
% Function to create a hammer
% x = x coordinate of circle
% y = y coordinate of circle
% r  = radius of circle
% n = number of points
% layer = GDS layer to place circle into
% p1 is a matrix of coordinates
%    -p1(1,:) = latch point
%    -p1(2,:) = bottom left point of hammer bottom

if ~isfield(h_ham,'ham_w')
    h_ham.ham_w = 400;          %Width across the hammer (where springs are)
end

if ~isfield(h_ham,'zif')
    h_ham.zif = 0;          %Width across the hammer (where springs are)
end

if ~isfield(h_ham,'anchor_w')
    anchor_w = 50;
else
    anchor_w = h_ham.anchor_w;
end

if ~isfield(h_ham,'piezo')
    h_ham.piezo = 0;          %Width across the hammer (where springs are)
end



h_etch.noetch = h_ham.noetch;

%Dummy fill points
dfill_pts = [];

if h_ham.chiplet == 1
    h_ham.ham_w = 200;          %Width across the hammer (where springs are)
    h_ham.bottom_l = 500;
    h_bottom.l = h_ham.bottom_l;
    h_bottom.d_flat = 40;
elseif(h_ham.large_dx==1)
    h_ham.ham_w = 600;          %Width across the hammer (where springs are)
    h_ham.bottom_l = 450;
    h_bottom.l = h_ham.bottom_l;
    h_bottom.d_flat = h_bottom.l/2.0;
else
    h_ham.bottom_l = 450;
    h_bottom.l = h_ham.bottom_l;
    h_bottom.d_flat = h_bottom.l/2.0;
end





%Bottom of hammer
h_ham.bottom_l = 200;
h_ham.bottom_w = h_ham.ham_w/2.0;

if h_ham.large_dx
    h_ham.num_springs = h_ham.num_springs*2;
end

% Anchors that the springs are attached to
h_anchor.x = h_ham.x;
h_anchor.y = h_ham.y;
h_anchor.w = anchor_w;
h_anchor.l = h_ham.num_springs*(h_ham.spring_gap+h_ham.spring_w)+h_ham.spring_gap;
h_anchor.layer = 6;
h_anchor.xnum = 2;
h_anchor.xspace = h_ham.ham_w + 2*h_ham.spring_l;
hammer_anchors = rect(h_anchor);
output_anchor_pt = [h_anchor.x h_anchor.y+h_anchor.l];
output_anchor_right_pt_t = [h_anchor.x+h_ham.ham_w + 2*h_ham.spring_l+2*anchor_w h_anchor.y+h_anchor.l];
output_anchor_right_pt_b = [h_anchor.x+h_ham.ham_w + 2*h_ham.spring_l+2*anchor_w h_anchor.y];

if h_ham.piezo      %Add more anchoring and add no-dummy layer
    additional_w = 500;
    additional_l = 400;
    h_temp.x = h_ham.x-additional_w;
    h_temp.y = h_ham.y - additional_l;
    h_temp.w = anchor_w+additional_w;
    h_temp.l = h_ham.num_springs*(h_ham.spring_gap+h_ham.spring_w)+h_ham.spring_gap +additional_l;
    h_temp.layer = 6;
    h_temp.xnum = 2;
    h_temp.xspace = h_ham.ham_w + 2*h_ham.spring_l;
    hammer_anchors2 = rect(h_temp);
    
    %Add no-dummy layer
    dummy_ring = 20;
    h_temp.x = h_ham.x-additional_w-dummy_ring;
    h_temp.y = h_ham.y-dummy_ring - additional_l;
    h_temp.w = anchor_w+additional_w+2*dummy_ring;
    h_temp.l = h_ham.num_springs*(h_ham.spring_gap+h_ham.spring_w)+h_ham.spring_gap+2*dummy_ring + additional_l;
    h_temp.layer = 8;
    h_temp.xnum = 2;
    h_temp.xspace = h_ham.ham_w + 2*h_ham.spring_l - 2*dummy_ring;
    df_piezo = rect(h_temp);
end

dfill_pts = [h_anchor.x + h_anchor.w h_anchor.y];
dfill_pts(2,:) = [h_anchor.x + h_anchor.w h_anchor.y-2*h_ham.displacement];


if h_ham.large_dx
    % Add etch holes to anchors
    h_etch.regions = cell(1,1);
    h_etch.r = 2;
    h_etch.undercut = 4;
    section.p0 = [h_anchor.x ,h_anchor.y];
    section.type = 'rect';
    section.w = h_anchor.w;
    section.l = h_anchor.l;
    h_etch.regions{1,1} = section;
    
    ham_anchor_holes_l = etch_hole(h_etch);
    
    % Add etch holes to anchors
    h_etch.regions = cell(1,1);
    h_etch.r = 2;
    h_etch.undercut = 4;
    section.p0 = [h_anchor.x+h_anchor.xspace+h_anchor.w,h_anchor.y];
    section.type = 'rect';
    section.w = h_anchor.w;
    section.l = h_anchor.l;
    h_etch.regions{1,1} = section;
    
    ham_anchor_holes_r = etch_hole(h_etch);
end

% Energy storing springs
h_springs.x = h_anchor.x + h_anchor.w;
h_springs.y = h_anchor.y + h_ham.spring_gap;
h_springs.w = h_ham.spring_l;
h_springs.l = h_ham.spring_w;
h_springs.layer = 6;
h_springs.xnum = 2;
h_springs.xspace = h_ham.ham_w;
h_springs.ynum = h_ham.num_springs;
h_springs.yspace = h_ham.spring_gap;
hammer_springs = rect(h_springs);




%Fillets for springs
h_fillet.xextent = 10;
h_fillet.yextent = 4;
h_fillet.d = .5;
h_fillet.layer = 6;
c = 0;
for i=0:h_ham.num_springs-1
    for j=0:1
        %Bottom left fillet
        c = c + 1;
        h_fillet.p0 = [h_springs.x + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+h_ham.spring_w+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p1 = [h_springs.x + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+h_ham.spring_w+h_fillet.yextent+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p2 = [h_springs.x+h_fillet.xextent + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+h_ham.spring_w+i*(h_ham.spring_gap+h_ham.spring_w)];
        ham_fillet{c} = fillet(h_fillet);
        %Top left fillet
        c = c + 1;
        h_fillet.p0 = [h_springs.x + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p1 = [h_springs.x + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y-h_fillet.yextent+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p2 = [h_springs.x+h_fillet.xextent + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+i*(h_ham.spring_gap+h_ham.spring_w)];
        ham_fillet{c} = fillet(h_fillet);
        %Bottom right fillet
        c = c + 1;
        h_fillet.p0 = [h_springs.x + h_ham.spring_l + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+h_ham.spring_w+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p1 = [h_springs.x + h_ham.spring_l + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+h_ham.spring_w+h_fillet.yextent+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p2 = [h_springs.x + h_ham.spring_l - h_fillet.xextent + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+h_ham.spring_w+i*(h_ham.spring_gap+h_ham.spring_w)];
        ham_fillet{c} = fillet(h_fillet);
        
        %Top right fillet
        c = c + 1;
        h_fillet.p0 = [h_springs.x + h_ham.spring_l  + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p1 = [h_springs.x + h_ham.spring_l  + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y-h_fillet.yextent+i*(h_ham.spring_gap+h_ham.spring_w)];
        h_fillet.p2 = [h_springs.x + h_ham.spring_l  - h_fillet.xextent + j*(h_ham.spring_l+h_ham.ham_w),h_springs.y+i*(h_ham.spring_gap+h_ham.spring_w)];
        ham_fillet{c} = fillet(h_fillet);
    end
end

% Hammer main body
h_hbody.x = h_anchor.x + h_anchor.w + h_ham.spring_l;
h_hbody.y = h_anchor.y;
h_hbody.w = h_ham.ham_w;
h_hbody.l = h_anchor.l;
h_hbody.layer = 6;
hammer_body = rect(h_hbody);
outpoint_body = [h_hbody.x h_hbody.y];


%Cut out anchors from hammer
if h_ham.large_dx
    ldx_anchor = 150;
    ldx_gap = 50;
    ldx_w = 5;
    ldx_p0_l = [h_hbody.x h_hbody.y+h_hbody.l/2];
    ldx_p1_l = ldx_p0_l + [ldx_anchor 0];
    ldx_p15_low = ldx_p0_l + [ldx_anchor-ldx_gap/2 0];
    ldx_p15_high = ldx_p15_low + [0 h_hbody.l/2];
    ldx_p2_l = ldx_p1_l + [0 h_hbody.l/2];
    
    temp = gds_element('path', 'xy',[ldx_p0_l;ldx_p1_l;ldx_p2_l],'width',ldx_w,'layer',2);
    str_name = sprintf('anc_[%d,%d]',round(ldx_p0_l(1)),round(ldx_p0_l(2)));
    cut_ldx = gds_structure(str_name,temp);
    
    %Cut out additional space to stop sticking
    temp = gds_element('path', 'xy',[ldx_p15_low;ldx_p15_high],'width',ldx_gap,'layer',2);
    str_name = sprintf('anc_gap[%d,%d]',round(ldx_p15_low(1)),round(ldx_p15_low(2)));
    cut_ldx_gap = gds_structure(str_name,temp);
    
    ldx_p0_r = [h_hbody.x+h_hbody.w h_hbody.y+h_hbody.l/2];
    ldx_p1_r = ldx_p0_r - [ldx_anchor 0];
    
    ldx_r_low = ldx_p0_r - [ldx_anchor-ldx_gap/2 0];
    ldx_r_high = ldx_r_low + [0 h_hbody.l/2];
    
    ldx_p2_r = ldx_p1_r + [0 h_hbody.l/2];
    
    temp = gds_element('path', 'xy',[ldx_p0_r;ldx_p1_r;ldx_p2_r],'width',ldx_w,'layer',2);
    str_name = sprintf('anc2_[%d,%d]',round(ldx_p0_r(1)),round(ldx_p0_r(2)));
    cut_ldx_right = gds_structure(str_name,temp);
    
    %Cut out additional space to stop sticking
    temp = gds_element('path', 'xy',[ldx_r_low;ldx_r_high],'width',ldx_gap,'layer',2);
    str_name = sprintf('anc_gap[%d,%d]',round(ldx_r_low(1)),round(ldx_r_low(2)));
    cut_ldx_r_gap = gds_structure(str_name,temp);
    
end

dfill_pts(3,:) = [h_hbody.x - 2*h_ham.displacement h_hbody.y - 2*h_ham.displacement];

%Etch holes for hammer main body
if h_ham.large_dx
    h_etch.regions = cell(1,1);
    h_etch.r = 2;
    h_etch.undercut = 6;
    section.p0 = [h_hbody.x, h_hbody.y - h_etch.r];
    section.type = 'rect';
    section.w = h_hbody.w;
    section.l = h_hbody.l/2 + h_etch.r;
    h_etch.regions{1,1} = section;
    
    ham_body_holes = etch_hole(h_etch);
else
    h_etch.regions = cell(1,1);
    h_etch.r = 2;
    h_etch.undercut = 6;
    section.p0 = [h_hbody.x, h_hbody.y - h_etch.r];
    section.type = 'rect';
    section.w = h_hbody.w;
    section.l = h_hbody.l + h_etch.r;
    h_etch.regions{1,1} = section;
    
    ham_body_holes = etch_hole(h_etch);
end



if h_ham.large_dx %large deflection hammer
    %Add extension to hammer
    h_rect = [];
    h_rect.x = ldx_p2_l(1) + ldx_w/2;
    h_rect.y = ldx_p2_l(2);
    h_rect.l = h_ham.displacement - h_ham.hhead_r;
    h_rect.w = ldx_p2_r(1) - ldx_p2_l(1) - ldx_w;
    h_rect.layer = 6;
    addition = rect(h_rect);
    
    %Etch holes in extension
    h_etch.regions = cell(1,1);
    h_etch.r = 2;
    h_etch.undercut = 6;
    section.p0 = [h_rect.x , h_rect.y - h_hbody.l/2-h_etch.undercut];
    section.type = 'rect';
    section.w = h_rect.w ;
    section.l = h_hbody.l/2 +  h_rect.l + h_etch.undercut;
    h_etch.regions{1,1} = section;
    
    ham_extra_holes = etch_hole(h_etch);
    
    
    if h_ham.no_head == 0
        % Hammer Head
        h_cir.x = h_hbody.x + h_hbody.w/2.0;
        h_cir.y = h_hbody.y + h_hbody.l + h_rect.l;
        h_cir.r = h_ham.hhead_r;
        h_cir.layer = 6;
        h_cir.n = 50;
        ham_head = circle(h_cir);
        ham_head_center = [h_cir.x h_cir.y - h_rect.l];
        
        %Etch holes for hammer head
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        h_etch.undercut = 6;
        section.p0 = [h_cir.x,h_cir.y];
        section.type = 'semicircle_t';
        section.r =  h_ham.hhead_r;
        h_etch.regions{1,1} = section;
        
        ham_head_holes = etch_hole(h_etch);
    end
else
    
        % Hammer Head
        h_cir.x = h_hbody.x + h_hbody.w/2.0;
        h_cir.y = h_hbody.y + h_hbody.l;
        h_cir.r = h_ham.hhead_r;
        h_cir.layer = 6;
        h_cir.n = 50;
        if h_ham.no_head == 0
            ham_head = circle(h_cir);
        end
        ham_head_center = [h_cir.x h_cir.y];
        
        %Etch holes for hammer head
        h_etch.regions = cell(1,1);
        h_etch.r = 2;
        h_etch.undercut = 6;
        section.p0 = [h_cir.x,h_cir.y];
        section.type = 'semicircle_t';
        section.r =  h_ham.hhead_r;
        h_etch.regions{1,1} = section;
        if h_ham.no_head == 0
            ham_head_holes = etch_hole(h_etch);
        end
end


% Bottom half of hammer. Where latch grabs on

h_bottom.p0 = [h_hbody.x h_hbody.y];

h_bottom.flat = 20;
if h_ham.stages>0
    h_bottom.w1 = h_ham.ham_w/2.0 + h_bottom.flat/2.0;
else
    h_bottom.w1 = h_ham.ham_w;
end
h_bottom.w2 = h_ham.ham_w;


h_bottom.layer = 6;
h_bottom.noetch = h_ham.noetch;

if h_ham.no_bottom == 1 %Only return the points and not the gds structure
    h_bottom.rp = 1;
else
    h_bottom.rp = 0;
end

[ham_bottom_part output_p1] = hammer_bottom(h_bottom);
    


if h_ham.dage
    dfill_pts(4,:) = output_p1(2,:) - [2*h_ham.displacement 2*h_ham.displacement + 400];
else
    dfill_pts(4,:) = output_p1(2,:) - [2*h_ham.displacement 2*h_ham.displacement];
end

dfill_pts(5,:) = dfill_pts(4,:) + [h_ham.ham_w + 4*h_ham.displacement 0];
dfill_pts(6,:) = dfill_pts(3,:) + [h_ham.ham_w + 4*h_ham.displacement 0];
dfill_pts(7,:) = dfill_pts(2,:) + [h_ham.ham_w + 2*h_ham.spring_l 0];

if h_ham.zif
    dfill_pts(8,:) = [dfill_pts(7,1) h_anchor.y+h_anchor.l+20];
    dfill_pts(9,:) = [dfill_pts(1,1) dfill_pts(8,2)];
else
    %dfill_pts(8,:) = [dfill_pts(7,1) h_anchor.y+h_anchor.l+h_ham.hhead_r+2*h_ham.displacement];
    dfill_pts(8,:) = [dfill_pts(7,1) h_anchor.y+h_anchor.l+20];
    dfill_pts(9,:) = [dfill_pts(1,1) dfill_pts(8,2)];
end


if h_ham.no_bottom == 1
    ham_bottom_part = 0;
    dfill_pts(5,:) = [];
    dfill_pts(4,:) = []; 
end

be = gds_element('boundary', 'xy',dfill_pts,'layer',8);
str_name = sprintf('dff_[%d,%d]%d',round(dfill_pts(1,1)),round(dfill_pts(1,2)),8);
f_df_ne = gds_structure(str_name,be);



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

%str_name = sprintf('Ham_[%d,%d]',round(h_ham.x),round(h_ham.x));

% Outputs a cell array of
out = b;
p1 = [output_p1;output_anchor_pt;outpoint_body;output_anchor_right_pt_t;output_anchor_right_pt_b;ham_head_center];


