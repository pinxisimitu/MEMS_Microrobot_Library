function out = Vapr_relay(h_relay)
% Boilerplate to make a new function

%Default values
if ~isfield(h_relay,'contact_head_type')
    h_relay.contact_head_type = 1;   % 1 square contact
                                     % 2 circle contact
                                     % 3 triangle contact
end

if ~isfield(h_relay,'contact_type')
    h_relay.contact_type = 1;   % 1 Hard contact, 
                                % 2 spring contact
end

if ~isfield(h_relay,'dc_offset')
    h_relay.dc_offset = 0;   % initial DC deflection of comb drive
end

if ~isfield(h_relay,'springw')
    h_relay.springw = 3;   % spring width
end

if ~isfield(h_relay,'drie_undercut')
    h_relay.drie_undercut = .5;   % undercut from the DRIE etch 
end





%How far resonator will move with 1.5 volts applied from ScUM
actuation_distance = 1;
p23_offset = 20;
head_r = 10;        %Radius (or height) of contact head
push1 = 250;


if h_relay.contact_type == 1
    % This is a hard contact (no cantilever)
    %Create path from p0 tp p2
    w1 = 30;
    p0 = h_relay.p0;
    p1 = p0 + [push1 0];
    p2 = p1 + [0 50];
    h_path.pts = [p0;p1;p2];
    
    h_path.w = w1;
    h_path.layer = 6;
    contact_path = m_path(h_path);
elseif h_relay.contact_type == 2
    % Soft contact with cantilever 
    cant_l = 150; %1/2 lenth of final spring
    cant_w = h_relay.springw;
    
    %Create path from p0 tp p2
    w1 = 30;
    p0 = h_relay.p0;
    p1 = p0 + [push1 0];

    h_path.pts = [p0;p1];
    
    h_path.w = w1;
    h_path.layer = 6;
    contact_path = m_path(h_path);
    
    %Make cantilever
    pcant = p1 + [-w1/2 w1/2];
    
    h_rect.y = pcant(2);
    h_rect.x = pcant(1);
    h_rect.etch = 0;
    h_rect.w = cant_w;
    h_rect.l = cant_l;              % Radius of curvature of edges (default is 0)
    cantilever_rect =rect(h_rect);    % Function to create a second rectangle GDS structure
    
    p2 = [p1(1) h_rect.y+cant_l+p23_offset-head_r];
    if h_relay.contact_head_type ==2
        p2 = [p1(1) h_rect.y+cant_l+p23_offset-2*head_r];
    end
end

% Resonator component that contacts contact_path
w2 = 20;


gap1 = 10;
arml = 100;

p3 = p2 - [w1/2+actuation_distance+h_relay.dc_offset+head_r+w2/2-2*h_relay.drie_undercut p23_offset];
p4 = p3 + [0 p23_offset+gap1+w2/2];
p5 = p4 + [arml 0];

% head containing arm contact
h_rect.y = p3(2);                  
h_rect.x = p3(1)-w2/2;                  
h_rect.w = w2;                  
h_rect.l = p4(2)-p3(2); 
h_rect.etch = 1;
h_rect.ETCH_R = 4;
h_rect.etch_layer = 2;
h_rect.UNDERCUT = 6;             % Radius of curvature of edges (default is 0)
p3_rect =rect(h_rect);    % Function to create a second rectangle GDS structure

%Add head to the contact
switch h_relay.contact_head_type
    case {1}    % Rectangular contact
        h_rect.y = p3(2);
        h_rect.x = p3(1)+w2/2;
        h_rect.w = head_r;
        h_rect.l = head_r;
        h_rect.etch = 0;
        contact_head =rect(h_rect);    % Function to create a second rectangle GDS structure
    case {2}    % Circular Contact
        h_cir.r = head_r;
        h_cir.layer = 6;
        h_cir.x = p3(1)+w2/2;
        h_cir.y = p3(2)+head_r;
        contact_head = circle(h_cir);
    case {3}    % Triangular contact
        tri_w = head_r;
        tri_height = head_r;
        
        pp0 = [p3(1)+w2/2 p3(2)];
        pp1 = [p3(1)+w2/2+tri_height p3(2)+tri_w/2];
        pp2 = [p3(1)+w2/2 p3(2)+tri_w];
        
        % Define points in the path [x,y]
        h_shape.pts = [pp0;pp1;pp2];       % Make sure each row is one point [x,y]
        h_shape.layer = 6;                     % Set layer of shape
        contact_head = m_shape(h_shape);             % Function to create the shape GDS structure
end
        




% long horizontal arm contact
h_rect.y = p4(2)-w2/2;                  
h_rect.x = p4(1)-w2/2;                  
h_rect.w = arml+w2;   
h_rect.etch = 1;               
h_rect.l = w2;              % Radius of curvature of edges (default is 0)
p4_rect =rect(h_rect);    % Function to create a second rectangle GDS structure

% Create resonator
res_l = 400;
center_contact = -(p5(2)-p3(2))+head_r;
p6 = p5 - [0 res_l/2-center_contact];
w3 = 134;
sq_h = 150;


%Left vertical rect
h_rect.y = p6(2);                  
h_rect.x = p6(1);                  
h_rect.w = w2;                  
h_rect.l = res_l; 
p6_rect =rect(h_rect);    % Function to create a second rectangle GDS structure

%Central square
p7 = p5 + [w2 -sq_h/2+center_contact];
h_rect.x = p7(1) - 2*h_rect.ETCH_R;                  
h_rect.y = p7(2);                  
h_rect.w = w3 + 4*h_rect.ETCH_R;                  
h_rect.l = sq_h; 
p7_rect =rect(h_rect);    % Function to create a second rectangle GDS structure

%Right Vertial rect
p8 = p6 + [w2+w3 0];
h_rect.y = p8(2);                  
h_rect.x = p8(1);                  
h_rect.w = w2;                  
h_rect.l = res_l; 
p8_rect =rect(h_rect);    % Function to create a second rectangle GDS structure

% Make central anchors
anchor = 100;
%bottom anchor
p9 = p6 + [w2 + (w3-anchor)/2 0];
h_rect.y = p9(2);                  
h_rect.x = p9(1);                  
h_rect.w = anchor;                  
h_rect.l = anchor; 
h_rect.etch = 0;
p9_rect =rect(h_rect);    % Function to create a second rectangle GDS structure

%top anchor
p10 = p9 + [0 res_l-anchor];
h_rect.y = p10(2);                  
h_rect.x = p10(1);                  
h_rect.w = anchor;                  
h_rect.l = anchor; 
h_rect.etch = 0;
p10_rect =rect(h_rect);    % Function to create a second rectangle GDS structure


% Add springs to bottom
springl = 300;
springw = h_relay.springw;


h_rect.x = p6(1)+w2/2-springw/2;                  
h_rect.y = p6(2)-springl;                  
h_rect.w = springw;                  
h_rect.l = springl; 
h_rect.etch = 0;
h_rect.xnum = 4;
h_rect.ynum = 1;
h_rect.xspace = (p8(1)-p6(1))/(3)-springw;
bot_springs =rect(h_rect);    % Function to create a second rectangle GDS structure

% Truss on bottom springs
truss_w = 10;
truss_l = 200;
p11 = [p7(1)+w3/2-truss_l/2 h_rect.y-truss_w];

h_rect.x = p11(1);                  
h_rect.y = p11(2);                  
h_rect.w = truss_l;                  
h_rect.l = truss_w;
h_rect.xnum = 1;
h_rect.ynum = 1;
bot_truss =rect(h_rect);    % Function to create a second rectangle GDS structure


% Add springs to the top
h_rect.x = p6(1)+w2/2-springw/2;                  
h_rect.y = p6(2)+res_l;                  
h_rect.w = springw;                  
h_rect.l = springl; 
h_rect.etch = 0;
h_rect.xnum = 4;
h_rect.ynum = 1;
h_rect.xspace = (p8(1)-p6(1))/(3)-springw;
top_springs =rect(h_rect);    % Function to create a second rectangle GDS structure

% Top Truss
p12 = [p7(1)+w3/2-truss_l/2 h_rect.y+springl];

h_rect.x = p12(1);                  
h_rect.y = p12(2);                  
h_rect.w = truss_l;                  
h_rect.l = truss_w;
h_rect.xnum = 1;
h_rect.ynum = 1;
top_truss =rect(h_rect);    % Function to create a second rectangle GDS structure


% Comb fingers
combw = 2;
combl = 50;
comb_gap = 2;
Ncomb = 50;
left_comb_to_comb = 6;
left_comb_offset = 3;

% Left combs
p13 = p8 + [w2 left_comb_offset];

h_rect.x = p13(1);                  
h_rect.y = p13(2);                  
h_rect.w = combl;                  
h_rect.l = combw; 
h_rect.etch = 0;
h_rect.xnum = 1;
h_rect.ynum = Ncomb;
h_rect.yspace = left_comb_to_comb;
left_combs =rect(h_rect);    % Function to create a second rectangle GDS structure

%Right combs
init_overlap = 30;

p14 = p13 + [2*combl-init_overlap -comb_gap-combw];
h_rect.x = p14(1)-combl;                  
h_rect.y = p14(2);                  
h_rect.w = combl;                  
h_rect.l = combw; 
h_rect.etch = 0;
h_rect.xnum = 1;
h_rect.ynum = Ncomb+1;
h_rect.yspace = left_comb_to_comb;
right_combs =rect(h_rect);    % Function to create a second rectangle GDS structure

%Anchor on the right side
anchor_extension = 5;

h_rect.x = p14(1);                  
h_rect.y = p14(2)-anchor_extension;                  
h_rect.w = anchor;                  
h_rect.l = res_l+2*anchor_extension; 
h_rect.etch = 0;
h_rect.xnum = 1;
h_rect.ynum = 1;
contact_anchor =rect(h_rect);    % Function to create a second rectangle GDS structure

% Add dummy exclude layer
buffer = 30;
h_rect.x = p11(1)-buffer;                  
h_rect.y = p11(2)- buffer;                  
h_rect.w = 2*buffer + truss_l;                  
h_rect.l = 2*buffer + res_l + 2*springl + 2*truss_w; 
h_rect.layer = 8;
h_rect.etch = 0;
h_rect.xnum = 1;
h_rect.ynum = 1;
dummy_relay =rect(h_rect);    % Function to create a second rectangle GDS structure

%Dummy exclude for comb contact
h_rect.x = p14(1) - combl-buffer;                  
h_rect.y = p14(2) - buffer;                  
h_rect.w = 2*buffer + combl + anchor;                  
h_rect.l = 2*buffer + res_l; 
h_rect.layer = 8;
h_rect.etch = 0;
h_rect.xnum = 1;
h_rect.ynum = 1;
dummy_relay_anchor =rect(h_rect);    % Function to create a second rectangle GDS structure

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