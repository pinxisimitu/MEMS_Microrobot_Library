function out = bond_pad_gripper(h_bpg)
% Boilerplate to make a new function

%Default values
if ~isfield(h_bpg,'rect_size')
    h_bpg.rect_size = 10;
end

%
if ~isfield(h_bpg,'c_gap')
    h_bpg.c_gap = 30;
end

if ~isfield(h_bpg,'cw') % Cantilever width
    h_bpg.cw = 24;
end

if ~isfield(h_bpg,'cl') % Cantilever length
    h_bpg.cl = 300;
end

if ~isfield(h_bpg,'psi') % gripper circle half angle
    h_bpg.psi = 35;
end


if ~isfield(h_bpg,'offset') % gripper circle half angle
    h_bpg.offset = 30;
end

if ~isfield(h_bpg,'hw') % gripper head width
    h_bpg.hw = 50;
end

if ~isfield(h_bpg,'theta') % angle of angled opening at top
    h_bpg.theta = 45;
end



% Add ring around the bond pad
p0 = [h_bpg.p0(1) h_bpg.p0(2)+h_bpg.ring/2];
p1 = p0 + [h_bpg.w-h_bpg.ring/2 0];
p2 = p1 + [0 h_bpg.l-h_bpg.ring];
p3 = p2 - [h_bpg.w-h_bpg.ring 0];
p4 = p3 - [0 h_bpg.l-h_bpg.ring];                          
h_path.pts = [p0;p1;p2;p3;p4];       

h_path.w = h_bpg.ring;                           % Width of path
h_path.layer = 6;                                % Set layer of path
ring_path = m_path(h_path);               % Function to create the path GDS structure

% Add two cantilevers
h_rect.x = h_bpg.p0(1) + h_bpg.w/2 - h_bpg.c_gap/2 - h_bpg.cw;
h_rect.y = h_bpg.p0(2) + h_bpg.ring;
h_rect.w = h_bpg.cw;
h_rect.l = h_bpg.cl;
h_rect.layer = 6;
h_rect.etch = 1;
h_rect.xnum = 2;
h_rect.xspace = h_bpg.c_gap;
cants = rect(h_rect);

left_head_p0 = [h_rect.x+h_bpg.cw h_rect.y+h_rect.l];
right_head_p0 = left_head_p0 + [h_bpg.c_gap 0];

h_rect.etch = 0;

%Add fillets around cantilevers
%left cantilever, left fillet
h_fillet.d = .5;
h_fillet.layer = 6;
h_fillet.p0 = [h_rect.x h_rect.y];
h_fillet.p1 = h_fillet.p0 - [h_bpg.dx 0];
h_fillet.p2 = h_fillet.p0 + [0 2*h_bpg.dx];
fill_ll = fillet(h_fillet);

%left cantilever, right fillet
h_fillet.p0 = [h_rect.x+h_bpg.cw h_rect.y];
h_fillet.p1 = h_fillet.p0 + [h_bpg.dx 0];
h_fillet.p2 = h_fillet.p0 + [0 2*h_bpg.dx];
fill_lr = fillet(h_fillet);

%right cantilever, left fillet
h_fillet.d = .5;
h_fillet.layer = 6;
h_fillet.p0 = [h_rect.x+h_bpg.cw+h_bpg.c_gap h_rect.y];
h_fillet.p1 = h_fillet.p0 - [h_bpg.dx 0];
h_fillet.p2 = h_fillet.p0 + [0 2*h_bpg.dx];
fill_rl = fillet(h_fillet);

%right cantilever, right fillet
h_fillet.p0 = [h_rect.x+h_bpg.cw+h_bpg.cw+h_bpg.c_gap h_rect.y];
h_fillet.p1 = h_fillet.p0 + [h_bpg.dx 0];
h_fillet.p2 = h_fillet.p0 + [0 2*h_bpg.dx];
fill_rr = fillet(h_fillet);




% Make gripper right head
%Create circle arc
angs = pi/180*linspace(-h_bpg.psi,h_bpg.psi,50);
cp = h_bpg.r*[cos(angs)' sin(angs)'];
cp = cp + [right_head_p0(1)-h_bpg.r*cosd(-h_bpg.psi),right_head_p0(2)-h_bpg.r*sind(-h_bpg.psi)];

p0 = cp;
p1 = cp(end,:) + [0 h_bpg.offset];
p2 = p1 + [h_bpg.hw h_bpg.hw*tand(h_bpg.theta)];
p3 = right_head_p0 + [h_bpg.hw 0];
p4 = right_head_p0;

h_shape.pts = [p0;p1;p2;p3;p4];       
h_shape.layer = 6;                     
head_r = m_shape(h_shape);             



% Make gripper left head
%Create circle arc
angs = pi/180*linspace(180+h_bpg.psi,180-h_bpg.psi,50);

cp = h_bpg.r*[cos(angs)' sin(angs)'];
%cp = flipud(cp);
cp = cp + [left_head_p0(1)-h_bpg.r*cosd(180+h_bpg.psi),left_head_p0(2)-h_bpg.r*sind(180+h_bpg.psi)];

p0 = cp;
p1 = cp(end,:) + [0 h_bpg.offset];
p2 = p1 + [-h_bpg.hw h_bpg.hw*tand(h_bpg.theta)];
p3 = left_head_p0 + [-h_bpg.hw 0];
p4 = left_head_p0;

% Define points in the path [x,y]
h_shape.pts = [p0;p1;p2;p3;p4];       
h_shape.layer = 6;                    
head_l = m_shape(h_shape);            

% Add fill around grippers

h_rect.x = h_bpg.p0(1) + h_bpg.ring;
h_rect.y = h_bpg.p0(2) + h_bpg.ring;
h_rect.w = h_bpg.w/2 - h_bpg.ring - h_bpg.cw - h_bpg.c_gap/2 - h_bpg.dx;
h_rect.l = h_bpg.cl - h_bpg.dx;
h_rect.layer = 6;
h_rect.xnum = 2;
h_rect.xspace = h_bpg.c_gap + 2*h_bpg.cw + 2*h_bpg.dx;
bottom_fill = rect(h_rect);

h_rect.x = h_bpg.p0(1) + h_bpg.ring;
h_rect.y = h_bpg.p0(2) + h_bpg.ring + h_bpg.cl - h_bpg.dx;
h_rect.w = h_bpg.w/2 - h_bpg.ring - h_bpg.hw - h_bpg.c_gap/2 - h_bpg.dx;
h_rect.l = h_bpg.l - 2*h_bpg.ring - (h_bpg.cl - h_bpg.dx);
h_rect.layer = 6;
h_rect.xnum = 2;
h_rect.xspace = h_bpg.c_gap + 2*h_bpg.hw + 2*h_bpg.dx;
top_fill = rect(h_rect);



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