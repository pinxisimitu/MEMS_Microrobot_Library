function out = mech_latch(h_mlatch)
% Function to create the mechanical latch
% p0 = Bottom right point on the latch (where it touches the gap stop)
% gap_stop_width = width of the gap stop
% arm_w = width of lever arm

% theta = angle by which to rotate 
% rot_p0 = point about which to rotate
% rp = return points, if rp==1, returns column vector of verticies

if ~isfield(h_mlatch,'theta')
    h_mlatch.theta = 0;
end
if ~isfield(h_mlatch,'rot_p0')
    h_mlatch.rot_p0 = [0 0];
end
if ~isfield(h_mlatch,'rp')
    h_mlatch.rp = 0;
end
if ~isfield(h_mlatch,'gap')
    h_mlatch.gap = 5;
end

h_mlatch.gap = 5;
h_mlatch.extension = 50; 
h_mlatch.len = h_mlatch.gap + 15;
h_mlatch.hold_w = 30;

%Create points in the latch
p1 = h_mlatch.p0 - [0 h_mlatch.gap+h_mlatch.t]; 
p2 = p1 - h_mlatch.orientation*[h_mlatch.gap_stop_width+(h_mlatch.arm_w+.5)+h_mlatch.t+h_mlatch.extension 0]; 
p3 = p2 + [0 h_mlatch.t];
p4 = p3 + h_mlatch.orientation*[h_mlatch.extension 0];
p5 = p4 + [0 h_mlatch.len]; 
p6 = p5 + h_mlatch.orientation*[h_mlatch.t 0];
p7 = p4 + h_mlatch.orientation*[h_mlatch.t 0];
p8 = p1 - [h_mlatch.orientation*h_mlatch.hold_w -h_mlatch.t];
p9 = p8 + [0 h_mlatch.gap];

points = [h_mlatch.p0;p1;p2;p3;p4;p5;p6;p7;p8;p9];

rot.pts = points;
rot.theta = h_mlatch.theta;
rot.p0 = h_mlatch.rot_p0;
points=rotate_pts(rot);
ca = gds_element('boundary', 'xy',points,'layer',h_mlatch.layer);
str_name = sprintf('Mech_latch_[%d,%d]',round(points(1,1)),round(points(1,2)));
g_mech_latch = gds_structure(str_name,ca);


%Create square pull tab
square_side = 120; 
y_dist = 100; %Distance from end of cantilever to the quare probe holder
rim_w = h_mlatch.t; 

h_rect.x = p3(1) + h_mlatch.orientation*h_mlatch.extension/2 - square_side/2;
h_rect.y = p2(2) - y_dist - square_side;
h_rect.w = square_side;
h_rect.l = square_side;
h_rect.layer = 6;
h_rect.theta = h_mlatch.theta;
h_rect.p0 = h_mlatch.rot_p0;
pull_tab = rect(h_rect);

%Dummy fill for latch
extra = 50;
h_rectdf.x = p3(1) + h_mlatch.orientation*h_mlatch.extension/2 - square_side/2 - extra;
h_rectdf.y = p2(2) - y_dist - square_side - extra;
h_rectdf.w = square_side + 2*extra;
h_rectdf.l = square_side + 2*extra + y_dist;
h_rectdf.layer = 8;
h_rectdf.theta = h_mlatch.theta;
h_rectdf.p0 = h_mlatch.rot_p0;
df_pull_tab = rect(h_rectdf);

%Center of square pull tab
p_cent = [h_rect.x+square_side/2 h_rect.y];
rot.pts = p_cent;
rot.theta = h_mlatch.theta;
rot.p0 = h_mlatch.rot_p0;
p_cent=rotate_pts(rot);

%cutout
h_rect.x = p3(1) +  h_mlatch.orientation*h_mlatch.extension/2 - square_side/2 + rim_w;
h_rect.y = p2(2) - y_dist - square_side + rim_w;
h_rect.w = square_side - 2*rim_w;
h_rect.l = square_side - 2*rim_w;
h_rect.layer = 2;
h_rect.theta = h_mlatch.theta;
h_rect.p0 = h_mlatch.rot_p0;
pull_tab_cut = rect(h_rect);

%Paths to pull tab
path_w = 6;
path_1 = midpt(points(4,:),points(5,:),0);
path_2 = midpt(points(4,:),points(5,:),.4);
path_3 = midpt(points(4,:),points(5,:),.8);


latch_path = gds_element('path', 'xy',[path_1;p_cent],'width', path_w,'layer',h_mlatch.layer);
str_name = sprintf('p_[%d,%d]',path_1(1,1),path_1(1,2));
g_path1 = gds_structure(str_name,latch_path);

latch_path = gds_element('path', 'xy',[path_2;p_cent],'width', path_w,'layer',h_mlatch.layer);
str_name = sprintf('p_[%d,%d]',path_2(1,1),path_2(1,2));
g_path2 = gds_structure(str_name,latch_path);

latch_path = gds_element('path', 'xy',[path_3;p_cent],'width', path_w,'layer',h_mlatch.layer);
str_name = sprintf('p_[%d,%d]',path_3(1,1),path_3(1,2));
g_path3 = gds_structure(str_name,latch_path);

%Cross struts on the pull tab
path_4 = [midpt(path_1,p_cent,.15);midpt(path_3,p_cent,.15)];
path_5 = [midpt(path_1,p_cent,.25);midpt(path_3,p_cent,.25)];
path_6 = [midpt(path_1,p_cent,.35);midpt(path_3,p_cent,.35)];
%path_7 = [midpt(path_1,p_cent,.45);midpt(path_3,p_cent,.45)];


latch_path = gds_element('path', 'xy',path_4,'width', path_w,'layer',h_mlatch.layer);
str_name = sprintf('p_[%d,%d]',path_4(1,1),path_4(1,2));
g_path4 = gds_structure(str_name,latch_path);

latch_path = gds_element('path', 'xy',path_5,'width', path_w,'layer',h_mlatch.layer);
str_name = sprintf('p_[%d,%d]',path_5(1,1),path_5(1,2));
g_path5 = gds_structure(str_name,latch_path);

latch_path = gds_element('path', 'xy',path_6,'width', path_w,'layer',h_mlatch.layer);
str_name = sprintf('p_[%d,%d]',path_6(1,1),path_6(1,2));
g_path6 = gds_structure(str_name,latch_path);

%latch_path = gds_element('path', 'xy',path_7,'width', path_w,'layer',h_mlatch.layer);
%str_name = sprintf('p_[%d,%d]',path_7(1,1),path_7(1,2));
%g_path7 = gds_structure(str_name,latch_path);



%% 
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
            fprintf('Empty Cell! Something went wrong with: %s!!\n',a(i).name)
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


if(h_mlatch.rp == 1)
    out = points;
else
    out = b;
end