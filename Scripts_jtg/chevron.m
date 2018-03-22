function out=chevron(h_chev)
% Function to crease a thermal Chevron actuator


h_chev.anchor = 100;
h_chev.gap = 40;

h_chev.center_w = 30;
h_chev.beam_w = 10;



%Create anchors
h_rect.x = h_chev.p0(1);
h_rect.y = h_chev.p0(2);
h_rect.w = h_chev.anchor;
h_rect.l = 2*h_chev.anchor + h_chev.beam_num*h_chev.beam_w + (h_chev.beam_num-1)*h_chev.gap;
h_chev.center_l = h_chev.beam_num*h_chev.beam_w + (h_chev.beam_num-1)*h_chev.gap + h_chev.center_w;
h_rect.xnum = 2;
h_rect.xspace = 2*h_chev.beam_l*cos(h_chev.theta*pi/180)+h_chev.center_w;
h_rect.layer = 6;
anchors = rect(h_rect);
h_rect.xnum = 1;

gap = 30;
df_x = h_chev.p0(1) - gap;
df_y = h_chev.p0(2) - gap;
df_w = 2*h_chev.beam_l*cos(h_chev.theta*pi/180)+h_chev.center_w + 2*(gap + h_chev.anchor);



%Create beams on left side
slop = 10;
h_rect.x = h_chev.p0(1) + h_chev.anchor - h_chev.beam_w;
h_rect.y = h_chev.p0(2) + h_chev.anchor;
h_rect.w = h_chev.beam_w;
h_rect.l = h_chev.beam_l+slop;
h_rect.ynum = h_chev.beam_num;
h_rect.yspace = -h_chev.beam_l + h_chev.gap;
h_rect.theta = -pi/180*(90-h_chev.theta);
h_rect.p0 = [h_rect.x+h_chev.beam_w h_rect.y];
h_rect.layer = 6;
beams_l = rect(h_rect);
h_rect.theta = 0;

%Create beams on right side
h_rect.x = h_chev.p0(1) + h_chev.anchor + h_rect.xspace;
h_rect.y = h_chev.p0(2) + h_chev.anchor;
h_rect.w = h_chev.beam_w;
h_rect.l = h_chev.beam_l+slop;
h_rect.ynum = h_chev.beam_num;
h_rect.yspace = -h_chev.beam_l + h_chev.gap;
h_rect.theta = pi/180*(90-h_chev.theta);
h_rect.p0 = [h_rect.x h_rect.y];
h_rect.layer = 6;
beams_r = rect(h_rect);
h_rect.theta = 0;
h_rect.ynum = 1;

%Create central thing
h_rect.x = h_chev.p0(1) + h_chev.anchor + h_rect.xspace/2 - h_chev.center_w/2;
h_rect.y = h_chev.p0(2) + h_chev.beam_l*sin(h_chev.theta*pi/180);
h_rect.w = h_chev.center_w;
h_rect.l = 2*h_chev.anchor + h_chev.beam_num*h_chev.beam_w + (h_chev.beam_num-1)*h_chev.gap;
h_chev.center_l = h_chev.beam_num*h_chev.beam_w + (h_chev.beam_num-1)*h_chev.gap + h_chev.center_w;
h_rect.xnum = 1;
h_rect.layer = 6;
central_rect = rect(h_rect);

df_l = gap + h_chev.beam_l*sin(h_chev.theta*pi/180) + 2*h_chev.anchor + h_chev.beam_num*h_chev.beam_w + (h_chev.beam_num-1)*h_chev.gap;

%Dummy fill
room_to_grow = 50;
h_rectdf.x = df_x;
h_rectdf.y = df_y;
h_rectdf.w = df_w;
h_rectdf.l = df_l + room_to_grow;
h_rectdf.layer = 8;
df_chev = rect(h_rectdf);



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


out = b;


