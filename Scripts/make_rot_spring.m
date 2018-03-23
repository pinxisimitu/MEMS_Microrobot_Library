function out=make_rot_spring(r_spring)
% pos = position
% N = number of meanders
% w = width of spring
% l = length of each segment
% theta1 = starting angle
% theta2 = ending angle

if ~isfield(r_spring,'pos')
    r_spring.pos = [0,0];
end

if ~isfield(r_spring,'N')
    r_spring.N = 3;
end

if ~isfield(r_spring,'layer')
    default_layer_properties;
    r_spring.layer = SOI;
end

if ~isfield(r_spring,'w')
    r_spring.w = 3;
end

if ~isfield(r_spring,'l')
    r_spring.l = 180;
end

if ~isfield(r_spring,'radius')
    r_spring.radius = 250+r_spring.l/2;
end

if ~isfield(r_spring,'theta1')
    r_spring.theta1 = 0;
end

if ~isfield(r_spring,'theta2')
    r_spring.theta2 = 90;
end

if ~isfield(r_spring,'rp')
    r_spring.rp = 0;
end


%% define points that make up the spring using a circular sine function
theta_range=(r_spring.theta2-r_spring.theta1)*pi/180;
freq=1/r_spring.N*theta_range;
amp = r_spring.l/2;
radius = r_spring.radius+amp;

%number of points per minor peak
step = 200;

%The angle of the major circle
thetaMajor = r_spring.theta1*pi/180:freq/step:r_spring.theta2*pi/180;

%The angle of the superimposed wave
%be careful to allow non-integer number of peaks
tmp = linspace(0,2*pi,step);
idx = mod((1:length(thetaMajor))-1,length(tmp))+1;
thetaMinor = tmp(idx)+pi;

%define the shape
r = radius + amp.*tanh(5*cos(thetaMinor));
x = r_spring.pos(1)+r.*cos(thetaMajor);
y = r_spring.pos(2)+r.*sin(thetaMajor);

% define rotary spring base
xy1(1,:)=[x(1)-10*cosd(r_spring.theta1),y(1)-10*sind(r_spring.theta1)];
xy1(2,:)=xy1(1,:)+[15*cosd(r_spring.theta1),15*sind(r_spring.theta1)];
xy1(3,:)=xy1(2,:)+[10*sind(r_spring.theta1),-10*cosd(r_spring.theta1)];
xy1(4,:)=xy1(1,:)+[10*sind(r_spring.theta1),-10*cosd(r_spring.theta1)];

xy2(1,:)=[x(end)+5*cosd(r_spring.theta2),y(end)+5*sind(r_spring.theta2)];
xy2(2,:)=xy2(1,:)-[15*cosd(r_spring.theta2),15*sind(r_spring.theta2)];
xy2(3,:)=xy2(2,:)-[10*sind(r_spring.theta2),-10*cosd(r_spring.theta2)];
xy2(4,:)=xy2(1,:)-[10*sind(r_spring.theta2),-10*cosd(r_spring.theta2)];



%% create a structure to hold elements
rs1 = gds_element('path', 'xy',[x ; y]','width',r_spring.w,'layer',r_spring.layer);
base1 = gds_element('boundary', 'xy',xy1,'layer',r_spring.layer);
base2 = gds_element('boundary', 'xy',xy2,'layer',r_spring.layer);
gs1 = gds_structure(['rot_spring_' num2str(round(r_spring.pos(1))) '_' num2str(round(r_spring.pos(2)))],rs1);
%gs2 = gds_structure(['rot_spring_base1' num2str(round(r_spring.pos(1))) '_' num2str(round(r_spring.pos(2)))],base1);
%gs3 = gds_structure(['rot_spring_base2' num2str(round(r_spring.pos(1))) '_' num2str(round(r_spring.pos(2)))],base2);


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

out=b;

if r_spring.rp==1
    % Create a top level structure
    ts = gds_structure('TOP');
    ts = add_ref(ts,out);  % add sref elements to top level
    
    % create a library to hold the structure
    glib = gds_library('lib', 'uunit',1e-6, 'dbunit',1e-9, ts,out);

    % Finally write the library to a file, looking first for other files names
    % the same thing and incrementing the appended digit if need be. 
    root='rotor_spring';
    ending = '.gds';
    str = sprintf('%s*%s',root,ending);
    a = dir([str]);
    suffix = [];
    if ~isempty(a)
        k = [];
        for ii=1:length(a)
            k = [k str2num(strtok(strtok(a(ii).name,'.'),root))];
        end
        if isempty(k)
            suffix = 1;
        else
            suffix = max(k)+1;
        end
    end
    % finally write the library to a file
    write_gds_library(glib, [root num2str(suffix) ending]);
end