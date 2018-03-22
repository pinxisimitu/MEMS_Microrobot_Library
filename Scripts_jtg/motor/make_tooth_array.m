function output = make_tooth_array(h_tooth)
% This function makes an array of rounded teeth
% h_tooth.pos: starting h_tooth.position of the teeth. Arrayed right to left
% h_tooth.toothW: width of the teeth
% h_tooth.toothL: length of the teeth out of their base
% h_tooth.toothS: space between the teeth
% h_tooth.Ntooth: number of teeth
% h_tooth.orientation: determines direction of teeth,
%   0=round edge facing up, else point down
% The origin (h_tooth.pos) will be at the rightmost corner of the rightmost tooth
global VREP_ignore

if ~isfield(h_tooth,'theta')
    h_tooth.theta = 0;
end

if ~isfield(h_tooth,'rotation_center')
    h_tooth.rotation_center = [0 0];
end

if ~isfield(h_tooth,'VREP_group')
    h_tooth.VREP_group = 'none';
end

output=gds_structure(['tooth_' num2str(round(h_tooth.pos(1))) '_' num2str(round(h_tooth.pos(2)))]);

for i=1:h_tooth.Ntooth
    if h_tooth.orientation == 0
        xy1=[h_tooth.pos ; h_tooth.pos+[-h_tooth.toothW 0];h_tooth.pos+[-h_tooth.toothW (h_tooth.toothL-h_tooth.toothW/2)];h_tooth.pos+[0 (h_tooth.toothL-h_tooth.toothW/2)];h_tooth.pos];
        
        %Add in rotation to teeth circular centers
        temp_pos = h_tooth.pos+[-h_tooth.toothW/2 h_tooth.toothL-h_tooth.toothW/2];
        rot.pts = temp_pos;
        rot.theta = h_tooth.theta;
        rot.p0 = h_tooth.rotation_center;
        temp_pos=rotate_pts(rot);
        
        tooth=arc(temp_pos,0,h_tooth.toothW/2,h_tooth.theta*180/pi,180+h_tooth.theta*180/pi);
        
        if ~VREP_ignore
            % Add VREP Code here...
            if ~isfield(h_tooth,'VREP_name')
                name = sprintf('T%d_%d%d',round(temp_pos(1,1)),round(temp_pos(1,2)),h_tooth.layer(1));
                name = strrep(name,'-','n');        % Replace any negative signs with 'n'
                vgroup = h_tooth.VREP_group;
            else
                name = h_tooth.VREP_name;
                vgroup = 'none';
            end
            
            cyl_name = sprintf('t_%d_%d');
            VREP_cylinder(temp_pos,h_tooth.toothW/2,name,0,vgroup)
        end
        
    else
        xy1=[h_tooth.pos; h_tooth.pos+[-h_tooth.toothW 0];h_tooth.pos+[-h_tooth.toothW -(h_tooth.toothL-h_tooth.toothW/2)];h_tooth.pos+[0 -(h_tooth.toothL-h_tooth.toothW/2)];h_tooth.pos];
        
        %Add in rotation to teeth circular centers
        temp_pos = h_tooth.pos+[-h_tooth.toothW/2 -(h_tooth.toothL-h_tooth.toothW/2)];
        rot.pts = temp_pos;
        rot.theta = h_tooth.theta;
        rot.p0 = h_tooth.rotation_center;
        temp_pos=rotate_pts(rot);
        
        tooth=arc(temp_pos,0,h_tooth.toothW/2,180+h_tooth.theta*180/pi,360+h_tooth.theta*180/pi);
        
        if ~VREP_ignore
            % Add VREP Code here...
            if ~isfield(h_tooth,'VREP_name')
                name = sprintf('T%d_%d%d',round(temp_pos(1,1)),round(temp_pos(1,2)),h_tooth.layer(1));
                name = strrep(name,'-','n');        % Replace any negative signs with 'n'
                vgroup = h_tooth.VREP_group;
            else
                name = h_tooth.VREP_name;
                vgroup = 'none';
            end
            
            cyl_name = sprintf('t_%d_%d');
            VREP_cylinder(temp_pos,h_tooth.toothW/2,name,0,vgroup)
        end
        
    end
    
    %Add in rotation
    rot.pts = xy1;
    rot.theta = h_tooth.theta;
    rot.p0 = h_tooth.rotation_center;
    xy1=rotate_pts(rot);
    
    output(end+1)=gds_element('boundary', 'xy',xy1,'layer',h_tooth.layer);
    output(end+1)=gdsii_arc(tooth,h_tooth.layer);
    h_tooth.pos = h_tooth.pos-[h_tooth.toothW+h_tooth.toothS,0];
end
