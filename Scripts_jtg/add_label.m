function out = add_label(h_label)
% Function to create text labels

if ~isfield(h_label,'layer')
    h_label.layer = 2;
end

if ~isfield(h_label,'data_type')
    h_label.data_type = 0;
end

if ~isfield(h_label,'angle')
    h_label.angle = 0;
end

if ~isfield(h_label,'height')
    h_label.height = 20;
end

str_name = sprintf('Label_%d_%d',h_label.p0(1),h_label.p0(2));

out = gds_structure(str_name,gdsii_boundarytext(h_label.text,h_label.p0,h_label.height,h_label.angle,h_label.layer,h_label.data_type));
