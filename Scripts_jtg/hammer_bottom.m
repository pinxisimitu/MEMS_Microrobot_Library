function [out p_out] = hammer_bottom(h_bottom)
% Function to create the bottom part of a hammer
% p0 = top left coordinate of hammer bottom
% l = length of hammer bottom piece
% w1 = width of hammer bottom piece (at bottom)
% w2 = width of hammer bottom piece (at top)
% flat = length of flat piece
% d_flat = distance to flat piece from bottom of bottom
% layer = GDS layer to place circle into
% downsample = the number by which to downsample the fillets

if ~isfield(h_bottom,'downsample')
    h_bottom.downsample = 1;
end

h_etch.noetch = h_bottom.noetch;

%Find various points that will make up the bottom piece

p0 = h_bottom.p0;
p1 = [p0(1),p0(2)-h_bottom.l];
p2 = [p1(1)+h_bottom.w1,p1(2)];
p3 = [p2(1),p2(2)+h_bottom.d_flat];
p4 = [p0(1)+h_bottom.w2,p0(2)];
p5 = [p0(1)+h_bottom.w1-h_bottom.flat,p0(2)];
p6 = [p3(1)-h_bottom.flat,p3(2)];
p65 = [p6(1)+h_bottom.flat/6,p6(2)];
p8 = [p6(1)+h_bottom.flat/3,p6(2)];
p85 = [p6(1)+h_bottom.flat/2,p6(2)];
p9 = [p6(1)+2*h_bottom.flat/3,p6(2)];
p_out = [p3;p1];


 
%Find central point for spline fit
p_temp = midpt(p4,p6,.8);
p7 = midpt(p5,p_temp,.5);

init_points = [p3' p9' p85' p8' p65' p6' p7' p4'];


interm_pts = fnplt(cscvn(init_points));

%Downsample the number of points.
%final_pts = [interm_pts(1,1) downsample(interm_pts(1,2:end-1),h_fillet.downsample) interm_pts(1,end);
%             interm_pts(2,1) downsample(interm_pts(2,2:end-1),h_fillet.downsample) interm_pts(2,end)];


%Does this work??
final_pts = [interm_pts(:,1),downsample(interm_pts(:,2:end-1),h_bottom.downsample),interm_pts(:,end)];

points = [final_pts';p0;p1;p2];
ca = gds_element('boundary', 'xy',points,'layer',h_bottom.layer);


str_name = sprintf('Hammer_bottom_[%d,%d]',round(p0(1)),round(p0(2)));
bottom_out = gds_structure(str_name,ca);

%Create the sections for the etch holes
h_etch.regions = cell(1,2);

%Same parameters for both sections of etch holes
h_etch.p0 = p1;
%h_etch.w = h_bottom.w1;
h_etch.l = h_bottom.d_flat;
h_etch.l = h_bottom.l;
h_etch.r = 2;
h_etch.undercut = 6;

%Code to ad etch holes to the structure 
%Etch Region1 (bottom rectangle)
section.p0 = p1;
section.type = 'rect';
section.w = h_bottom.w1;
section.l = h_bottom.d_flat;

h_etch.regions{1,1} = section;

%Etch Region2 ()
section2.p0 = [p1(1),p1(2)+h_bottom.d_flat - h_etch.undercut];
section2.type = 'rcurve';
section2.w = h_bottom.w1;
section2.l = h_bottom.l - h_bottom.d_flat + 2*h_etch.undercut;
section2.rcurve = final_pts';

h_etch.regions{1,2} = section2;

etch_holes = etch_hole(h_etch);


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

% Outputs a cell array of GDS structures if return points is not enabled
if h_bottom.rp==0
    out = b;        % Returns GDS structures. 
else
    out = 0;        % Do not return GDS structures
end