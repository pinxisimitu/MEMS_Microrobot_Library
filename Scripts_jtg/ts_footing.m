function out = ts_footing(h_tsf)
% Function to make a test structure to measure footing with no 

%Default values
if ~isfield(h_tsf,'beam_l')
    h_tsf.beam_l = 50;
end

if ~isfield(h_tsf,'release_min')    % Minimum oxide undercut during release step
    h_tsf.release_min = 4;
end

if ~isfield(h_tsf,'release_max')    % Maximum oxide undercut during release step
    h_tsf.release_max = 10;
end



%Make an array of lines that will measure the HF release irrespective of
%fooing. 

beam_gap = 80;
release_step = .5;
sbeam_w = 3;            %Supporting beam width
sbeam_l = 200;          %Supporting beam length

lp = h_tsf.p0(1);
num_arms = floor((h_tsf.release_max-h_tsf.release_min)/release_step)+1;
for i = 1:num_arms
    %Releasable anchors
    h_rect.x = lp;
    h_rect.y = h_tsf.p0(2);
    h_rect.w = 2*h_tsf.release_min + (i-1)*2*release_step;
    h_rect.l = h_tsf.beam_l;
    h_rect.layer = 6;
    str = sprintf('line_%d = rect(h_rect);',i);
    eval(str);
    lp = h_rect.x + h_rect.w + beam_gap;
    %Attaching beams
    h_rect.x = h_rect.x + h_rect.w/2 - sbeam_w/2;  
    h_rect.l = -1*sbeam_l;
    h_rect.w = sbeam_w;
    str = sprintf('line2_%d = rect(h_rect);',i);
    eval(str);
end



%Create solid block for beams to attach to
block_w = 150;
block_l = (num_arms+1)*(beam_gap+h_tsf.release_max);

h_rect.x = h_tsf.p0(1) - num_arms*(h_tsf.release_max - h_tsf.release_min)/2;
h_rect.y = h_tsf.p0(2) - sbeam_l - block_w;
h_rect.w = block_l;
h_rect.l = block_w;
h_rect.layer = 6;
block = rect(h_rect);

%Create arms at the bottom of the block
lp = h_tsf.p0(1);
for i = 1:num_arms
    %Releasable anchors
    h_rect.x = lp;
    h_rect.y = h_tsf.p0(2) - block_w - 2*sbeam_l;
    h_rect.w = 2*h_tsf.release_min + 2*(i-1)*release_step;
    h_rect.l = -h_tsf.beam_l;
    h_rect.layer = 6;
    lp = h_rect.x + h_rect.w + beam_gap;
    %Attaching beams
    h_rect.x = h_rect.x + h_rect.w/2 - sbeam_w/2;  
    h_rect.l = sbeam_l;
    h_rect.w = sbeam_w;
    str = sprintf('line4_%d = rect(h_rect);',i);
    eval(str);
end

%Add connecting beam to intersect all other beams (on the top)
h_rect.x = h_tsf.p0(1)+2.5;
h_rect.y = h_tsf.p0(2) - sbeam_l/2;
h_rect.w = block_l;
h_rect.l = sbeam_w;
h_rect.layer = 6;
conbeam = rect(h_rect);

%Add square hold for probe tip.
s = 100;
h_rect.x = h_tsf.p0(1)+release_step+block_l;
h_rect.y = h_tsf.p0(2) - sbeam_l/2-s/2;
h_rect.w = s;
h_rect.l = s;
h_rect.layer = 6;
holder = rect(h_rect);

h_etch.regions = cell(1,1);
h_etch.r = 2;
h_etch.undercut = 3;
section.p0 = [h_tsf.p0(1)+release_step+block_l h_tsf.p0(2) - sbeam_l/2-s/2];
section.type = 'rect';
section.w = s;
section.l = s;
h_etch.regions{1,1} = section;
eh_holder = etch_hole(h_etch);

%Add in more not-SOI

h_rect.x = h_tsf.p0(1)+release_step+block_l + 2*h_etch.undercut + 2*h_etch.r;
h_rect.y = h_tsf.p0(2) - sbeam_l/2-s/2 + 2*h_etch.undercut+ 2*h_etch.r;
h_rect.w = s-4*h_etch.undercut - 4*h_etch.r;
h_rect.l = s-4*h_etch.undercut - 4*h_etch.r;
h_rect.layer = 2;
holder_not_SOI = rect(h_rect);



%Add connecting beam to intersect all other beams (on the bottom)
h_rect.x = h_tsf.p0(1) + 2.5;
h_rect.y = h_tsf.p0(2) - 3*sbeam_l/2 - block_w;
h_rect.w = block_l;
h_rect.l = sbeam_w;
h_rect.layer = 6;
conbeam2 = rect(h_rect);

%Add square hold for probe tip.
s = 100;
h_rect.x = h_tsf.p0(1)+release_step+block_l;
h_rect.y = h_tsf.p0(2) - 3*sbeam_l/2-block_w-s/2;
h_rect.w = s;
h_rect.l = s;
h_rect.layer = 6;
holder2 = rect(h_rect);

h_etch.regions = cell(1,1);
h_etch.r = 2;
h_etch.undercut = 3;
section.p0 = [h_tsf.p0(1)+release_step+block_l h_tsf.p0(2) - 3*sbeam_l/2-block_w-s/2];
section.type = 'rect';
section.w = s;
section.l = s;
h_etch.regions{1,1} = section;
eh_holder2 = etch_hole(h_etch);

%Add in more not-SOI

h_rect.x = h_tsf.p0(1)+release_step+block_l + 2*h_etch.undercut + 2*h_etch.r;
h_rect.y = h_tsf.p0(2) - 3*sbeam_l/2 - block_w - s/2 + 2*h_etch.undercut+ 2*h_etch.r;
h_rect.w = s-4*h_etch.undercut - 4*h_etch.r;
h_rect.l = s-4*h_etch.undercut - 4*h_etch.r;
h_rect.layer = 2;
holder_not_SOI2 = rect(h_rect);


%Add text labels
letter_height = 20;
lp = h_tsf.p0 - [letter_height sbeam_l+block_w/2+letter_height/2];
for i = 1:num_arms
    h_label.p0 = lp;
    h_label.text = sprintf('%.1f',h_tsf.release_min+release_step*(i-1));
    h_label.height = letter_height;
    str = sprintf('lab_%d = add_label(h_label);',i);
    eval(str);  
    lp = lp + [2*h_tsf.release_min + 2*i*release_step + beam_gap 0];
end

%Create squares with etch holes in them
lp = h_tsf.p0 - [0 block_w + 2*sbeam_l];
for ii = 1:num_arms

    %p0 = [h_tsf.p0(1) h_tsf.p0(2) - block_w - 2*sbeam_l];
    undercut = h_tsf.release_min+release_step*(ii-1);
    h_etch.r = 2;
    if ii > num_arms/2
        num_holes = 4;
    else
        num_holes = 5;
    end
    rpoints = zeros(num_holes^2,2);
    tw = (num_holes-1)*(2*undercut + 2*h_etch.r)/sqrt(2)+2*h_etch.r;
    p0 = lp - [tw/2 tw+2*h_etch.r+undercut+ (ii-1)*release_step]; 
    xpoints = zeros(num_holes,num_holes);
    ypoints = zeros(num_holes,num_holes);
    for i=1:num_holes
        for j=1:num_holes
            xpoints(i,j) = (i-1)*(2*undercut + 2*h_etch.r)/sqrt(2) + p0(1) + h_etch.r + undercut;
            ypoints(i,j) = (j-1)*(2*undercut + 2*h_etch.r)/sqrt(2) + p0(2) + h_etch.r + undercut;
        end
    end
    
    
    %Add rect with etch holes
    h_rect.x = p0(1) + h_etch.r + undercut - (h_etch.r + undercut);
    h_rect.y = p0(2) + h_etch.r + undercut - (h_etch.r + undercut);
    h_rect.w = tw + 2*undercut;
    h_rect.l = tw + 2*undercut;
    h_rect.layer = 6;
    str = sprintf('b_%d = rect(h_rect);',ii);
    eval(str);
    
    raw_pts = [reshape(xpoints,[num_holes^2 1]) reshape(ypoints,[num_holes^2 1])];
    h_etch.pts = raw_pts;
    h_etch.regions = cell(1,1);
    section.type = 'raw_points';
    h_etch.regions{1,1} = section;
    h_etch.layer = 2;
    str = sprintf('h_%d = etch_hole(h_etch);',ii);
    eval(str);
    h_etch.layer = 11;
    str = sprintf('hh_%d = etch_hole(h_etch);',ii);
    eval(str);
    lp = lp + [2*h_tsf.release_min + 2*ii*release_step + beam_gap 0];
end



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