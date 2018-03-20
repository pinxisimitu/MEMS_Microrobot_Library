% define a bar with width and length as the dimensions
% width and length will be approximate according to the number of holes
% that closely round up to the given width and length
% hole dimensions 8x8 microns
% hole pitch 6 um
% offset allows for an offset on any edge of the bar, each entry in offset 
% read clockwise from the origin, pos
% rot rotates the bar. an unrotated bar has the origin on the center of the
% leftmost face. roatations are implemented CCW from this orientation
function [gsbar,dim] = make_bar(pos,width,length,offset,rot,label,layers)
    % hole dimension
    holeS=8;
    % hole spaceing
    lineS=6;
    % number of holes and lines according to given width dimension
    nholesW=round((width-lineS-offset(2)-offset(4))/(holeS+lineS));
    nlinesW=nholesW+1;
    % actual width dimension
    widthAct=nholesW*holeS+nlinesW*lineS+offset(2)+offset(4);
    % number of holes and lines according to given length dimension
    nholesL=round((length-lineS-offset(1)-offset(3))/(holeS+lineS));
    nlinesL=nholesL+1;
    % actual length dimension
    lengthAct=nholesL*holeS+nlinesL*lineS+offset(1)+offset(3);

    % define bar according to dimensions and rotation
    gsbar=gds_structure(label);
    xy1=pos+[widthAct/2*sind(rot),-widthAct/2*cosd(rot)];
    xy2=xy1+[lengthAct*cosd(rot),lengthAct*sind(rot)];
    xy3=xy2+[-widthAct*sind(rot),widthAct*cosd(rot)];
    xy4=xy3+[-lengthAct*cosd(rot),-lengthAct*sind(rot)];
    xy=[xy1 ; xy2 ; xy3 ; xy4 ; xy1];
    dim=[widthAct,lengthAct,rot];

    gsbar(end+1)=gds_element('boundary', 'xy',xy,'layer',layers(1));
    
    %set position to begin hole array
    holePos=xy1+[(offset(1)+lineS)*cosd(rot)-(offset(2)+lineS)*sind(rot),(offset(1)+lineS)*sind(rot)+(offset(2)+lineS)*cosd(rot)];
    % define holes in the bar
    for i=1:nholesW
        for j=1:nholesL
            xy1=holePos+[(j-1)*(holeS+lineS)*cosd(rot),(j-1)*(holeS+lineS)*sind(rot)];
            xy2=xy1+[holeS*cosd(rot),holeS*sind(rot)];
            xy3=xy2+[-holeS*sind(rot),holeS*cosd(rot)];
            xy4=xy3+[-holeS*cosd(rot),-holeS*sind(rot)];
            xy=[xy1 ; xy2 ; xy3 ; xy4 ; xy1];
            gsbar(end+1)=gds_element('boundary', 'xy',xy,'layer',layers(2));
        end
        holePos=holePos+[-(holeS+lineS)*sind(rot),(holeS+lineS)*cosd(rot)];
    end
end