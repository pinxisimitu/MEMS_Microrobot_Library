%pos - origin, set at the front center face of the shuttle
%fingerW - GCA finger width
%fingerL - GCA finger overlap length
%fingerSupportL - support at the base of GCA fingers
%gap1 - front gap for GCA
%gap2 - back gap for GCA
%N - number of fingers
%travel - maximum shuttle travel
%armW - width of angled arm
%armL - length of angled arm
%armAngle - angle of the angled arm
%springW - support spring width
%springL - support spring length
%gstopGap - gap stop gap
%label - iw name
function gs=make_iw_motor(pos,toothW,toothL,toothS,shuttleG,fingerW,fingerL,fingerSupportL,gap1,gap2,N,travel,armW,armL,armAngle,springW,springL,gstopGap,label)
    % edge offset for rectangle shapes. see make_bar function
    offset=[0,0,0,0];

    % minimum anchor width
    anchorW=50;

    % default stucture spacing (spacing between structures, spacing around
    % springs, etc.
    space=5;
    
    % routing width
    routeW=30;
    
    % stator anchor length
    anchorL=N/2*(2*fingerW+gap1+gap2)-gap2+space;
    
    % gap stop parameters
    % gapstopL, stop length
    % gstopW, stop width
    % gstopBumpW, stop bump width
    % gstopBumpL, stop bump length
    % gstopSpace, pitch between stop bumps
    gstopL=110;
    gstopW=10;
    gstopBumpW=2;
    gstopBumpL=5;
    gstopSpace=15;
    
    % mid-bar length that exceeds the anchor length, to give breadth for
    % support springs on GCA
    midBarOL=50;    
    
    % width and length of the rotor midbar, holding the rotor fingers
    % default array has 2 rows of fingers
    barW=20;
    barL=N/2*(2*fingerW+gap1+gap2)+midBarOL+gstopW+gstopBumpW+anchorW;
    
    % width and length for shuttle
    % gap between pawls and shuttle
    shuttleW=38;
    shuttleL=4*anchorW+2*(fingerL+2*fingerSupportL)+2*(barW+springL)+4*space+3*routeW+travel;
    
    %% make shuttle
    % tooth dimensions can be altered here
    [gs,dimShuttle]=make_bar(pos,shuttleW,shuttleL,offset,180,'shuttle');
    nTooth=dimShuttle(2)/(toothS+toothW);

    gsTemp=make_tooth_array(pos+[0,dimShuttle(1)/2],toothW,toothL,toothS,nTooth,0,'1');
    gs=join_gds_structures(gs,gsTemp);
    gsTemp=make_tooth_array(pos+[0,-dimShuttle(1)/2],toothW,toothL,toothS,nTooth,1,'2');
    gs=join_gds_structures(gs,gsTemp);

    %% make angled arms
    nToothPawl=2;
    pawlW=nToothPawl*toothW+toothS;
    % length of pawl head (dimension along tooth array)
    pawlL=toothW*nToothPawl+toothS*(nToothPawl-1)+0.5;
    
    % define pawl origins according to motor dimensions, offset
    % appropriately to align with teeth on the shuttle
    pawlX=anchorW+fingerL+fingerSupportL+barW/2-armW/sind(armAngle)/2-armL*cosd(armAngle)-pawlL-(pawlL-armW/sind(armAngle));
    pawlX=pawlX-mod(pawlX,toothW+toothS)-toothW/2-toothS/2;
    pawlY=dimShuttle(1)/2+2*toothL+shuttleG;
    
    % offset in the y-dimension for rear pawls
    rearPawlX=anchorW+fingerL+fingerSupportL+barW+2*(springL+anchorW)+3*routeW+4*space+barW/2-armW/sind(armAngle)/2-armL*cosd(armAngle)-pawlL-(pawlL-armW/sind(armAngle));
    sepErr=mod(rearPawlX,toothW+toothS);
    spaceRoute=space-sepErr/4;
    rearPawlX=rearPawlX-sepErr;
    
    % define angled arm gds structs
    pawlPos=pos-[pawlX -pawlY];
    gsTemp=make_angled_arm(pawlPos,armW,armL,nToothPawl,pawlW,pawlL,toothW,toothL,toothS,armAngle,0,'arm1');
    gs=join_gds_structures(gs,gsTemp);
    pawlPos=pos-[pawlX pawlY];
    gsTemp=make_angled_arm(pawlPos,armW,armL,nToothPawl,pawlW,pawlL,toothW,toothL,toothS,armAngle,1,'arm1');
    gs=join_gds_structures(gs,gsTemp);
    pawlPos=pos-[rearPawlX -pawlY];
    gsTemp=make_angled_arm(pawlPos,armW,armL,nToothPawl,pawlW,pawlL,toothW,toothL,toothS,armAngle,0,'arm1');
    gs=join_gds_structures(gs,gsTemp);
    pawlPos=pos-[rearPawlX pawlY];
    gsTemp=make_angled_arm(pawlPos,armW,armL,nToothPawl,pawlW,pawlL,toothW,toothL,toothS,armAngle,1,'arm1');
    gs=join_gds_structures(gs,gsTemp);

    %% make gap closer arrays
    % define GCA origin
    gcaX=pawlX+armL*cosd(armAngle)+pawlL-armW/sind(armAngle)/2;
    gcaY=pawlY+pawlW+armL*sind(armAngle);
    reargcaX=rearPawlX+armL*cosd(armAngle)+pawlL-armW/sind(armAngle)/2;
    gcaPos=pos-[gcaX -gcaY];
    gsTemp=make_GCA_array(gcaPos,fingerW,fingerL,fingerSupportL,gap1,gap2,N,springW,springL,gstopGap,anchorW,space,anchorL,90,0,'blah');
    gs=join_gds_structures(gs,gsTemp);
    gcaPos=pos-[reargcaX -gcaY];
    gsTemp=make_GCA_array(gcaPos,fingerW,fingerL,fingerSupportL,gap1,gap2,N,springW,springL,gstopGap,anchorW,space,anchorL,90,1,'blah');
    gs=join_gds_structures(gs,gsTemp);
    gcaPos=pos-[gcaX gcaY];
    gsTemp=make_GCA_array(gcaPos,fingerW,fingerL,fingerSupportL,gap1,gap2,N,springW,springL,gstopGap,anchorW,space,anchorL,270,0,'blah');
    gs=join_gds_structures(gs,gsTemp);
    gcaPos=pos-[reargcaX gcaY];
    gsTemp=make_GCA_array(gcaPos,fingerW,fingerL,fingerSupportL,gap1,gap2,N,springW,springL,gstopGap,anchorW,space,anchorL,270,1,'blah');
    gs=join_gds_structures(gs,gsTemp);
    
    
    % uncomment to test function on its own
    %glib=gds_library(strcat('iw_',label), 'uunit',1e-6, 'dbunit',1e-9,gs);
    %write_gds_library(glib, strcat('iw_',label,'.gds'));


end
