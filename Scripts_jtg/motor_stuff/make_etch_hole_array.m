function make_etch_hole_array(pos,width,space,Nrow,Ncol)
    gs_etch_holes=gds_structure('etch_hole_arr');
    orig=pos-[Ncol*width-(Ncol+1)*space];
    xy1=[x,y;x,y+width;x+length_tot,y+width;x+length_tot,y;x,y];
    gs(end+1)=gds_element('boundary', 'xy',xy1, 'layer',1);
    for i=1:N-1
        if mod(i,2)==0
            y=y+width+gap2;
            x=x-length_support;
            xy1=[x,y;x,y+width;x+length_tot,y+width;x+length_tot,y;x,y];
            gs(end+1) = gds_element('boundary', 'xy',xy1, 'layer',1);
        else
            y=y+width+gap1;
            x=x+length_support;
            xy1=[x,y;x,y+width;x+length_tot,y+width;x+length_tot,y;x,y];
            gs(end+1) = gds_element('boundary', 'xy',xy1, 'layer',1);
        end
    end
    %create a library to hold the structure
    glib = gds_library('comb_array', 'uunit',1e-6, 'dbunit',1e-9, gs);
    %finally write the library to a file
    write_gds_library(glib, 'comb_array.gds');
end
