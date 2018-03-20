% define arc element with center at pos
% iradius = inner radius of arc
% width = width of the arc
% a1 = starting angle
% a2 = ending angle
function arc=arc(pos,iradius,width,a1,a2)
    arc.c=pos;
    arc.r=iradius;
    arc.w=width;
    arc.a1=a1;
    arc.a2=a2;
    arc.e=0.001;
    arc.pap=2;
end
