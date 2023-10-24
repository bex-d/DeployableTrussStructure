
function myPlot(x,ix)

NN = length(x);
NT = length(ix);
    
    colours = ['k','k','k','k','r','r','m','g','b']; %vertical, horizontal, horizontal ring, outside ring diagonals, inside ring diagonals
    
for i=1:NT
    
    n1 = ix(i,1);
    n2 = ix(i,2);
    t  = ix(i,3);
    
    xp = x([n1,n2],1);
    yp = x([n1,n2],2);
    zp = x([n1,n2],3);
    
    col = colours(t);
    
    plot3(xp,yp,zp,col);
    hold all;

end
         
axis equal;

xmin = min(x(1:NN,1));
ymin = min(x(1:NN,2));
zmin = min(x(1:NN,3));
xmax = max(x(1:NN,1));
ymax = max(x(1:NN,2));
zmax = max(x(1:NN,3));
dx = (xmax - xmin) * .3;
dy = (ymax - ymin) * .3;
dz = (zmax - zmin) * .3;

axis([xmin-dx xmax+dx ymin-dy ymax+dy zmin-dz zmax+dz]);

end


