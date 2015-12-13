function [p,t]=pacman(h)

fd=@(p) ddiff(dcircle(p,0,0,1),drectangle(p,0,1,0,1));
fh = @huniform;
[p,t]=distmesh2d(fd,fh,h,[-3,-3;3,3],[0,0;0,1;1,0]);

end