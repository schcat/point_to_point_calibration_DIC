load('speckle_map_test_2.mat');
x = [1:4:255];
x1 = [1:1:255];
mycolor = interp1(x,jet,x1);
[X, Y] = meshgrid(1:1920,1:1080);
%surf(X, Y, mean_map_u);
%view(0,270)
pcolor(X, Y, mean_map_v);
shading interp; 
colorbar;
colormap(mycolor);
axis equal
axis([0 1920,0 1080])
xlabel('x£¨pixel£©');ylabel('y£¨pixel£©');
set(gca,'YDir','reverse')