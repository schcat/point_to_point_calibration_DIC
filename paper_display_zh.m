value = load('map_x.txt');
x = [1:4:255];
x1 = [1:1:255];
mycolor = interp1(x,jet,x1);
seed_x = [1:1:1920];
aid_x = repmat(seed_x,1080,1);
seed_y = [1:1:1080]';
aid_y = repmat(seed_y,1,1920);
value = value - aid_x;
[X, Y] = meshgrid(1:1920,1:1080);
%surf(X, Y, mean_map_u);
%view(0,270)
pcolor(X, Y, value);
shading interp; 
colorbar; colormap(mycolor);
axis equal
axis([0 1920,0 1080])
xlabel('x£¨pixel£©');ylabel('y£¨pixel£©');
set(gca,'YDir','reverse')