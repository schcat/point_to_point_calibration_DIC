load('speckle_map_test.mat');
value = load('map_y.txt');
x = [1:4:255];
x1 = [1:1:255];
mycolor = interp1(x,jet,x1);

seed_x = [1:1:1920];
aid_x = repmat(seed_x,1080,1);
seed_y = [1:1:1080]';
aid_y = repmat(seed_y,1,1920);
value = value - aid_y;

result = mean_map_v - value;

[X, Y] = meshgrid(1:1920,1:1080);
pcolor(X, Y, result);

shading interp; 
colorbar;
colormap(mycolor);
axis equal
axis([0 1920,0 1080])
xlabel('x£¨pixel£©');ylabel('y£¨pixel£©');
set(gca,'YDir','reverse')
