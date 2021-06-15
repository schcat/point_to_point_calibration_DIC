x = [0,1];
y = [4210.15,3379.05];
h = fill([x,fliplr(x)],[y,0*ones(1,length(y))],'r');
set(h,'edgealpha',0,'facealpha',0.3)

%使用 fill 函数绘制置信范围，同时使用 plot 函数绘制数据点，以此方式创建含有置信范围的绘图
x = 0:0.2:10;                     
y = besselj(0, x);
xconf = [x x(end:-1:1)] ;%一个来回         
yconf = [y+0.15 y(end:-1:1)-0.15];%0.15就是条带宽度，换成矩阵就会有不同的宽度
p = fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColor为填充颜色，EdgeColor为边框颜色
hold on
plot(x,y,'ro')
