x = [0,1];
y = [4210.15,3379.05];
h = fill([x,fliplr(x)],[y,0*ones(1,length(y))],'r');
set(h,'edgealpha',0,'facealpha',0.3)

%ʹ�� fill �����������ŷ�Χ��ͬʱʹ�� plot �����������ݵ㣬�Դ˷�ʽ�����������ŷ�Χ�Ļ�ͼ
x = 0:0.2:10;                     
y = besselj(0, x);
xconf = [x x(end:-1:1)] ;%һ������         
yconf = [y+0.15 y(end:-1:1)-0.15];%0.15����������ȣ����ɾ���ͻ��в�ͬ�Ŀ��
p = fill(xconf,yconf,'r','FaceColor',[1 0.8 0.8],'EdgeColor','none');%FaceColorΪ�����ɫ��EdgeColorΪ�߿���ɫ
hold on
plot(x,y,'ro')
