%求解单应矩阵并优化,每一幅图像取49个点

function f=homography(m,M);     %M为世界坐标，m为图像焦点坐标
p=ones(98,9);
for i=1:49      %每幅图像取49个点
    p(2*i-1,:)=[M(:,i)' 0 0 0 -m(1,i)*M(:,i)'];
    p(2*i,:)=  [0 0 0 M(:,i)' -m(2,i)*M(:,i)'];
end;
[U S V]=svd(p); %正交分解
H=V(:,9);       %取V的最后一列
H=H/H(9);       %归一化
options = optimset('LargeScale','off','Algorithm','levenberg-marquardt');% 使用lsqnonlin进行非线性最小二乘求解
x = lsqnonlin( @fun1, reshape(H,1,9) , [],[],options, m, M);
H=reshape(x,3,3);		  % 将矩阵H恢复为3*3
H=H';
f=H;
    