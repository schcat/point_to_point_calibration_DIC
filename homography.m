%��ⵥӦ�����Ż�,ÿһ��ͼ��ȡ49����

function f=homography(m,M);     %MΪ�������꣬mΪͼ�񽹵�����
p=ones(98,9);
for i=1:49      %ÿ��ͼ��ȡ49����
    p(2*i-1,:)=[M(:,i)' 0 0 0 -m(1,i)*M(:,i)'];
    p(2*i,:)=  [0 0 0 M(:,i)' -m(2,i)*M(:,i)'];
end;
[U S V]=svd(p); %�����ֽ�
H=V(:,9);       %ȡV�����һ��
H=H/H(9);       %��һ��
options = optimset('LargeScale','off','Algorithm','levenberg-marquardt');% ʹ��lsqnonlin���з�������С�������
x = lsqnonlin( @fun1, reshape(H,1,9) , [],[],options, m, M);
H=reshape(x,3,3);		  % ������H�ָ�Ϊ3*3
H=H';
f=H;
    