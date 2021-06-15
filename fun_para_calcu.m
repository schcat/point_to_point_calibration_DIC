function [A, Rm] = fun_para_calcu() 
n = 20;
% read world coordinate and key points coordinate
M=load('board.txt');       % world coordinate
M=[M';ones(1,49)];
m_all=load('speckle_ori.txt');       % key points coordinate
m_one=ones(3,49,n);
for i=1:1:n
    m_temp = m_all((i-1)*49+1:i*49,:);
    m_one(:,:,i) = [m_temp';ones(1,49)];
end

% solve homography matrix and get b
H = ones(3,3,n);   % homography matrix for each image
for i=1:1:n
    H(:,:,i)=homography(m_one(:,:,i),M);
end

V=ones(2*n,6);           % V*b=0
for i=1:n               % get V
    V(2*i-1,:)=[H(1,1,i)*H(1,2,i) H(1,1,i)*H(2,2,i)+H(2,1,i)*H(1,2,i) H(2,1,i)*H(2,2,i) ...
                H(3,1,i)*H(1,2,i)+H(1,1,i)*H(3,2,i) H(3,1,i)*H(2,2,i)+H(2,1,i)*H(3,2,i) H(3,1,i)*H(3,2,i)];
    p1=[H(1,1,i)^2 H(1,1,i)*H(2,1,i)+H(2,1,i)*H(1,1,i) H(2,1,i)^2 H(3,1,i)*H(1,1,i)+H(1,1,i)*H(3,1,i) H(3,1,i)*H(2,1,i)+H(2,1,i)*H(3,1,i) H(3,1,i)^2];
    p2=[H(1,2,i)^2 H(1,2,i)*H(2,2,i)+H(2,2,i)*H(1,2,i) H(2,2,i)^2 H(3,2,i)*H(1,2,i)+H(1,2,i)*H(3,2,i) H(3,2,i)*H(2,2,i)+H(2,2,i)*H(3,2,i) H(3,2,i)^2];
    V(2*i,:)=p1-p2;
end;

[u s v]=svd(V);        % get b by orthogonal decomposition
b=v(:,6);

% get internal parameter matrix K
cy=(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2);
lamda=b(6)-(b(4)^2+(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2)*(b(2)*b(4)-b(1)*b(5)))/b(1);
kx=sqrt(lamda/b(1));
ky=sqrt((lamda*b(1))/(b(1)*b(3)-b(2)^2));
gama=-(b(2)*kx^2*ky)/lamda;
cx=(gama*cy)/kx-(b(4)*kx^2)/lamda;
%A=[kx gama cx;0 ky cy;0 0 1];      % get K
A=[kx 0 cx;0 ky cy;0 0 1]      % get K
mp=ones(3,49,n);

% get external parameter matrix RT
R=[];
Rm=[];
for i=1:n
    q1=(norm(inv(A)*H(:,1,i))+norm(inv(A)*H(:,2,i)))/2;
    r1=(1/q1)*inv(A)*H(:,1,i);
    r2=(1/q1)*inv(A)*H(:,2,i);
    r3=cross(r1,r2);
    RL=[r1 r2 r3];
    [u s v]=svd(RL);
    RL=u*v';                       % rotation matrix orthogonalization
    p=(1/q1)*inv(A)*H(:,3,i);
    
    RT=[r1 r2 p];
    R=[R;RT];

    sq = sqrt(RL(3,2)^2+RL(3,3)^2);
    Q1 = atan2(RL(3,2),RL(3,3));
    Q2 = atan2(-RL(3,1),sq);
    Q3 = atan2(RL(2,1),RL(1,1));
    
    [Q1 Q2 Q3]';
    p;
    R_new=[Q1,Q2,Q3,p'];
%    R_new=[rotationVectors',p'];
    Rm=[Rm , R_new];
    
end

% reprojection error of closed solutions
mp = m_one;
res = ones(3,49,n);
for i=1:n
    RT=R([(i-1)*3+1 : (i-1)*3+3],:);
    x=A*RT*M;
    x=[x(1,:)./x(3,:) ; x(2,:)./x(3,:); x(3,:)./x(3,:)];	% make sure the last column is 1
    res(:,:,i) = mp(:,:,i)-x;
end

res_sum = zeros(1,98*n);
for i = 1:n
    for j = 1:49
        res_sum(1,2*((i-1)*49+j)-1) = res(1,j,i);
        res_sum(1,2*((i-1)*49+j)) = res(2,j,i);
    end
end
res_sum.*res_sum;
xxx=sum(res_sum.*res_sum);
initres = sqrt(sum(res_sum.*res_sum,2)/(n*49))
