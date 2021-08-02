function [A, Rm] = fun_para_readtxt() 
n = 20;
% read world coordinate and key points coordinate
M=load('board.txt');       % world coordinate
M=[M';ones(1,49)];
m_all=load('../result/train_2_6/speckle_ori_70.txt');       % key points coordinate
m_one=ones(3,49,n);
for i=1:1:n
    m_temp = m_all((i-1)*49+1:i*49,:);
    m_one(:,:,i) = [m_temp';ones(1,49)];
end

% get internal parameter matrix K
% get external parameter matrix RT
R=[];
Rm=[];
A = load('../result/train_2_6/para_k.txt');
r_read = load('../result/train_2_6/para_r.txt');
t_read = load('../result/train_2_6/para_t.txt');

for i=1:n
    r_vect = r_read(i,:);
    t_vect = t_read(i,:);
    RL = Rodrigues(r_vect);
    
    RT=[RL(:,1:2) t_vect'];
    R=[R;RT];

    sq = sqrt(RL(3,2)^2+RL(3,3)^2);
    Q1 = atan2(RL(3,2),RL(3,3));
    Q2 = atan2(-RL(3,1),sq);
    Q3 = atan2(RL(2,1),RL(1,1));
    
    R_new=[Q1,Q2,Q3,t_vect];
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
