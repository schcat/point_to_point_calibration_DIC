%��ͶӰ���ͼ�?��map������У�����ͽ������ѱ궨���궨
function fun_test(file);
%��ȡmap

%[A, Rm] = fun_para_readtxt();
%para=[Rm,A(1,1),A(1,3),A(2,2),A(2,3)];
%para = load('para_temp_25.txt');
para = load(file);

res_camera_w =1920;
res_camera_h = 1080;
n=20;

if 0

R=[];
for i=1:n               %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½?¸ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
      R_new=para([(i-1)*6+1 : (i-1)*6+6] );
%      RL = vision.internal.calibration.rodriguesVectorToMatrix(R_new([1:3])');
      Q1=R_new(1);
      Q2=R_new(2);
      Q3=R_new(3);
      TL=R_new([4:6])';
      RL=[cos(Q2)*cos(Q3)    sin(Q1)*sin(Q2)*cos(Q3)-cos(Q1)*sin(Q3)    cos(Q1)*sin(Q2)*cos(Q3)+sin(Q1)*sin(Q3); %ï¿½ï¿½?ªï¿½ï¿½ï¿½ï¿½
          cos(Q2)*sin(Q3)    sin(Q1)*sin(Q2)*sin(Q3)+cos(Q1)*cos(Q3)    cos(Q1)*sin(Q2)*sin(Q3)-sin(Q1)*cos(Q3);
          -sin(Q2)           sin(Q1)*cos(Q2)                            cos(Q1)*cos(Q2)];
      RT=[RL(:,1:2) , TL];
      R=[R;RT];
end
A=[para(n*6+1) 0 para(n*6+2); 0 para(n*6+3) para(n*6+4); 0,0,1];%ï¿½?½ï¿½ï¿½ï¿½?²ï¿½ï¿½ï¿?

%ï¿½?±ï¿½?¼?¼ï¿½ï¿½ï¿½ï¿½ï¿½?¤ï¿½ï¿½?¶??

value = imread('../image/speckle_pattern_4000_pad_0111_20_15000_2000.png');
value = value(200:1800,200:1800);

tic
for k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(value, A, RT, k);
    imwrite(mat2gray(result), ['../image/train_2_2/simulate_speckle_patent_proj_test_',num2str(k),'.png']);
end
roi = imread('../image/ROI_2000_2000_200.png');
roi = roi(200:1800,200:1800);
for k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(roi, A, RT, k);
    imwrite(mat2gray(result), ['../image/train_2_2/ROI_proj_test_',num2str(k),'.png']);
end
toc
%DICï¿½ï¿½ï¿½ï¿½
pathname = '../image/train_2_2/';
radius = 70; %DIC subset radius
spacing = 0;
% ï¿½ï¿½ï¿?ROIï¿½ï¿½padï¿½?±ä»¯ï¿½ï¿½ï¿½?µï¿½ï¿½ï¿½ncorr_auto_initseeds.mï¿½ï¿½ï¿½ï¿½ï¿?pad
parpool(20);
mat_name = ['../result/train_2_2_70/speckle_map_test_'];
parfor k=1:1:n
tic
    filename_ref = ['simulate_speckle_patent_proj_test_',num2str(k),'.png'];
    filename_roi = ['ROI_proj_test_',num2str(k),'.png'];
    filename_cur = ['cali_',num2str(k),'.png'];
    displacements(k) = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
 %   save(['results/ori_ref_mat/speckle_map_',num2str(k),'.mat'],'displacements');
toc
end
delete(gcp('nocreate'));
end

if 0 

mean_map_u = zeros(res_camera_h,res_camera_w);
mean_map_v = zeros(res_camera_h,res_camera_w);
plot_u = zeros(res_camera_h,res_camera_w,n);
plot_v = zeros(res_camera_h,res_camera_w,n);
for k = 1:1:n
    str=strcat('../result/train_2_2_70/speckle_map_test_',num2str(k),'.mat');
    load(str);
    plot_u(:,:,k) = displacements.plot_u;
    plot_v(:,:,k) = displacements.plot_v;
    plot_u_temp(:,:,k) = displacements.plot_u;
    plot_v_temp(:,:,k) = displacements.plot_v;
end
for k = 1:1:n
    for v_loc = 21:1:res_camera_h-20
        for u_loc = 21:1:res_camera_w-20
            temp_u_x = isequal(plot_u_temp(v_loc, u_loc, k),0) + isequal(plot_u_temp(v_loc, u_loc - 1, k),0) + isequal(plot_u_temp(v_loc, u_loc + 1, k),0);
            temp_u_y = isequal(plot_u_temp(v_loc, u_loc, k),0) + isequal(plot_u_temp(v_loc - 1, u_loc, k),0) + isequal(plot_u_temp(v_loc + 1, u_loc, k),0);
            temp_v_x = isequal(plot_v_temp(v_loc, u_loc, k),0) + isequal(plot_v_temp(v_loc, u_loc - 1, k),0) + isequal(plot_v_temp(v_loc, u_loc + 1, k),0);
            temp_v_y = isequal(plot_v_temp(v_loc, u_loc, k),0) + isequal(plot_v_temp(v_loc - 1, u_loc, k),0) + isequal(plot_v_temp(v_loc + 1, u_loc, k),0);
            if temp_u_x > 0 && temp_u_x < 3 || temp_u_y > 0 && temp_u_y < 3 || temp_v_x > 0 && temp_v_x < 3 || temp_v_y > 0 && temp_v_y < 3
                plot_u(v_loc-20:v_loc+20, u_loc-20:u_loc+20, k) = 0;
                plot_v(v_loc-20:v_loc+20, u_loc-20:u_loc+20, k) = 0;
            end
        end
    end
end
count_sum = 0;
for v_loc = 1:1:res_camera_h
    for u_loc = 1:1:res_camera_w
        count = 0;
        bias_sum_u = 0;
        bias_sum_v = 0;
        for k = 1:1:n
            u_bias = plot_u(v_loc,u_loc,k);
            v_bias = plot_v(v_loc,u_loc,k);
            if u_bias~=0 || v_bias~=0
                count = count+1;
                bias_sum_u = bias_sum_u + u_bias;
                bias_sum_v = bias_sum_v + v_bias;
            end
        end
        if count~=0
            mean_map_u(v_loc,u_loc) = bias_sum_u/count;
            mean_map_v(v_loc,u_loc) = bias_sum_v/count;
        end
        count_sum = count_sum + count;
    end
end
save(['../result/train_2_2_70/speckle_map_test_70.mat'],'mean_map_u','mean_map_v');

plot_u = zeros(res_camera_h,res_camera_w);
plot_v = zeros(res_camera_h,res_camera_w);
%pathname = 'image/ori_0315/';
load('../result/train_2_2_70/speckle_map_test_70.mat');
%imshow(mean_map_u);

end

if 0

%����У��
for k=1:1:n
    I = imread(['../image/train_2_2/cali_',num2str(k),'.png']);
    I_correct = zeros(res_camera_h,res_camera_w);
    for i = 1:1:res_camera_h
        for j = 1:1:res_camera_w
            u_bias = mean_map_u(i,j);
            v_bias = mean_map_v(i,j);
            if u_bias == 0 || v_bias == 0 || j + u_bias < 1 || i + v_bias < 1 || j + u_bias > 1920 || i + v_bias > 1080
                I_correct(i,j) = I(i,j);
            else 
                uu = j + u_bias;
                vv = i + v_bias;
                %˫���Բ�ֵ
                value_1 = I(floor(vv),floor(uu));
                value_2 = I(floor(vv),ceil(uu));
                value_3 = I(ceil(vv),ceil(uu));
                value_4 = I(ceil(vv),floor(uu));
                result_1 = (ceil(vv) - vv) * value_1 + (vv - floor(vv)) * value_4;
                result_2 = (ceil(vv) - vv) * value_2 + (vv - floor(vv)) * value_3;
                result_value = (ceil(uu) - uu) * result_1 + (uu - floor(uu)) * result_2;
            
                I_correct(i,j) = result_value;
            end
        end
    end
%    imshow(I_correct/max(I_correct(:)));
    imwrite(I_correct/max(I_correct(:)), ['../image/train_2_2/correct_',num2str(k),'.png']);
end

end

%��ȡ���Ƶ�
height_ref = 1000;
width_ref = 1000;
pad = 100;
pathname = '../image/';
filename_ref = 'speckle_pattern_4000_pad_0111_20_15000_1000.png';
filename_roi = 'ROI_1000_1000_100.png';
radius = 70; %DIC subset radius
spacing = 0;
% ���?ROI��pad�б仯���ǵ���ncorr_auto_initseeds.m�����?pad

parpool(20);
mat_name = ['../result/train_2_2_70/speckle_correct_'];
parfor k=1:1:n
tic
    filename_cur = ['train_2_2/correct_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
%    save(['results/ori_ref_mat/speckle_correct_',num2str(k),'.mat'],'displacements');
toc
end
delete(gcp('nocreate'));

fid = fopen('../result/train_2_2_70/speckle_correct.txt','w');
for k = 1:1:n
    str=strcat('../result/train_2_2_70/speckle_correct_',num2str(k),'.mat');
    load(str);
    for v_loc = 200:100:800
        for u_loc = 200:100:800
       % 19.2 = (1728-2*96)/(140+1)��ȥ���߿����Ƶ��ļ��?
            u_bias = displacements.plot_u(v_loc,u_loc);
            v_bias = displacements.plot_v(v_loc,u_loc);
            uu = u_loc + u_bias;
            vv = v_loc + v_bias;
            fprintf(fid, '%f %f\n',uu,vv);
        end
    end
end
fclose(fid);

M=load('board.txt');       %��ȡ��������
M=[M';ones(1,49)];
m_all=load('../result/train_2_2_70/speckle_correct.txt');       %��ȡ�ؼ�������
m_one=ones(3,49,n);
for i=1:1:n
    m_temp = m_all((i-1)*49+1:i*49,:);
    m_one(:,:,i) = [m_temp';ones(1,49)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ������ⵥ�?���󲢽��?b                                      %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

H = ones(3,3,n);   %ÿ��ͼһ����Ӧ����
for i=1:1:n
    H(:,:,i)=homography(m_one(:,:,i),M);
end

V=ones(2*n,6);           %V*b=0
for i=1:n               %��V
    V(2*i-1,:)=[H(1,1,i)*H(1,2,i) H(1,1,i)*H(2,2,i)+H(2,1,i)*H(1,2,i) H(2,1,i)*H(2,2,i) ...
                H(3,1,i)*H(1,2,i)+H(1,1,i)*H(3,2,i) H(3,1,i)*H(2,2,i)+H(2,1,i)*H(3,2,i) H(3,1,i)*H(3,2,i)];
    p1=[H(1,1,i)^2 H(1,1,i)*H(2,1,i)+H(2,1,i)*H(1,1,i) H(2,1,i)^2 H(3,1,i)*H(1,1,i)+H(1,1,i)*H(3,1,i) H(3,1,i)*H(2,1,i)+H(2,1,i)*H(3,1,i) H(3,1,i)^2];
    p2=[H(1,2,i)^2 H(1,2,i)*H(2,2,i)+H(2,2,i)*H(1,2,i) H(2,2,i)^2 H(3,2,i)*H(1,2,i)+H(1,2,i)*H(3,2,i) H(3,2,i)*H(2,2,i)+H(2,2,i)*H(3,2,i) H(3,2,i)^2];
    V(2*i,:)=p1-p2;
end;

[u s v]=svd(V);        %�������ֽ���b
b=v(:,6);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                        ���·ֽ��ڲ�������                                            %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

cy=(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2);
lamda=b(6)-(b(4)^2+(b(2)*b(4)-b(1)*b(5))/(b(1)*b(3)-b(2)^2)*(b(2)*b(4)-b(1)*b(5)))/b(1);
kx=sqrt(lamda/b(1));
ky=sqrt((lamda*b(1))/(b(1)*b(3)-b(2)^2));
gama=-(b(2)*kx^2*ky)/lamda;
cx=(gama*cy)/kx-(b(4)*kx^2)/lamda;
%A=[kx gama cx;0 ky cy;0 0 1];      %����ڲ�������?
A=[kx 0 cx;0 ky cy;0 0 1]      %����ڲ�������?
mp=ones(3,49,n);                %ͼ��������󣨺ϳ����?��ʽ������ѭ��������
mp = m_one;
%para1=[kx gama cx ky cy ];
%A = [2500, 0, 512; 0, 2500, 512; 0, 0, 1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           ����������������                                        %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
inv(A)*H(:,1,1);
inv(A)*H(:,2,1);
q1=(norm(inv(A)*H(:,1,1))+norm(inv(A)*H(:,2,1)))/2;
q1=(norm(inv(A)*H(:,1,1),2)+norm(inv(A)*H(:,2,1),2))/2;

R=[];
Rm=[];
for i=1:n
    q1=(norm(inv(A)*H(:,1,i))+norm(inv(A)*H(:,2,i)))/2;      %����������������
    r1=(1/q1)*inv(A)*H(:,1,i);
    r2=(1/q1)*inv(A)*H(:,2,i);
    r3=cross(r1,r2);
    RL=[r1 r2 r3];
    [u s v]=svd(RL);
    RL=u*v';                                             %��ת����������
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ��ʽ�����ͶӰ���            %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

res = ones(3,49,n);
for i=1:n
    RT=R([(i-1)*3+1 : (i-1)*3+3],:);
    x=A*RT*M;
    x=[x(1,:)./x(3,:) ; x(2,:)./x(3,:); x(3,:)./x(3,:)];		% ��֤���һ���?1
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

end
