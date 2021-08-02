%�����Ż������������?

function f=fun_mapopt(para);
global COUNT;
if mod(COUNT,125) == 0
    fun_save_para(para, COUNT);
%    fun_test(['para_temp_',num2str(COUNT/125),'.txt']);
%    fun_test_simp(['para_temp_',num2str(COUNT/125),'.txt']);
%    fun_test_fullpixel_simp(['para_temp_',COUNT,'.txt']);
end
COUNT = COUNT + 1;
res_camera_w = 1920;
res_camera_h = 1080;
n=20;
R=[];
for i=1:n               %���������ָ������������?
      R_new=para([(i-1)*6+1 : (i-1)*6+6] );
%      RL = vision.internal.calibration.rodriguesVectorToMatrix(R_new([1:3])');
      Q1=R_new(1);
      Q2=R_new(2);
      Q3=R_new(3);
      TL=R_new([4:6])';    
      RL=[cos(Q2)*cos(Q3)    sin(Q1)*sin(Q2)*cos(Q3)-cos(Q1)*sin(Q3)    cos(Q1)*sin(Q2)*cos(Q3)+sin(Q1)*sin(Q3); %��ת����
          cos(Q2)*sin(Q3)    sin(Q1)*sin(Q2)*sin(Q3)+cos(Q1)*cos(Q3)    cos(Q1)*sin(Q2)*sin(Q3)-sin(Q1)*cos(Q3); 
          -sin(Q2)           sin(Q1)*cos(Q2)                            cos(Q1)*cos(Q2)];
      RT=[RL(:,1:2) , TL];
      R=[R;RT];
end
A=[para(n*6+1) 0 para(n*6+2); 0 para(n*6+3) para(n*6+4); 0,0,1];%�ֽ���ڲ���?

%�ѱ�׼ͼ�����Ĥ��Ͷ�?

value = imread('../image/speckle_pattern_4000_pad_0111_20_15000_2000.png');
value = value(200:1800,200:1800);
%tic
parfor k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(value, A, RT, k);
    imwrite(mat2gray(result), ['../image/train_2_6/simulate_speckle_patent_proj_',num2str(k),'.png']);
end
roi = imread('../image/ROI_2000_2000_200.png');
roi = roi(200:1800,200:1800);
parfor k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(roi, A, RT, k);
    imwrite(mat2gray(result), ['../image/train_2_6/ROI_proj_',num2str(k),'.png']);
end
%toc
%DIC����
pathname = '../image/train_2_6/';
radius = 70; %DIC subset radius
spacing = 9;
% ���?ROI��pad�б仯���ǵ���ncorr_auto_initseeds.m�����?pad
tic
mat_name = ['../result/train_2_6/speckle_map_'];
parfor k=1:1:n
    filename_ref = ['simulate_speckle_patent_proj_',num2str(k),'.png'];
    filename_roi = ['ROI_proj_',num2str(k),'.png'];
    filename_cur = ['cali_',num2str(k),'.png'];
    displacements(k) = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
 %   save(['results/ori_ref_mat/speckle_map_',num2str(k),'.mat'],'displacements');
end
toc
mean_map_u = zeros(res_camera_h,res_camera_w);
mean_map_v = zeros(res_camera_h,res_camera_w);
plot_u = zeros(res_camera_h,res_camera_w,n);
plot_v = zeros(res_camera_h,res_camera_w,n);
plot_u_temp = zeros(res_camera_h,res_camera_w,n);
plot_v_temp = zeros(res_camera_h,res_camera_w,n);
%tic
for k = 1:1:n
    str=strcat('../result/train_2_6/speckle_map_',num2str(k),'.mat');
    load(str);
    plot_u(:,:,k) = displacements.plot_u;
    plot_v(:,:,k) = displacements.plot_v;
    plot_u_temp(:,:,k) = displacements.plot_u;
    plot_v_temp(:,:,k) = displacements.plot_v;
    plot_coef(:,:,k) = displacements.plot_corrcoef;
%    edge_u(:,:,k) = edge(plot_u(:,:,k));
%    edge_v(:,:,k) = edge(plot_v(:,:,k));
%    imwrite(edge_u(:,:,k), ['image/ori_0315/edge_u_',num2str(k),'.png']);
%    imwrite(edge_v(:,:,k), ['image/ori_0315/edge_v_',num2str(k),'.png']);
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
    for v_loc = 1:1:res_camera_h
        for u_loc = 1:1:res_camera_w
            if v_loc < 21 || v_loc > res_camera_h - 20 || u_loc < 21 || u_loc > res_camera_w - 20
                plot_u(v_loc, u_loc, k) = 0;
                plot_v(v_loc, u_loc, k) = 0;
            end
        end
    end
    %    imwrite(plot_u(:,:,k), ['image/ori_0315/delete_edge_u_',num2str(k),'.png']);
%    imwrite(plot_v(:,:,k), ['image/ori_0315/delete_edge_v_',num2str(k),'.png']);
end
%toc
count_sum = 0;
for v_loc = 1:1:res_camera_h
    for u_loc = 1:1:res_camera_w
        count = 0;
        bias_sum_u = 0;
        bias_sum_v = 0;
        for k = 1:1:n
            u_bias = plot_u(v_loc,u_loc,k);
            v_bias = plot_v(v_loc,u_loc,k);
	    coef = plot_coef(v_loc, u_loc, k);
            if (u_bias~=0 || v_bias~=0) && coef < 0.065
%            if u_bias~=0 || v_bias~=0
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
save(['../result/train_2_6/speckle_map.mat'],'mean_map_u','mean_map_v');

res_bias = zeros(1,res_camera_h*res_camera_w*n/50);
for v_loc = 10:10:res_camera_h
    for u_loc = 10:10:res_camera_w
        for k = 1:1:n
            u_bias = plot_u(v_loc,u_loc,k);
            v_bias = plot_v(v_loc,u_loc,k);
	    coef = plot_coef(v_loc, u_loc, k);
            if (u_bias~=0 || v_bias~=0) && coef < 0.065
%            if u_bias~=0 || v_bias~=0
                res_bias(1,2*((v_loc*res_camera_w/10+u_loc)*n/10+k)-1) = u_bias - mean_map_u(v_loc,u_loc);
                res_bias(1,2*((v_loc*res_camera_w/10+u_loc)*n/10+k)) = v_bias - mean_map_v(v_loc,u_loc);
            end
        end
    end
end
size(res_bias);
count_sum;
reproj=sqrt(sum(res_bias.*res_bias,2)/(res_camera_h*res_camera_w*n/100));
f = res_bias;

