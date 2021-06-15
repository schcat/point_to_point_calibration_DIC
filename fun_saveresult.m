%保存结果，包括最终结果map和n幅图的map误差

function fun_saveresult(para);
res_camera = 1024;
n=15;
R=[];
for i=1:n               %下面用来恢复出外参数矩阵
      R_new=para([(i-1)*6+1 : (i-1)*6+6] );
%      RL = vision.internal.calibration.rodriguesVectorToMatrix(R_new([1:3])');
      Q1=R_new(1);
      Q2=R_new(2);
      Q3=R_new(3);
      RL=[cos(Q2)*cos(Q3)    sin(Q1)*sin(Q2)*cos(Q3)-cos(Q1)*sin(Q3)    cos(Q1)*sin(Q2)*cos(Q3)+sin(Q1)*sin(Q3); %旋转矩阵
          cos(Q2)*sin(Q3)    sin(Q1)*sin(Q2)*sin(Q3)+cos(Q1)*cos(Q3)    cos(Q1)*sin(Q2)*sin(Q3)-sin(Q1)*cos(Q3); 
          -sin(Q2)           sin(Q1)*cos(Q2)                            cos(Q1)*cos(Q2)]
      TL=R_new([4:6])'
      RT=[RL(:,1:2) , TL];
      R=[R;RT];
end
A=[para(n*6+1) 0 para(n*6+2); 0 para(n*6+3) para(n*6+4); 0,0,1];%分解出内参数

%把标准图像和掩膜做投影

value = imread('image/speckle_pattern_2048_40_pad_12_100000.png');
value = value(41:2088,41:2088);
for k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(value, A, RT, k);
    imwrite(mat2gray(result), ['image/ori_0203/simulate_speckle_patent_proj_',num2str(k),'.png']);
end
roi = imread('image/ROI_2048_2048_40.png');
roi = roi(41:2088,41:2088);
for k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(roi, A, RT, k);
    imwrite(mat2gray(result), ['image/ori_0203/ROI_proj_',num2str(k),'.png']);
end
%DIC计算
pathname = 'image/ori_0203/';
radius = 10; %DIC subset radius
% 如果ROI的pad有变化，记得在ncorr_auto_initseeds.m里面改pad
for k=1:1:n
    filename_ref = ['simulate_speckle_patent_proj_',num2str(k),'.png'];
    filename_roi = ['ROI_proj_',num2str(k),'.png'];
    filename_cur = ['simulate_speckle_patent_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius);
    save(['results/ori_ref_mat/speckle_map_',num2str(k),'.mat'],'displacements');
end

mean_map_u = zeros(res_camera);
mean_map_v = zeros(res_camera);
plot_u = zeros(res_camera,res_camera,n);
plot_v = zeros(res_camera,res_camera,n);
for k = 1:1:n
    str=strcat('results/ori_ref_mat/speckle_map_',num2str(k),'.mat');
    load(str);
    plot_u(:,:,k) = displacements.plot_u;
    plot_v(:,:,k) = displacements.plot_v;
end
count_sum = 0;
for v_loc = 1:1:res_camera
    for u_loc = 1:1:res_camera
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
save(['results/speckle_map.mat'],'mean_map_u','mean_map_v');

for k = 1:1:n
    u_bias_err = zeros(res_camera);
    v_bias_err = zeros(res_camera);
    for v_loc = 1:1:res_camera
        for u_loc = 1:1:res_camera
            u_bias = plot_u(v_loc,u_loc,k);
            v_bias = plot_v(v_loc,u_loc,k);
            if u_bias~=0 || v_bias~=0
                u_bias_err(v_loc,u_loc) = u_bias - mean_map_u(v_loc,u_loc);
                v_bias_err(v_loc,u_loc) = v_bias - mean_map_v(v_loc,u_loc);
            end
        end
    end
    save(['results/ori_ref_mat/speckle_map_err_',num2str(k),'.mat'],'u_bias_err','v_bias_err');
end

