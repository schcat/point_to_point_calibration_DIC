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
    imwrite(mat2gray(result), ['../image_twi/train_2_2/simulate_speckle_patent_proj_test_',num2str(k),'.png']);
end
roi = imread('../image/ROI_2000_2000_200.png');
roi = roi(200:1800,200:1800);
for k=1:1:n
    RT=R([(k-1)*3+1 : (k-1)*3+3],:);
    result = fun_projection(roi, A, RT, k);
    imwrite(mat2gray(result), ['../image_twi/train_2_2/ROI_proj_test_',num2str(k),'.png']);
end
toc
%DICï¿½ï¿½ï¿½ï¿½
pathname = '../image_twi/train_2_2/';
radius = 70; %DIC subset radius
spacing = 0;
% ï¿½ï¿½ï¿?ROIï¿½ï¿½padï¿½?±ä»¯ï¿½ï¿½ï¿½?µï¿½ï¿½ï¿½ncorr_auto_initseeds.mï¿½ï¿½ï¿½ï¿½ï¿?pad
parpool(20);
mat_name = ['../result/train_2_2_70_0080/speckle_map_test_'];
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
    str=strcat('../result/train_2_2_70_0080/speckle_map_test_',num2str(k),'.mat');
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
save(['../result/train_2_2_70_0080/speckle_map_test_70.mat'],'mean_map_u','mean_map_v');

plot_u = zeros(res_camera_h,res_camera_w);
plot_v = zeros(res_camera_h,res_camera_w);
%pathname = 'image/ori_0315/';
%imshow(mean_map_u);

end

if 0

fid = fopen('/home/wsco/cnn/calibration/camera_calibration_distance_calculate/build/applications/camera_calibration/points_correct.txt','w');
points=load('/home/wsco/cnn/calibration/camera_calibration_distance_calculate/build/applications/camera_calibration/points_with_index.txt');
load('../result/train_2_1/speckle_map_test_2.mat');
mean_map_u_1 = mean_map_u;
mean_map_v_1 = mean_map_v;
load('../result/train_2_2_70/speckle_map_test_70.mat');
mean_map_u_2 = mean_map_u;
mean_map_v_2 = mean_map_v;
for k=1:1:36030
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
%    imshow(I_correct/max(I_correct(:)));
    imwrite(I_correct/max(I_correct(:)), ['../image_twi/train_2_2/correct_',num2str(k),'.png']);
end

end

if 0

%��ȡ���Ƶ�
height_ref = 1000;
width_ref = 1000;
pad = 100;
pathname = '../image_twi/';
filename_ref = 'speckle_pattern_4000_pad_0111_20_15000_1000.png';
filename_roi = 'ROI_1000_1000_100.png';
radius = 70; %DIC subset radius
spacing = 0;
% ���?ROI��pad�б仯���ǵ���ncorr_auto_initseeds.m�����?pad

parpool(20);
mat_name = ['../result/train_2_2_70_0080/speckle_correct_'];
parfor k=1:1:n
tic
    filename_cur = ['train_2_2/correct_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
%    save(['results/ori_ref_mat/speckle_correct_',num2str(k),'.mat'],'displacements');
toc
delete(gcp('nocreate'));

end

if 1

for k = 1:1:n
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

end

