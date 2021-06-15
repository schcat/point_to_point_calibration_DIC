function  patent2_InAndOut_test(file)

global COUNT;
n = 20;

res_camera_w =1920;
res_camera_h = 1080;
height_ref = 1000;
width_ref = 1000;
pad = 100;

if 1

load(file);
for k=1:1:n
    I = imread(['image/ori/test_2_1/cali_',num2str(k),'.png']);
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
    imwrite(I_correct/max(I_correct(:)), ['image/ori/test_2_1/correct_',num2str(k),'.png']);
end

pathname = 'image/';
filename_ref = 'speckle_pattern_4000_pad_0111_20_15000_1000.png';
filename_roi = 'ROI_1000_1000_100.png';
radius = 40; %DIC subset radius
spacing = 0;

parpool(20);
mat_name = ['test_2/speckle_correct_'];
parfor k=1:1:n
tic
    filename_cur = ['ori/test_2_1/correct_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
toc
end
delete(gcp('nocreate'));

fid = fopen('speckle_correct.txt','w');
for k = 1:1:n
    str=strcat('results/ori_ref_mat/test_2/speckle_correct_',num2str(k),'.mat');
    load(str);
    for v_loc = 200:100:800
        for u_loc = 200:100:800
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

if 0
[A, Rm] = fun_para_readtxt();
end

if 0
parpool(20);
COUNT = 0;
%para = load('para_temp.txt');
para=[Rm,A(1,1),A(1,3),A(2,2),A(2,3)];%�Ż��ڲ������������ƽ�ƾ��������ת�ǣ�
%options = optimset('Algorithm','levenberg-marquardt','InitDamping',1e2,'Display','iter');
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','FunctionTolerance',1e-10,'Display','iter');

[x,nor,res,exitflag,output] = lsqnonlin( @fun_mapopt, para, [],[],options);
A=[x(n*6+1) 0 x(n*6+2); 0 x(n*6+3) x(n*6+4); 0,0,1];
disp('map�����Ż������ڲ���Ϊ');
disp(A);
res;
sqrt(nor/(49*n))
exitflag;
output;
delete(gcp('nocreate'));
end
