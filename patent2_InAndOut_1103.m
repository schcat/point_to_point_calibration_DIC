function  patent2_InAndOut
%      ��m�ļ���������ڲ�������������Ż���
%      �����ƹ�������ͼ��ÿ��ͼ��ȡ10���㡣�����ÿ��ͼ��ĵ�Ӧ����Ȼ����������?�������һ��?
%      ��V*b=0���?b����Ȼ����b�ֽ���ڲ�������?A�������?A�͵�Ӧ����������������
global COUNT;
num = 60;
n = 58;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 �Ա�׼ͼƬ-����ͼƬ����DIC���?                                          %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if 1
height_ref = 1000;
width_ref = 1000;
pad = 100;
% pathname = '../image/';
pathname = '../image/';
filename_ref = 'speckle_pattern_4000_pad_0111_20_15000_1000.png';
filename_roi = 'ROI_1000_1000_100.png';
radius = num; %DIC subset radius
spacing = 0;
% ���?ROI��pad�б仯���ǵ���ncorr_auto_initseeds.m�����?pad
nn=1
mat_name = ['../result/train_dist_2/speckle_'];
tic
parpool(20);
parfor k=1:1:20
    filename_cur = ['2_select_num/cali_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
%    save(['results/ori_ref_mat/speckle_',num2str(k),'.mat'],'displacements');
end
toc
tic
parfor k=21:1:40
    filename_cur = ['2_select_num/cali_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
%    save(['results/ori_ref_mat/speckle_',num2str(k),'.mat'],'displacements');
end
toc
tic
parfor k=41:1:58
    filename_cur = ['2_select_num/cali_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
%    save(['results/ori_ref_mat/speckle_',num2str(k),'.mat'],'displacements');
end
delete(gcp('nocreate'));
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 ��DIC���������ȡ���Ƶ�?                                              %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if 1
for k = 1:1:n
fid = fopen(['../result/train_dist_2/speckle_ori_',num2str(k),'.txt'],'w');
    str=strcat('../result/train_dist_2/speckle_',num2str(k),'.mat');
    load(str);
    for v_loc = 40:20:160
        for u_loc = 40:20:160
       % 19.2 = (1728-2*96)/(140+1)��ȥ���߿����Ƶ��ļ��?
            u_bias = displacements.plot_u(v_loc * 5, u_loc * 5);
            v_bias = displacements.plot_v(v_loc * 5, u_loc * 5);
            uu = u_loc * 5 + u_bias;
            vv = v_loc * 5 + v_bias;
            fprintf(fid, '%f %f\n',uu,vv);
        end
    end
fclose(fid);
end
end

if 0
[A, Rm] = fun_para_readtxt();
end

if 0
parpool(36);
COUNT = 0;
%para = load('../result/train_num_28/para_temp_108.txt');
para=[Rm,A(1,1),A(1,3),A(2,2),A(2,3)];%�Ż��ڲ������������ƽ�ƾ��������ת�ǣ�
%options = optimset('Algorithm','levenberg-marquardt','InitDamping',1e2,'Display','iter');
options = optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt','Display','iter');

[x,nor,res,exitflag,output] = lsqnonlin( @fun_mapopt, para, [],[],options);
fun_save_para(x, 13500);
A=[x(n*6+1) 0 x(n*6+2); 0 x(n*6+3) x(n*6+4); 0,0,1];
disp('map�����Ż������ڲ���Ϊ');
disp(A);
res;
exitflag;
output;
delete(gcp('nocreate'));
fun_test(['para_temp_108.txt']);
end
if 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           ������������ѱ�?�����Ż��������ڲ�����������ת�Ƕȣ�ƽ��������                     %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

para=[Rm,A(1,1),A(1,3),A(2,2),A(2,3)];%�Ż��ڲ������������ƽ�ƾ��������ת�ǣ�
options = optimset('Algorithm','levenberg-marquardt','Display','iter');

[x,nor,res] = lsqnonlin( @fun2, para, [],[],options, mp, M, n);
A=[x(n*6+1) 0 x(n*6+2); 0 x(n*6+3) x(n*6+4); 0,0,1];
disp('�����ѱ궨���Ż������ڲ���Ϊ');
disp(A);
sqrt(nor/(49*n))
end
