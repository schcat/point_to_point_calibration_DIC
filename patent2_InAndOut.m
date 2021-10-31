function  patent2_InAndOut
%      ��m�ļ���������ڲ�������������Ż���
%      �����ƹ�������ͼ��ÿ��ͼ��ȡ10���㡣�����ÿ��ͼ��ĵ�Ӧ����Ȼ����������?�������һ��?
%      ��V*b=0���?b����Ȼ����b�ֽ���ڲ�������?A�������?A�͵�Ӧ����������������
global COUNT;
num = 70;
n = 36;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 �Ա�׼ͼƬ-����ͼƬ����DIC���?                                          %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if 0
height_ref = 1000;
width_ref = 1000;
pad = 100;
pathname = '../image/';
filename_ref = 'speckle_pattern_4000_pad_0111_20_15000_1000.png';
filename_roi = 'ROI_1000_1000_100.png';
radius = num; %DIC subset radius
spacing = 0;
% ���?ROI��pad�б仯���ǵ���ncorr_auto_initseeds.m�����?pad

mat_name = ['../result/train_2_14/speckle_'];
tic
parpool(32);
parfor k=1:1:n
    filename_cur = ['../image/train_2_14/cali_',num2str(k),'.png'];
    displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);
%    save(['results/ori_ref_mat/speckle_',num2str(k),'.mat'],'displacements');
end
delete(gcp('nocreate'));
toc
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 ��DIC���������ȡ���Ƶ�?                                              %                           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
if 0
fid = fopen(['../result/train_2_14/speckle_ori_',num2str(num),'.txt'],'w');
for k = 1:1:n
    str=strcat('../result/train_2_14/speckle_',num2str(k),'.mat');
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
end
fclose(fid);
end

if 1
[A, Rm] = fun_para_readtxt();
end

if 1
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
