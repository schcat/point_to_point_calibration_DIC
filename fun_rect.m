function fun_rect

res_camera_w =1920;
res_camera_h = 1080;
n=20;

load('results/speckle_map_test.mat');

for k=1:1:n
    I = imread(['image/ori/train_2/cali_',num2str(k),'.png']);
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
    imwrite(I_correct/max(I_correct(:)), ['image/ori/train_2/correct_',num2str(k),'.png']);
end

