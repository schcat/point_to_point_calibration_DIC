res_camera_w =1920;
res_camera_h = 1080;
n=1428;
point = zeros(n*4,2);
load('../result/train_1_1/speckle_map_test_1.mat');
load('../result/train_2_1/speckle_map_test_2.mat');
point_file = load('~/cnn/implementation/point_distance/build/speckle_test_ori.txt');

for i=1:n
    point(i*4-3,:) = point_file(i,1:2);
    point(i*4-2,:) = point_file(i,3:4);
    point(i*4-1,:) = point_file(i,5:6);
    point(i*4,:) = point_file(i,7:8);
end

%整数像素最近邻像素坐标
container_coor_1 = zeros(n*4,2);
container_coor_2 = zeros(n*4,2);
container_coor_3 = zeros(n*4,2);
container_coor_4 = zeros(n*4,2);
%整数像素最近邻像素欧氏距离
container_dist_1 = ones(n*4)*100;
container_dist_2 = ones(n*4)*100;
container_dist_3 = ones(n*4)*100;
container_dist_4 = ones(n*4)*100;
%整数像素最近邻像素值
container_ori_1 = ones(n*4,2)*500;
container_ori_2 = ones(n*4,2)*500;
container_ori_3 = ones(n*4,2)*500;
container_ori_4 = ones(n*4,2)*500;

point_ori = zeros(n*4,2);

for i = 1:1:res_camera_h
    for j = 1:1:res_camera_w

        u_bias_1 = mean_map_u_1(i,j);
        v_bias_1 = mean_map_v_1(i,j);

		if u_bias_1 == 0 || v_bias_1 == 0 || j + u_bias_1 < 1 || i + v_bias_1 < 1 || j + u_bias_1 > 1920 || i + v_bias_1 > 1080
			uu_1 = j;
			vv_1 = i;
        else 
			uu_1 = j + u_bias_1;
			vv_1 = i + v_bias_1;
		end
        u_bias_2 = mean_map_u_2(i,j);
        v_bias_2 = mean_map_v_2(i,j);

		if u_bias_2 == 0 || v_bias_2 == 0 || j + u_bias_2 < 1 || i + v_bias_2 < 1 || j + u_bias_2 > 1920 || i + v_bias_2 > 1080
			uu_2 = j;
			vv_2 = i;
        else 
			uu_2 = j + u_bias_2;
			vv_2 = i + v_bias_2;
		end
		for k=1:1:n*4
			if mod(k,4) == 1 || mod(k,4) == 2
				dist =  sqrt((point(k,1) - uu_1)^2 + (point(k,2) - vv_1)^2);
				if dist < 5
					if point(k,1) > uu_1 && point(k,2) > vv_1 % left top
						if dist < container_dist_1(k)
							container_coor_1(k,1) = uu_1;
							container_coor_1(k,2) = vv_1;
							container_ori_1(k,1) = j;
							container_ori_1(k,2) = i;
							container_dist_1(k) = dist;
						end
					elseif point(k,1) < uu_1 && point(k,2) > vv_1 % right top
						if dist < container_dist_2(k)
							container_coor_2(k,1) = uu_1;
							container_coor_2(k,2) = vv_1;
							container_ori_2(k,1) = j;
							container_ori_2(k,2) = i;
							container_dist_2(k) = dist;
						end
					elseif point(k,1) < uu_1 && point(k,2) < vv_1 % right botton
						if dist < container_dist_3(k)
							container_coor_3(k,1) = uu_1;
							container_coor_3(k,2) = vv_1;
							container_ori_3(k,1) = j;
							container_ori_3(k,2) = i;
							container_dist_3(k) = dist;
						end
					elseif point(k,1) > uu_1 && point(k,2) < vv_1 % left botton
						if dist < container_dist_4(k)
							container_coor_4(k,1) = uu_1;
							container_coor_4(k,2) = vv_1;
							container_ori_4(k,1) = j;
							container_ori_4(k,2) = i;
							container_dist_4(k) = dist;
						end
					end
				end
			elseif mod(k,4) == 3 || mod(k,4) == 0
				dist =  sqrt((point(k,1) - uu_2)^2 + (point(k,2) - vv_2)^2);
				if dist < 5
					if point(k,1) > uu_2 && point(k,2) > vv_2 % left top
						if dist < container_dist_1(k)
							container_coor_1(k,1) = uu_2;
							container_coor_1(k,2) = vv_2;
							container_ori_1(k,1) = j;
							container_ori_1(k,2) = i;
							container_dist_1(k) = dist;
						end
					elseif point(k,1) < uu_2 && point(k,2) > vv_2 % right top
						if dist < container_dist_2(k)
							container_coor_2(k,1) = uu_2;
							container_coor_2(k,2) = vv_2;
							container_ori_2(k,1) = j;
							container_ori_2(k,2) = i;
							container_dist_2(k) = dist;
						end
					elseif point(k,1) < uu_2 && point(k,2) < vv_2 % right botton
						if dist < container_dist_3(k)
							container_coor_3(k,1) = uu_2;
							container_coor_3(k,2) = vv_2;
							container_ori_3(k,1) = j;
							container_ori_3(k,2) = i;
							container_dist_3(k) = dist;
						end
					elseif point(k,1) > uu_2 && point(k,2) < vv_2 % left botton
						if dist < container_dist_4(k)
							container_coor_4(k,1) = uu_2;
							container_coor_4(k,2) = vv_2;
							container_ori_4(k,1) = j;
							container_ori_4(k,2) = i;
							container_dist_4(k) = dist;
						end
					end
				end
			end
		end

	end
end
for k=1:1:n*4
	if container_dist_1(k)<100 && container_dist_2(k)<100 && container_dist_3(k)<100 && container_dist_4(k)<100
		% two bilinear interpolation
		result_1 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(k,1) - container_coor_4(k,1)) * container_ori_4(k,1) + ...
			(point(k,1) - container_coor_4(k,1))/(container_coor_3(k,1) - container_coor_4(k,1)) * container_ori_3(k,1);
		y_2 = (container_coor_3(k,2) - point(k,2))/(container_coor_3(k,2) - container_coor_4(k,2)) * container_ori_4(k,2) + ...
			(point(k,2) - container_coor_4(k,2))/(container_coor_3(k,2) - container_coor_4(k,2)) * container_ori_3(k,2);
		result_2 = (container_coor_2(k,1) - point(k,1))/(container_coor_2(k,1) - container_coor_1(k,1)) * container_ori_1(k,1) + ...
			(point(k,1) - container_coor_1(k,1))/(container_coor_2(k,1) - container_coor_1(k,1)) * container_ori_2(k,1);
		y_1 = (container_coor_2(k,2) - point(k,2))/(container_coor_2(k,2) - container_coor_1(k,2)) * container_ori_1(k,2) + ...
			(point(k,2) - container_coor_1(k,2))/(container_coor_2(k,2) - container_coor_1(k,2)) * container_ori_2(k,2);
		result_value = (y_2 - point(k,2))/(y_2 - y_1) * result_1 + (point(k,2) - y_1)/(y_2 - y_1) * result_2;
		point_ori(k,1) = result_value;

		result_1 = (container_coor_4(k,2) - point(k,2))/(container_coor_4(k,2) - container_coor_1(k,2)) * container_ori_1(k,2) + ...
			(point(k,2) - container_coor_1(k,2))/(container_coor_4(k,2) - container_coor_1(k,2)) * container_ori_4(k,2);
		y_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_1(k,1) + ...
			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_4(k,1);
		result_2 = (container_coor_3(k,2) - point(k,2))/(container_coor_3(k,2) - container_coor_2(k,2)) * container_ori_2(k,2) + ...
			(point(k,2) - container_coor_2(k,2))/(container_coor_3(k,2) - container_coor_2(k,2)) * container_ori_3(k,2);
		y_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_2(k,1) + ...
			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_3(k,1);
		result_value = (y_2 - point(k,1))/(y_2 - y_1) * result_1 + (point(k,1) - y_1)/(y_2 - y_1) * result_2;
		point_ori(k,2) = result_value;

%		result_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_1(k,2) + ...
%			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_4(k,2);
%		y_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_1(k,2) + ...
%			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_4(k,2);
%		result_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_2(k,2) + ...
%			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_3(k,2);
%		y_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_2(k,2) + ...
%			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_3(k,2);
%		result_value = (y_2 - point(k,2))/(y_2 - y_1) * result_1 + (point(k,2) - y_1)/(y_2 - y_1) * result_2;
%		point_ori(k,2) = result_value;
	end
end

fid = fopen('~/cnn/implementation/point_distance/build/speckle_test_speckle_novel_undistort.txt','w');
for k = 1:1:n*4
	if mod(k,4) == 0
		fprintf(fid, '%f %f\n',point_ori(k,1),point_ori(k,2));
	else
		fprintf(fid, '%f %f ',point_ori(k,1),point_ori(k,2));
	end
end
fclose(fid);
