res_camera_w =1920;
res_camera_h = 1080;
n=20;
point = zeros(n*4,2);
load('../result/train_2_2_70/speckle_map_test_70.mat');
point_file = load('~/cnn/implementation/point_distance/build/speckle_points_97.txt');

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

        u_bias = mean_map_u(i,j);
        v_bias = mean_map_v(i,j);

		if u_bias == 0 || v_bias == 0 || j + u_bias < 1 || i + v_bias < 1 || j + u_bias > 1920 || i + v_bias > 1080
			uu = j;
			vv = i;
        else 
			uu = j + u_bias;
			vv = i + v_bias;
                
%			value_1 = I(floor(vv),floor(uu));
%			value_2 = I(floor(vv),ceil(uu));
%			value_3 = I(ceil(vv),ceil(uu));
%			value_4 = I(ceil(vv),floor(uu));
%			result_1 = (ceil(vv) - vv) * value_1 + (vv - floor(vv)) * value_4;
%			result_2 = (ceil(vv) - vv) * value_2 + (vv - floor(vv)) * value_3;
%			result_value = (ceil(uu) - uu) * result_1 + (uu - floor(uu)) * result_2;

		for k=1:1:n*4
			dist =  sqrt((point(k,1) - uu)^2 + (point(k,2) - vv)^2);
			if dist < 3
				if point(k,1) > uu && point(k,2) > vv % left top
					if dist < container_dist_1(k)
						container_coor_1(k,1) = uu;
						container_coor_1(k,2) = vv;
						container_ori_1(k,1) = j;
						container_ori_1(k,2) = i;
						container_dist_1(k) = dist;
					end
				elseif point(k,1) < uu && point(k,2) > vv % right top
					if dist < container_dist_2(k)
						container_coor_2(k,1) = uu;
						container_coor_2(k,2) = vv;
						container_ori_2(k,1) = j;
						container_ori_2(k,2) = i;
						container_dist_2(k) = dist;
					end
				elseif point(k,1) < uu && point(k,2) < vv % right botton
					if dist < container_dist_3(k)
						container_coor_3(k,1) = uu;
						container_coor_3(k,2) = vv;
						container_ori_3(k,1) = j;
						container_ori_3(k,2) = i;
						container_dist_3(k) = dist;
					end
				elseif point(k,1) > uu && point(k,2) < vv % left botton
					if dist < container_dist_4(k)
						container_coor_4(k,1) = uu;
						container_coor_4(k,2) = vv;
						container_ori_4(k,1) = j;
						container_ori_4(k,2) = i;
						container_dist_4(k) = dist;
					end
				end
			end
		end

            
	end
end
for k=1:1:n*4
	if container_dist_1(k)<100 && container_dist_2(k)<100 && container_dist_3(k)<100 && container_dist_4(k)<100
		% two bilinear interpolation
		result_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_1(k,1) + ...
			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_4(k,1);
		y_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_coor_1(k,2) + ...
			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_coor_4(k,2);
		result_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(x,y,1) - container_coor_2(k,1)) * container_ori_2(k,1) + ...
			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_3(k,1);
		y_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_coor_2(k,2) + ...
			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_coor_3(k,2);
		result_value = (y_2 - point(k,2))/(y_2 - y_1) * result_1 + (point(k,2) - y_1)/(y_2 - y_1) * result_2;
		point_ori(k,1) = result_value;

		result_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_1(k,2) + ...
			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_ori_4(k,2);
		y_1 = (container_coor_4(k,1) - point(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_coor_1(k,2) + ...
			(point(k,1) - container_coor_1(k,1))/(container_coor_4(k,1) - container_coor_1(k,1)) * container_coor_4(k,2);
		result_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(x,y,1) - container_coor_2(k,1)) * container_ori_2(k,2) + ...
			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_ori_3(k,2);
		y_2 = (container_coor_3(k,1) - point(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_coor_2(k,2) + ...
			(point(k,1) - container_coor_2(k,1))/(container_coor_3(k,1) - container_coor_2(k,1)) * container_coor_3(k,2);
		result_value = (y_2 - point(k,2))/(y_2 - y_1) * result_1 + (point(k,2) - y_1)/(y_2 - y_1) * result_2;
		point_ori(k,2) = result_value;
	end
end

fid = fopen('~/cnn/implementation/point_distance/build/speckle_points_ori_97.txt','w');
for k = 1:1:n*4
	if mod(k,4) == 0
		fprintf(fid, '%f %f\n',point_ori(k,1),point_ori(k,2));
	else
		fprintf(fid, '%f %f ',point_ori(k,1),point_ori(k,2));
	end
end
fclose(fid);
