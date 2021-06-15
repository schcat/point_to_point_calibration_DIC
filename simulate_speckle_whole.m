%此程序用来生成仿真的圆点图像，角度比较大，是30度和20度
clear all;
res_h = 2128;
res_w = 2128;
res_result = 1024; %定义"相机"分辨率
%value = imread('simulate/speckle_pattern_4096_pad_0918.png');
value = imread('image/speckle_pattern_2048_40_pad_8_350000.tif');

dis = 0.058027;
RT = zeros(3,4,24);
A = zeros(3,3,24);

B = [-62;-62;270];

A(:,:,1) = [0.9986295348013184, 0, 0.05233595535100499;
            0, 1, 0;
            -0.05233595535100499, 0, 0.9986295348013184];
A(:,:,2) = [1, 0, 0;
            0, 0.9986295348013184, -0.05233595535100499;
            0, 0.05233595535100499, 0.9986295348013184];
A(:,:,3) = [0.9986298478876809, 0.001370152112319112, 0.052312040593122;
            0.001370152112319112, 0.9986298478876809, -0.052312040593122;
            -0.052312040593122, 0.052312040593122, 0.9972596957753618];
A(:,:,4) = [0.9975640503428961, 0, -0.06975647255614205;
            0, 1, 0;
            0.06975647255614205, 0, 0.9975640503428961];
A(:,:,5) = [1, 0, 0;
            0, 0.9975640503428961, 0.06975647255614205;
            0, -0.06975647255614205, 0.9975640503428961];
A(:,:,6) = [0.9975650396394263, 0.002434960360573612, -0.06969980385361313;
            0.002434960360573612, 0.9975650396394263, 0.06969980385361313;
            0.06969980385361313, -0.06969980385361313, 0.9951300792788528];
A(:,:,7) = [0.9961946982214861, 0, 0.08715574126471727;
            0, 1, 0;
            -0.08715574126471727, 0, 0.9961946982214861];
A(:,:,8) = [1, 0, 0;
            0, 0.9961946982214861, -0.08715574126471727;
            0, 0.08715574126471727, 0.9961946982214861];
A(:,:,9) = [0.9961971128340786, 0.003802887165921364, -0.08704510572254609;
            0.003802887165921364, 0.9961971128340786, 0.08704510572254609;
            0.08704510572254609, -0.08704510572254609, 0.9923942256681573];
A(:,:,10) = [0.9945218955549954, 0, -0.1045284614911124;
            0, 1, 0;
            0.1045284614911124, 0, 0.9945218955549954];
A(:,:,11) = [1, 0, 0;
            0, 0.9945218955549954, 0.1045284614911124;
            0, -0.1045284614911124, 0.9945218955549954];
A(:,:,12) = [0.9945269008179672, 0.005473099182032758, -0.1043373793745694;
            0.005473099182032758, 0.9945269008179672, 0.1043373793745694;
            0.1043373793745694, -0.1043373793745694, 0.9890538016359345];
A(:,:,13) = [0.9925461518953035, 0, 0.1218693413366346;
            0, 1, 0;
            -0.1218693413366346, 0, 0.9925461518953035];
A(:,:,14) = [1, 0, 0;
            0, 0.9925461518953035, -0.1218693413366346;
            0, 0.1218693413366346, 0.9925461518953035];
A(:,:,15) = [0.9925554210907472, 0.007444578909252706, 0.1215660902893287;
            0.007444578909252706, 0.9925554210907472, -0.1215660902893287;
            -0.1215660902893287, 0.1215660902893287, 0.9851108421814946];
A(:,:,16) = [0.9902680690730484, 0, -0.1391730986014763;
            0, 1, 0;
            0.1391730986014763, 0, 0.9902680690730484];
A(:,:,17) = [1, 0, 0;
            0, 0.9902680690730484, 0.1391730986014763;
            0, -0.1391730986014763, 0.9902680690730484];
A(:,:,18) = [0.9902838746855357, 0.009716125314464306, -0.1387207426691331;
            0.009716125314464306, 0.9902838746855357, 0.1387207426691331;
            0.1387207426691331, -0.1387207426691331, 0.9805677493710714];
A(:,:,19) = [0.9876883410143023, 0, 0.1564344623937302;
            0, 1, 0;
            -0.1564344623937302, 0, 0.9876883410143023];
A(:,:,20) = [1, 0, 0;
            0, 0.9876883410143023, -0.1564344623937302;
            0, 0.1564344623937302, 0.9876883410143023];
A(:,:,21) = [0.9877136454372361, 0.01228635456276394, 0.1557908858330423;
            0.01228635456276394, 0.9877136454372361, -0.1557908858330423;
            -0.1557908858330423, 0.1557908858330423, 0.9754272908744721];
A(:,:,22) = [0.9848077535291954, 0, -0.1736481747349499;
            0, 1, 0;
            0.1736481747349499, 0, 0.9848077535291954];
A(:,:,23) = [1, 0, 0;
            0, 0.9848077535291954, 0.1736481747349499;
            0, -0.1736481747349499, 0.9848077535291954];
A(:,:,24) = [0.9848462991395015, 0.01515370086049844, -0.1727661205834581;
            0.01515370086049844, 0.9848462991395015, 0.1727661205834581;
            0.1727661205834581, -0.1727661205834581, 0.9696925982790031];

for num = 1:1:24
    RT(:,:,num) = cat(2,A(:,:,num),B);
end

cx = 512;
cy = 512;
f = 2500;
par = [f,0,cx;
       0,f,cy;
       0,0,1];
%k1 = 0.5; k2 = -5; k3 = -10; k4 = -0.1; k5 = 0.5; k6 = 1;
%p1 = 0.001; p2 = 0.0012;
%s1 = -0.006; s2 = 0.005; s3 = -0.007; s4 = 0.006;

k1 = 0; k2 = 0; k3 = 0; k4 = 0; k5 = 0; k6 = 0;
p1 = 0; p2 = 0;
s1 = 0; s2 = 0; s3 = 0; s4 = 0;

%k1 = 0.5; k2 = -8; k3 = -16; k4 = -0.1; k5 = 4; k6 = 8;
%p1 = 0.04; p2 = 0.03;
%s1 = -0.08; s2 = 0.05; s3 = -0.07; s4 = 0.06;

for num = 1:1:24
    
    camera_loc = zeros(res_h,res_w);
    camera_loc_norm = zeros(res_h,res_w);
    camera_loc_distort = zeros(res_h,res_w);
    pix_loc = zeros(res_h,res_w,2);
    result = ones(res_result)*255;

    rt = RT(:,:,num);

    for x = 1:1:res_h
        for y = 1:1:res_w
            %世界坐标系到相机坐标系
            camera_loc(x,y,1) = rt(1,1)*x*dis + rt(1,2)*y*dis + rt(1,3)*0 + rt(1,4)*1;
            camera_loc(x,y,2) = rt(2,1)*x*dis + rt(2,2)*y*dis + rt(2,3)*0 + rt(2,4)*1;
            camera_loc(x,y,3) = rt(3,1)*x*dis + rt(3,2)*y*dis + rt(3,3)*0 + rt(3,4)*1;  
            
            %归一化像平面(z=1)
            camera_loc_norm(x,y,1) = camera_loc(x,y,1)/camera_loc(x,y,3);
            camera_loc_norm(x,y,2) = camera_loc(x,y,2)/camera_loc(x,y,3);
            
            %加畸变
            r = camera_loc_norm(x,y,1)^2 + camera_loc_norm(x,y,2)^2;
            %径向畸变
            k_component = (1 + k1*r + k2*(r^2) + k3*(r^3)) / (1 + k4*r + k5*(r^2) + k6*(r^3));
            k_component_x = camera_loc_norm(x,y,1) * k_component;
            k_component_y = camera_loc_norm(x,y,2) * k_component;
            %切向畸变
            p_component_x = 2*p1*camera_loc_norm(x,y,1)*camera_loc_norm(x,y,2) + p2*(r + 2*camera_loc_norm(x,y,1)^2);
            p_component_y = p1*(r + 2*camera_loc_norm(x,y,2)^2) + 2*p2*camera_loc_norm(x,y,1)*camera_loc_norm(x,y,2);
            %薄棱镜畸变
            s_component_x = s1*r + s2*r^2;
            s_component_y = s3*r + s4*r^2;
            
            camera_loc_distort(x,y,1) = k_component_x + p_component_x + s_component_x;
            camera_loc_distort(x,y,2) = k_component_y + p_component_y + s_component_y;
        
            %到像素坐标系
            pix_loc(x,y,1) = par(1,1)*camera_loc_distort(x,y,1) + par(1,2)*camera_loc_distort(x,y,2) + par(1,3)*1;
            pix_loc(x,y,2) = par(2,1)*camera_loc_distort(x,y,1) + par(2,2)*camera_loc_distort(x,y,2) + par(2,3)*1;
        end
    end
    %整数像素最近邻像素坐标
    container_coor_1 = zeros(res_result,res_result,2);
    container_coor_2 = zeros(res_result,res_result,2);
    container_coor_3 = zeros(res_result,res_result,2);
    container_coor_4 = zeros(res_result,res_result,2);
    %整数像素最近邻像素欧氏距离
    container_dist_1 = ones(res_result,res_result)*100;
    container_dist_2 = ones(res_result,res_result)*100;
    container_dist_3 = ones(res_result,res_result)*100;
    container_dist_4 = ones(res_result,res_result)*100;
    %整数像素最近邻像素值
    container_value_1 = ones(res_result,res_result)*500;
    container_value_2 = ones(res_result,res_result)*500;
    container_value_3 = ones(res_result,res_result)*500;
    container_value_4 = ones(res_result,res_result)*500;
    
    for x = 2:1:res_h-1
        for y = 2:1:res_w-1
            xx=x;
            yy=y;
            if pix_loc(xx,yy,1)>1 && pix_loc(xx,yy,1)<res_result && pix_loc(xx,yy,2)>1 && pix_loc(xx,yy,2)<res_result
                %对于每个源像素投影到目标图像上的位置，计算到周围四个整数像素的距离，为每个整数像素存储其周围四个区域内的最小欧氏距离像素，
                %这样在下一步每个像素都有四个近邻像素用来做双线性插值
                
                %左上
                if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2<...
                        container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) && ...
                        pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && pix_loc(x,y,2)-floor(pix_loc(xx,yy,2))>0
                    container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2;
                    container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                    container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                    container_value_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = value(x,y);
                end
                %右上
                if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
                        container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
                        pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
                    container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
                    container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                    container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                    container_value_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
                end
                %右下
                if (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
                        container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
                        ceil(pix_loc(xx,yy,1))-pix_loc(x,y,1)>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
                    container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
                    container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                    container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                    container_value_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
                end
                %左下
                if (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2<...
                        container_dist_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) && ...
                        ceil(pix_loc(xx,yy,1))-pix_loc(x,y,1)>0 && pix_loc(x,y,2)-floor(pix_loc(xx,yy,2))>0
                    container_dist_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2;
                    container_coor_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                    container_coor_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                    container_value_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = value(x,y);
                end
            end
        end
    end
    
    for x = 1:1:res_result
        for y = 1:1:res_result
            if container_value_1(x,y)<500 && container_value_2(x,y)<500 && container_value_3(x,y)<500 && container_value_4(x,y)<500
                %双线性插值
                result_1 = (container_coor_4(x,y,1) - x)/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_value_1(x,y) + ...
                    (x - container_coor_1(x,y,1))/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_value_4(x,y);
                y_1 = (container_coor_4(x,y,1) - x)/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_coor_1(x,y,2) + ...
                    (x - container_coor_1(x,y,1))/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_coor_4(x,y,2);
                result_2 = (container_coor_3(x,y,1) - x)/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_value_2(x,y) + ...
                    (x - container_coor_2(x,y,1))/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_value_3(x,y);
                y_2 = (container_coor_3(x,y,1) - x)/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_coor_2(x,y,2) + ...
                    (x - container_coor_2(x,y,1))/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_coor_3(x,y,2);
                result_value = (y_2 - y)/(y_2 - y_1) * result_1 + (y - y_1)/(y_2 - y_1) * result_2;
                result(x,y) = result_value;
            end
        end
    end
    imshow(result/255);
    imwrite(mat2gray(result), ['image/ori_whole/simulate_speckle_whole_',num2str(num),'.tif']);
end

