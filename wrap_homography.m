function result = wrap_homography(value, h, res_h, res_w, res_result_h, res_result_w)
    result = zeros(res_h,res_w);
    pix_loc = zeros(res_h,res_w,2);
    for x = 1:1:res_h
        for y = 1:1:res_w
            pix_loc(x,y,1) = (h(1,1)*x + h(1,2)*y + h(1,3)) / (h(3,1)*x + h(3,2)*y + h(3,3));
            pix_loc(x,y,2) = (h(2,1)*x + h(2,2)*y + h(2,3)) / (h(3,1)*x + h(3,2)*y + h(3,3));
        end
    end
    %整数像素最近邻像素坐标
    container_coor_1 = zeros(res_result_h,res_result_w,2);
    container_coor_2 = zeros(res_result_h,res_result_w,2);
    container_coor_3 = zeros(res_result_h,res_result_w,2);
    container_coor_4 = zeros(res_result_h,res_result_w,2);
    %整数像素最近邻像素欧氏距离
    container_dist_1 = ones(res_result_h,res_result_w)*100;
    container_dist_2 = ones(res_result_h,res_result_w)*100;
    container_dist_3 = ones(res_result_h,res_result_w)*100;
    container_dist_4 = ones(res_result_h,res_result_w)*100;
    %整数像素最近邻像素值
    container_value_1 = ones(res_result_h,res_result_w)*500;
    container_value_2 = ones(res_result_h,res_result_w)*500;
    container_value_3 = ones(res_result_h,res_result_w)*500;
    container_value_4 = ones(res_result_h,res_result_w)*500;
    
    for x = 2:1:res_h-1
        for y = 2:1:res_w-1
            xx=x;
            yy=y;
            if pix_loc(xx,yy,1)>1 && pix_loc(xx,yy,1)<res_result_h && pix_loc(xx,yy,2)>1 && pix_loc(xx,yy,2)<res_result_w
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
end