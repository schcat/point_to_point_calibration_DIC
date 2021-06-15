%�˳����������ɷ����Բ��ͼ�񣬽ǶȱȽϴ���30�Ⱥ�20��
function f=fun_projection(value, A, RT, k)

[res_h,res_w] = size(value);
res_camera_w = 1920;
res_camera_h = 1080;
dis = 0.03; %60mm/2000pix
    
camera_loc = zeros(res_h,res_w);
camera_loc_norm = zeros(res_h,res_w);
pix_loc = zeros(res_h,res_w);
result = zeros(res_camera_h,res_camera_w);

for x = 1:1:res_h
    for y = 1:1:res_w
        %��������ϵ���������ϵ
        camera_loc(x,y,1) = RT(1,1)*x*dis + RT(1,2)*y*dis + RT(1,3)*1;
        camera_loc(x,y,2) = RT(2,1)*x*dis + RT(2,2)*y*dis + RT(2,3)*1;
        camera_loc(x,y,3) = RT(3,1)*x*dis + RT(3,2)*y*dis + RT(3,3)*1;  

        %��һ����ƽ��(z=1)
        camera_loc_norm(x,y,1) = camera_loc(x,y,1)/camera_loc(x,y,3);
        camera_loc_norm(x,y,2) = camera_loc(x,y,2)/camera_loc(x,y,3);

        %����������ϵ
        pix_loc(x,y,2) = A(1,1)*camera_loc_norm(x,y,1) + A(1,2)*camera_loc_norm(x,y,2) + A(1,3)*1;
        pix_loc(x,y,1) = A(2,1)*camera_loc_norm(x,y,1) + A(2,2)*camera_loc_norm(x,y,2) + A(2,3)*1;
    end
end
%���������������������
container_coor_1 = zeros(res_camera_h,res_camera_w,2);
container_coor_2 = zeros(res_camera_h,res_camera_w,2);
container_coor_3 = zeros(res_camera_h,res_camera_w,2);
container_coor_4 = zeros(res_camera_h,res_camera_w,2);
%�����������������ŷ�Ͼ���
container_dist_1 = ones(res_camera_h,res_camera_w)*100;
container_dist_2 = ones(res_camera_h,res_camera_w)*100;
container_dist_3 = ones(res_camera_h,res_camera_w)*100;
container_dist_4 = ones(res_camera_h,res_camera_w)*100;
%�����������������ֵ
container_value_1 = ones(res_camera_h,res_camera_w)*500;
container_value_2 = ones(res_camera_h,res_camera_w)*500;
container_value_3 = ones(res_camera_h,res_camera_w)*500;
container_value_4 = ones(res_camera_h,res_camera_w)*500;

for x = 2:1:res_h-1
    for y = 2:1:res_w-1
        xx=x;
        yy=y;
        if pix_loc(xx,yy,1)>1 && pix_loc(xx,yy,1)<res_camera_h && pix_loc(xx,yy,2)>1 && pix_loc(xx,yy,2)<res_camera_w
            %����ÿ��Դ����ͶӰ��Ŀ��ͼ���ϵ�λ�ã����㵽��Χ�ĸ��������صľ��룬Ϊÿ���������ش洢����Χ�ĸ������ڵ���Сŷ�Ͼ������أ�
            %��������һ��ÿ�����ض����ĸ���������������˫���Բ�ֵ

            %����
            if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2<...
                    container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) && ...
                    pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && pix_loc(x,y,2)-floor(pix_loc(xx,yy,2))>0
                container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2;
                container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                container_value_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = value(x,y);
            end
            %����
            if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
                    container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
                    pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
                container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
                container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                container_value_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
            end
            %����
            if (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
                    container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
                    ceil(pix_loc(xx,yy,1))-pix_loc(x,y,1)>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
                container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
                container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                container_value_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
            end
            %����
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

for x = 1:1:res_camera_h
    for y = 1:1:res_camera_w
        if container_value_1(x,y)<500 && container_value_2(x,y)<500 && container_value_3(x,y)<500 && container_value_4(x,y)<500
            %˫���Բ�ֵ
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
f = result;

