%´Ë³ÌÐòÓÃÀ´Éú³É·ÂÕæµÄÔ²µãÍ¼Ïñ£¬½Ç¶È±È½Ï´ó£¬ÊÇ30¶ÈºÍ20¶È
clear all;
% res_h = 2128;
% res_w = 2128;
n = 20;
res_h = 2000;
res_w = 2000;
%res_result = 1024; %¶¨Òå"Ïà»ú"·Ö±æÂÊ
res_camera_h = 1080; %¶¨Òå"Ïà»ú"·Ö±æÂÊ
res_camera_w = 1920; %¶¨Òå"Ïà»ú"·Ö±æÂÊ
value = imread('simulate_synthetic.png');

[A, Rm] = fun_para_readtxt();
para=[Rm,A(1,1),A(1,3),A(2,2),A(2,3)];
R=[];
for i=1:n               %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ö¸ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
      R_new=para([(i-1)*6+1 : (i-1)*6+6] );
%      RL = vision.internal.calibration.rodriguesVectorToMatrix(R_new([1:3])');
      Q1=R_new(1);
      Q2=R_new(2);
      Q3=R_new(3);
      TL=R_new([4:6])';    
      RL=[cos(Q2)*cos(Q3)    sin(Q1)*sin(Q2)*cos(Q3)-cos(Q1)*sin(Q3)    cos(Q1)*sin(Q2)*cos(Q3)+sin(Q1)*sin(Q3);
          cos(Q2)*sin(Q3)    sin(Q1)*sin(Q2)*sin(Q3)+cos(Q1)*cos(Q3)    cos(Q1)*sin(Q2)*sin(Q3)-sin(Q1)*cos(Q3); 
          -sin(Q2)           sin(Q1)*cos(Q2)                            cos(Q1)*cos(Q2)];
      RT=[RL(:,1:2) , TL];
      R=[R;RT];
end
A=[para(n*6+1) 0 para(n*6+2); 0 para(n*6+3) para(n*6+4); 0,0,1];
RT = zeros(3,3,20);
for k=1:1:n
    RT(:,:,k)=R([(k-1)*3+1 : (k-1)*3+3],:);
end

dis = 0.03;
f = zeros(res_camera_h,res_camera_w);

for k = 5:6
camera_loc = zeros(res_h,res_w);
camera_loc_norm = zeros(res_h,res_w);
pix_loc = zeros(res_h,res_w);
result = zeros(res_camera_h,res_camera_w);

for x = 1:1:res_h
    for y = 1:1:res_w
        %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ïµï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ï?
        camera_loc(x,y,1) = RT(1,1,k)*x*dis + RT(1,2,k)*y*dis + RT(1,3,k)*1;
        camera_loc(x,y,2) = RT(2,1,k)*x*dis + RT(2,2,k)*y*dis + RT(2,3,k)*1;
        camera_loc(x,y,3) = RT(3,1,k)*x*dis + RT(3,2,k)*y*dis + RT(3,3,k)*1;  

        %ï¿½ï¿½Ò»ï¿½ï¿½ï¿½ï¿½Æ½ï¿½ï¿½(z=1)
        camera_loc_norm(x,y,1) = camera_loc(x,y,1)/camera_loc(x,y,3);
        camera_loc_norm(x,y,2) = camera_loc(x,y,2)/camera_loc(x,y,3);

        %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ïµ
        pix_loc(x,y,2) = A(1,1)*camera_loc_norm(x,y,1) + A(1,2)*camera_loc_norm(x,y,2) + A(1,3)*1;
        pix_loc(x,y,1) = A(2,1)*camera_loc_norm(x,y,1) + A(2,2)*camera_loc_norm(x,y,2) + A(2,3)*1;
    end
end
%ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿?
container_coor_1 = zeros(res_camera_h,res_camera_w,2);
container_coor_2 = zeros(res_camera_h,res_camera_w,2);
container_coor_3 = zeros(res_camera_h,res_camera_w,2);
container_coor_4 = zeros(res_camera_h,res_camera_w,2);
%ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Å·ï¿½Ï¾ï¿½ï¿½ï¿?
container_dist_1 = ones(res_camera_h,res_camera_w)*100;
container_dist_2 = ones(res_camera_h,res_camera_w)*100;
container_dist_3 = ones(res_camera_h,res_camera_w)*100;
container_dist_4 = ones(res_camera_h,res_camera_w)*100;
%ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ö?
container_value_1 = ones(res_camera_h,res_camera_w)*500;
container_value_2 = ones(res_camera_h,res_camera_w)*500;
container_value_3 = ones(res_camera_h,res_camera_w)*500;
container_value_4 = ones(res_camera_h,res_camera_w)*500;

for x = 2:1:res_h-1
    for y = 2:1:res_w-1
        xx=x;
        yy=y;
        if pix_loc(xx,yy,1)>1 && pix_loc(xx,yy,1)<res_camera_h && pix_loc(xx,yy,2)>1 && pix_loc(xx,yy,2)<res_camera_w
            %ï¿½ï¿½ï¿½ï¿½Ã¿ï¿½ï¿½Ô´ï¿½ï¿½ï¿½ï¿½Í¶Ó°ï¿½ï¿½Ä¿ï¿½ï¿½Í¼ï¿½ï¿½ï¿½Ïµï¿½Î»ï¿½Ã£ï¿½ï¿½ï¿½ï¿½ãµ½ï¿½ï¿½Î§ï¿½Ä¸ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ØµÄ¾ï¿½ï¿½ë£¬ÎªÃ¿ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ø´æ´¢ï¿½ï¿½ï¿½ï¿½Î§ï¿½Ä¸ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Úµï¿½ï¿½ï¿½Ð¡Å·ï¿½Ï¾ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ø£ï¿½
            %ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ò»ï¿½ï¿½Ã¿ï¿½ï¿½ï¿½ï¿½ï¿½Ø¶ï¿½ï¿½ï¿½ï¿½Ä¸ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½ï¿½Ë«ï¿½ï¿½ï¿½Ô²ï¿½Öµ

            %ï¿½ï¿½ï¿½ï¿½
            if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2<...
                    container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) && ...
                    pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && pix_loc(x,y,2)-floor(pix_loc(xx,yy,2))>0
                container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2;
                container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                container_value_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = value(x,y);
            end
            %ï¿½ï¿½ï¿½ï¿½
            if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
                    container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
                    pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
                container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
                container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                container_value_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
            end
            %ï¿½ï¿½ï¿½ï¿½
            if (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
                    container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
                    ceil(pix_loc(xx,yy,1))-pix_loc(x,y,1)>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
                container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
                container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
                container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
                container_value_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
            end
            %ï¿½ï¿½ï¿½ï¿½
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
            %Ë«ï¿½ï¿½ï¿½Ô²ï¿½Öµ
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
f = f + result;
end
imwrite(mat2gray(f), ['simulate_speckle_patent_proj.png']);

% 
% 
% cx = 512;
% cy = 512;
% f = 2500;
% par = [f,0,cx;
%        0,f,cy;
%        0,0,1];
% k1 = 0.25; k2 = -2.5; k3 = -5; k4 = -0.1; k5 = 0.25; k6 = 0.5;
% p1 = 0.001; p2 = 0.0012;
% s1 = -0.006; s2 = 0.005; s3 = -0.007; s4 = 0.006;
% 
% for num = 1:1:20
%     
%     camera_loc = zeros(res_h,res_w);
%     camera_loc_norm = zeros(res_h,res_w);
%     camera_loc_distort = zeros(res_h,res_w);
%     image_loc = zeros(res_h,res_w);
%     pix_loc = zeros(res_h,res_w);
%     result = ones(res_result_h,res_result_w)*255;
% 
%     rt = RT(:,:,num);
% 
%     for x = 1:1:res_h
%         for y = 1:1:res_w
%             %ÊÀ½ç×ø±êÏµµ½Ïà»ú×ø±êÏµ
%             camera_loc(x,y,1) = rt(1,1)*x*dis + rt(1,2)*y*dis + rt(1,3)*0 + rt(1,4)*1;
%             camera_loc(x,y,2) = rt(2,1)*x*dis + rt(2,2)*y*dis + rt(2,3)*0 + rt(2,4)*1;
%             camera_loc(x,y,3) = rt(3,1)*x*dis + rt(3,2)*y*dis + rt(3,3)*0 + rt(3,4)*1;  
%             
%             %¹éÒ»»¯ÏñÆ½Ãæ(z=1)
%             camera_loc_norm(x,y,1) = camera_loc(x,y,1)/camera_loc(x,y,3);
%             camera_loc_norm(x,y,2) = camera_loc(x,y,2)/camera_loc(x,y,3);
%             
%             %¼Ó»û±ä
%             r = camera_loc_norm(x,y,1)^2 + camera_loc_norm(x,y,2)^2;
%             %¾¶Ïò»û±ä
%             k_component = (1 + k1*r + k2*(r^2) + k3*(r^3)) / (1 + k4*r + k5*(r^2) + k6*(r^3));
%             k_component_x = camera_loc_norm(x,y,1) * k_component;
%             k_component_y = camera_loc_norm(x,y,2) * k_component;
%             %ÇÐÏò»û±ä
%             p_component_x = 2*p1*camera_loc_norm(x,y,1)*camera_loc_norm(x,y,2) + p2*(r + 2*camera_loc_norm(x,y,1)^2);
%             p_component_y = p1*(r + 2*camera_loc_norm(x,y,2)^2) + 2*p2*camera_loc_norm(x,y,1)*camera_loc_norm(x,y,2);
%             %±¡Àâ¾µ»û±ä
%             s_component_x = s1*r + s2*r^2;
%             s_component_y = s3*r + s4*r^2;
%             
%             camera_loc_distort(x,y,1) = k_component_x + p_component_x + s_component_x;
%             camera_loc_distort(x,y,2) = k_component_y + p_component_y + s_component_y;
%         
%             %µ½ÏñËØ×ø±êÏµ
%             pix_loc(x,y,1) = par(1,1)*camera_loc_distort(x,y,1) + par(1,2)*camera_loc_distort(x,y,2) + par(1,3)*1;
%             pix_loc(x,y,2) = par(2,1)*camera_loc_distort(x,y,1) + par(2,2)*camera_loc_distort(x,y,2) + par(2,3)*1;
%         end
%     end
%     %ÕûÊýÏñËØ×î½üÁÚÏñËØ×ø±ê
%     container_coor_1 = zeros(res_result_h,res_result_w,2);
%     container_coor_2 = zeros(res_result_h,res_result_w,2);
%     container_coor_3 = zeros(res_result_h,res_result_w,2);
%     container_coor_4 = zeros(res_result_h,res_result_w,2);
%     %ÕûÊýÏñËØ×î½üÁÚÏñËØÅ·ÊÏ¾àÀë
%     container_dist_1 = ones(res_result_h,res_result_w)*100;
%     container_dist_2 = ones(res_result_h,res_result_w)*100;
%     container_dist_3 = ones(res_result_h,res_result_w)*100;
%     container_dist_4 = ones(res_result_h,res_result_w)*100;
%     %ÕûÊýÏñËØ×î½üÁÚÏñËØÖµ
%     container_value_1 = ones(res_result_h,res_result_w)*500;
%     container_value_2 = ones(res_result_h,res_result_w)*500;
%     container_value_3 = ones(res_result_h,res_result_w)*500;
%     container_value_4 = ones(res_result_h,res_result_w)*500;
%     
%     for x = 2:1:res_h-1
%         for y = 2:1:res_w-1
%             xx=x;
%             yy=y;
%             if pix_loc(xx,yy,1)>1 && pix_loc(xx,yy,1)<res_result_h && pix_loc(xx,yy,2)>1 && pix_loc(xx,yy,2)<res_result_w
%                 %¶ÔÓÚÃ¿¸öÔ´ÏñËØÍ¶Ó°µ½Ä¿±êÍ¼ÏñÉÏµÄÎ»ÖÃ£¬¼ÆËãµ½ÖÜÎ§ËÄ¸öÕûÊýÏñËØµÄ¾àÀë£¬ÎªÃ¿¸öÕûÊýÏñËØ´æ´¢ÆäÖÜÎ§ËÄ¸öÇøÓòÄÚµÄ×îÐ¡Å·ÊÏ¾àÀëÏñËØ£¬
%                 %ÕâÑùÔÚÏÂÒ»²½Ã¿¸öÏñËØ¶¼ÓÐËÄ¸ö½üÁÚÏñËØÓÃÀ´×öË«ÏßÐÔ²åÖµ
%                 
%                 %×óÉÏ
%                 if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2<...
%                         container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) && ...
%                         pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && pix_loc(x,y,2)-floor(pix_loc(xx,yy,2))>0
%                     container_dist_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2;
%                     container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
%                     container_coor_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
%                     container_value_1(floor(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = value(x,y);
%                 end
%                 %ÓÒÉÏ
%                 if (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
%                         container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
%                         pix_loc(x,y,1)-floor(pix_loc(xx,yy,1))>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
%                     container_dist_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-floor(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
%                     container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
%                     container_coor_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
%                     container_value_2(floor(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
%                 end
%                 %ÓÒÏÂ
%                 if (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2<...
%                         container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) && ...
%                         ceil(pix_loc(xx,yy,1))-pix_loc(x,y,1)>0 && ceil(pix_loc(xx,yy,2))-pix_loc(x,y,2)>0
%                     container_dist_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-ceil(pix_loc(xx,yy,2)))^2;
%                     container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
%                     container_coor_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
%                     container_value_3(ceil(pix_loc(xx,yy,1)),ceil(pix_loc(xx,yy,2))) = value(x,y);
%                 end
%                 %×óÏÂ
%                 if (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2<...
%                         container_dist_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) && ...
%                         ceil(pix_loc(xx,yy,1))-pix_loc(x,y,1)>0 && pix_loc(x,y,2)-floor(pix_loc(xx,yy,2))>0
%                     container_dist_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = (pix_loc(x,y,1)-ceil(pix_loc(xx,yy,1)))^2+(pix_loc(x,y,2)-floor(pix_loc(xx,yy,2)))^2;
%                     container_coor_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),1) = pix_loc(x,y,1);
%                     container_coor_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2)),2) = pix_loc(x,y,2);
%                     container_value_4(ceil(pix_loc(xx,yy,1)),floor(pix_loc(xx,yy,2))) = value(x,y);
%                 end
%             end
%         end
%     end
%     
%     for x = 1:1:res_result_h
%         for y = 1:1:res_result_w
%             if container_value_1(x,y)<500 && container_value_2(x,y)<500 && container_value_3(x,y)<500 && container_value_4(x,y)<500
%                 %Ë«ÏßÐÔ²åÖµ
%                 result_1 = (container_coor_4(x,y,1) - x)/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_value_1(x,y) + ...
%                     (x - container_coor_1(x,y,1))/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_value_4(x,y);
%                 y_1 = (container_coor_4(x,y,1) - x)/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_coor_1(x,y,2) + ...
%                     (x - container_coor_1(x,y,1))/(container_coor_4(x,y,1) - container_coor_1(x,y,1)) * container_coor_4(x,y,2);
%                 result_2 = (container_coor_3(x,y,1) - x)/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_value_2(x,y) + ...
%                     (x - container_coor_2(x,y,1))/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_value_3(x,y);
%                 y_2 = (container_coor_3(x,y,1) - x)/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_coor_2(x,y,2) + ...
%                     (x - container_coor_2(x,y,1))/(container_coor_3(x,y,1) - container_coor_2(x,y,1)) * container_coor_3(x,y,2);
%                 result_value = (y_2 - y)/(y_2 - y_1) * result_1 + (y - y_1)/(y_2 - y_1) * result_2;
%                 result(x,y) = result_value;
%             end
%         end
%     end
%     imshow(result/255);
%     %imwrite(mat2gray(result), ['image/ori/simulate_speckle_patent_',num2str(num),'_0417.png']);
%     imwrite(mat2gray(result), ['simulate_speckle_patent_',num2str(num),'_0417.png']);
% end
% 
