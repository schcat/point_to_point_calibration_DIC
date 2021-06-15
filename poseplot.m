% R_arm = load('para_r_arm.txt');
% T_arm = load('para_t_arm.txt');
% R_cam = load('para_r.txt');
% T_cam = load('para_t.txt');
% R_dic = load('para_r_dic.txt');
% T_dic = load('para_t_dic.txt');
% 
% RT_arm = zeros(4,4,20);
% for i = 1:1:20
%     RT_arm(:,:,i) = [Rodrigues(R_arm(i,:)),T_arm(i,:)';0,0,0,1];
% end
% 
% inv(Rodrigues(R_arm(1,:)))*[0,0,1]';
% 
% RT_cam = zeros(4,4,20);
% for i = 1:1:20
%     RT_cam(:,:,i) = [Rodrigues(R_cam(i,:)),T_cam(i,:)'/1000;0,0,0,1];
% end
% 
% RT_dic = zeros(4,4,20);
% for i = 1:1:20
%     RT_dic(:,:,i) = [Rodrigues(R_dic(i,:)),T_dic(i,:)'/1000;0,0,0,1];
% end
% 
% vec_cam_ori = [0,0,0.01];
% 
% scatter3(0,0,0,20,'r')
% hold on
% quiver3(0,0,0,0,0,0.01,'-b','LineWidth',2)
% hold on
% quiver3(0,0,0,0,0.01,0,'-r','LineWidth',2)
% hold on
% quiver3(0,0,0,0.01,0,0,'-g','LineWidth',2)
% 
% for i = 1:1:20
%     vec_cam = RT_cam(1:3,1:3,i)*vec_cam_ori';
%     vec_cam_x = RT_cam(1:3,1:3,i)*Rodrigues([0,pi/2,0])*vec_cam_ori';
%     vec_cam_y = -RT_cam(1:3,1:3,i)*Rodrigues([pi/2,0,0])*vec_cam_ori';
% %    scatter3(T_cam(i,1),T_cam(i,2),T_cam(i,3),10,'r')
%     hold on
%     quiver3(RT_cam(1,4,i),RT_cam(2,4,i),RT_cam(3,4,i),vec_cam(1),vec_cam(2),vec_cam(3),'-b','LineWidth',2)
%     hold on
%     quiver3(RT_cam(1,4,i),RT_cam(2,4,i),RT_cam(3,4,i),vec_cam_x(1),vec_cam_x(2),vec_cam_x(3),'-r','LineWidth',2)
%     hold on
%     quiver3(RT_cam(1,4,i),RT_cam(2,4,i),RT_cam(3,4,i),vec_cam_y(1),vec_cam_y(2),vec_cam_y(3),'-g','LineWidth',2)
% %      if i>1
% %          hold on
% %          plot3([RT_cam(1,4,i-1),RT_cam(1,4,i)],[RT_cam(2,4,i-1),RT_cam(2,4,i)],[RT_cam(3,4,i-1),RT_cam(3,4,i)],'g');
% %      end
% end
axis equal
imagePoints = zeros(49,2,17,2);
imagePoints_ori_1 = load('speckle_correct_1.txt');
imagePoints_ori_2 = load('speckle_correct_2.txt');
imagePoints_ori;
for i = 1:1:17
    imagePoints(:,:,i,1) = imagePoints_ori_1((i-1)*49+1:(i-1)*49+49,:);
    imagePoints(:,:,i,2) = imagePoints_ori_2((i-1)*49+1:(i-1)*49+49,:);
end
worldPoints = load('board.txt');
worldPoints;
[cameraParams, imagesUsed, estimationErrors] = estimateCameraParameters(imagePoints, worldPoints,'NumRadialDistortionCoefficients',3,'EstimateTangentialDistortion',true);

showExtrinsics(cameraParams)