m_1 = load('/home/wsco/cnn/calibration/speckle_camera_stereo_calibration_points/build/speckle_ori_1_1_s.txt');
m_2 = load('/home/wsco/cnn/calibration/speckle_camera_stereo_calibration_points/build/speckle_ori_2_1_s.txt');
%m_1 = load('/home/wsco/cnn/calibration/speckle_camera_stereo_calibration_points/build/speckle_correct_1_1_s.txt');
%m_2 = load('/home/wsco/cnn/calibration/speckle_camera_stereo_calibration_points/build/speckle_correct_2_1_s.txt');
n=14;
m_1_atom = zeros(n,49,2);
m_2_atom = zeros(n,49,2);
points_assem = zeros(84*n,8);
points_assem_atom = zeros(9,84,8);
for k = 1:1:n
    m_1_atom(k,:,:) = m_1((k-1)*49+1:k*49,:,:);
    m_2_atom(k,:,:) = m_2((k-1)*49+1:k*49,:,:);
end
for k = 1:1:n
    for i = 1:1:7
        for j = 1:1:6
            points_assem_atom(k,6*(i-1)+j,:) = [m_1_atom(k,7*(i-1)+j,1),m_1_atom(k,7*(i-1)+j,2),...
                                                m_1_atom(k,7*(i-1)+j+1,1),m_1_atom(k,7*(i-1)+j+1,2),...
                                                m_2_atom(k,7*(i-1)+j,1),m_2_atom(k,7*(i-1)+j,2),...
                                                m_2_atom(k,7*(i-1)+j+1,1),m_2_atom(k,7*(i-1)+j+1,2)];
        end
    end
    for i = 1:1:6
        for j = 1:1:7
            points_assem_atom(k,42+7*(i-1)+j,:) = [m_1_atom(k,7*(i-1)+j,1),m_1_atom(k,7*(i-1)+j,2),...
                                                m_1_atom(k,7*i+j,1),m_1_atom(k,7*i+j,2),...
                                                m_2_atom(k,7*(i-1)+j,1),m_2_atom(k,7*(i-1)+j,2),...
                                                m_2_atom(k,7*i+j,1),m_2_atom(k,7*i+j,2)];
        end
    end
end
for k = 1:1:n
    points_assem((k-1)*84+1:k*84,:) = points_assem_atom(k,:,:);
end
fid = fopen('/home/wsco/cnn/implementation/point_distance/build/speckle_test_ori_correct.txt','w');
for k = 1:1:84*n
    fprintf(fid, '%f %f %f %f %f %f %f %f\n',points_assem(k,1),points_assem(k,2),points_assem(k,3),points_assem(k,4),...
                                            points_assem(k,5),points_assem(k,6),points_assem(k,7),points_assem(k,8));
end
fclose(fid);
