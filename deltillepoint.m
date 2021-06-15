I = imread('distance/2_matlab_speckle_result/deltille_out/correct_97_matlab_distort.png');
imshow(I)
hold on
a=load('distance/2_matlab_speckle_result/97_all.orpc');
plot(a(:,1)+1,a(:,2)+1,'s','markersize',5);
set(gca,'YDir','reverse')
grid on