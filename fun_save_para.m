function fun_save_para(para, count)
fid = fopen(['para_temp_',num2str(count/125),'.txt'],'w');
para_length = size(para,2);
for k = 1:1:para_length
    fprintf(fid, '%f ',para(1,k));
end
fclose(fid);
