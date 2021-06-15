fid = fopen('board_fullpixel.txt','w');
for i = 10:10:790
    for j = 10:10:790
        uu = i*0.06;
        vv = j*0.06;
        fprintf(fid, '%f %f\n',uu,vv);
    end
end
fclose(fid);