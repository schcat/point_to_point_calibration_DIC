%width = 1616;
%height = 1077;
width = 1024;
height = 1024;
pad = 20;
width_expand = width + 2*pad;
height_expand = height + 2*pad;
I = [zeros(pad,width_expand);zeros(height,pad),ones(height,width)*255,zeros(height,pad);zeros(pad,width_expand)];
imwrite(I,'image/ROI_1024_1024_20.png');
