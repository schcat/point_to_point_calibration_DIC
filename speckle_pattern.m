%此程序用来生成散斑图
%r是散斑区域像素数，pad是边缘宽度，k_num是散斑数量，D是散斑半径
%相关论文Camera calibration using synthetic random speckle pattern and digital image correlation
r = 2048;
%pad = 256;
pad = 40;
value = zeros(r+pad*2);
x_pix = r;
y_pix = r;
k_num = 350000;
D = 2;

for k = 1:k_num
    xi = rand(1) * x_pix + pad;
    yi = rand(1) * y_pix + pad;
    Ii = rand(1) * 255;
    xi_ro = round(xi);
    yi_ro = round(yi);
    for x = xi_ro - 4*D:1:xi_ro + 4*D
        for y = yi_ro - 4*D:1:yi_ro + 4*D
            value(x,y) = value(x,y) + Ii * exp(-((x-xi)^2+(y-yi)^2)/(D^2));
            if value(x,y) > 255
                value(x,y) = 255;
            end
        end
    end
end

imshow(value/256);
%imwrite(mat2gray(value), 'image/speckle_pattern_2048_40_pad_12_350000.png');
imwrite(mat2gray(value), 'image/speckle_pattern_2048_40_pad_8_350000.tif');