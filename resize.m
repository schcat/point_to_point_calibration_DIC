value = imread('image/speckle_pattern_2048_40_pad_8_350000.png');
re = imresize(value, 0.5, "bilinear");
imwrite(re, 'image/speckle_pattern_1024_20_pad_8_350000.png');