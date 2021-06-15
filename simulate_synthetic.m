value = imread('image/speckle_pattern_4000_pad_0111_20_15000_2000.png');
alpha = 0.6;
value = alpha * value + (1-alpha) * 255;
R = 12
for i = 1:1:7
    for j = 1:1:7
        for k_i = i*200-200:1:i*200+200
            for k_j = j*200-200:1:j*200+200
                dist = sqrt((k_i - i*200)^2 + (k_j - j*200)^2);
                if dist < R || dist == R
                    value(200+k_i,200+k_j) = 0;
                end
            end
        end
    end
end
imwrite(mat2gray(value), 'simulate_synthetic.png');