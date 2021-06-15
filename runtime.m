pathname = 'image/ori_0315/';
radius = 40; %DIC subset radius
spacing = 0;
filename_ref = ['simulate_speckle_patent_proj_1.png'];
filename_roi = ['ROI_proj_1.png'];
filename_cur = ['temp_1.jpg'];
tic
displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing);
toc
save(['results/ori_ref_mat/speckle_map_16.mat'],'displacements');