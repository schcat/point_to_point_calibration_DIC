function displacements = fun_dic(pathname, filename_ref, filename_roi ,filename_cur, radius, spacing, mat_name, k);

obj.reference = struct('imginfo',{},'roi',{});  
% Current image(s) info
obj.current = struct('imginfo',{},'roi',{});    
% DIC data            
obj.data_dic = struct('displacements',struct('plot_u_dic',{}, ...                                                % U plot after DIC (can be regular or eulerian)
                                             'plot_v_dic',{}, ...                                                % V plot after DIC (can be regular or eulerian)
                                             'plot_corrcoef_dic',{}, ...                                         % Correlaton coef plot after DIC (can be regular or eulerian)
                                             'roi_dic',{}, ...                                                   % ROI after DIC
                                             'plot_u_ref_formatted',{}, ...                                      % formatted U plot WRT reference image used for displaying displacements
                                             'plot_v_ref_formatted',{}, ...                                      % formatted V plot WRT reference image used for displaying displacements
                                             'roi_ref_formatted',{}, ...                                         % ROI after formatting used for plotting displacements
                                             'plot_u_cur_formatted',{}, ...                                      % formatted U plot WRT current image used for displaying displacements
                                             'plot_v_cur_formatted',{}, ...                                      % formatted V plot WRT current image used for displaying displacements
                                             'roi_cur_formatted',{}), ...                                        % ROI after formatting used for plotting displacements
                      'dispinfo',struct('type',{}, ...                                                           % Type of DIC: either regular or backward
                                        'radius',{}, ...                                                         % Radius for DIC
                                        'spacing',{}, ...                                                        % Spacing for DIC                                                    
                                        'cutoff_diffnorm',{}, ...                                                % Cutoff for norm of the difference vector
                                        'cutoff_iteration',{}, ...                                               % Cutoff for the number of iterations
                                        'total_threads',{}, ...                                                  % Number of threads for computation
                                        'stepanalysis',struct('enabled',{},'type',{},'auto',{},'step',{}), ...   % Indicates whether or not to implement step analysis for high strain
                                        'subsettrunc',{}, ...                                                    % Indicates whether or not to implement subset truncation for DIC analysis                
                                        'imgcorr',struct('idx_ref',{},'idx_cur',{}), ...                         % Image correspondences
                                        'pixtounits',{}, ...                                                     % Ratio of "units" to pixels. Assumes pixels are square
                                        'units',{}, ...                                                          % String to display units
                                        'cutoff_corrcoef',{}, ...                                                % Correlation coefficient cutoff for each formatted displacement plot
                                        'lenscoef',{}), ...                                                      % Radial lens distortion coefficient
                      'strains',struct('plot_exx_ref_formatted',{}, ...                                          % Exx Green-Lagragian strain plot
                                       'plot_exy_ref_formatted',{}, ...                                          % Exy Green-Lagragian strain plot
                                       'plot_eyy_ref_formatted',{}, ...                                          % Exy Green-Lagragian strain plot
                                       'roi_ref_formatted',{}, ...                                               % ROI used for plotting strains 
                                       'plot_exx_cur_formatted',{}, ...                                          % Exx Eulerian-Almansi strain plot 
                                       'plot_exy_cur_formatted',{}, ...                                          % Exy Eulerian-Almansi strain plot 
                                       'plot_eyy_cur_formatted',{}, ...                                          % Exy Eulerian-Almansi strain plot 
                                       'roi_cur_formatted',{}), ...                                              % ROI used for plotting strains                          
                      'straininfo',struct('radius',{}, ...                                                       % Strain radius used for calculating strains
                                          'subsettrunc',{}));

type_ref = 'file';
type_cur = 'file';
type_roi = 'file';

%obj.reference(1).imginfo = ncorr_class_img;
ref_prelim = ncorr_class_img;
ref_prelim.set_img(type_ref,struct('img',imread(fullfile(pathname,filename_ref)),'name',filename_ref,'path',pathname)); 
obj.reference(1).imginfo = ref_prelim;

mask_img_prelim = ncorr_class_img;
mask_img_prelim.set_img(type_roi,struct('img',imread(fullfile(pathname,filename_roi)),'name',filename_roi,'path',pathname)); 
mask_load_prelim = im2bw(mask_img_prelim.get_gs());
roi_prelim = ncorr_class_roi;
roi_prelim.set_roi('load',struct('mask',mask_load_prelim,'cutoff',2000));
obj.reference(1).roi = roi_prelim;

cur_prelim = ncorr_class_img;
cur_prelim.set_img(type_cur,struct('img',imread(fullfile(pathname,filename_cur)),'name',filename_cur,'path',pathname)); 
obj.current(1).imginfo = cur_prelim;

imgs = [obj.reference obj.current];

cutoff_diffnorm = 1e-5;
cutoff_iteration = 50;
total_threads = 1;
enabled_stepanalysis = false;
subsettrunc = false;
num_img = 1;
total_imgs = 1;
params_init = [];

[displacements,rois_dic,seedinfo,outstate] = ncorr_alg_dicanalysis(imgs,radius,spacing,cutoff_diffnorm,cutoff_iteration,total_threads,enabled_stepanalysis,subsettrunc,num_img,total_imgs,params_init);
if spacing > 0
    displacements.plot_u = imresize(displacements.plot_u,spacing+1,'bicubic');
    displacements.plot_v = imresize(displacements.plot_v,spacing+1,'bicubic');
    displacements.plot_corrcoef = imresize(displacements.plot_corrcoef,spacing+1,'bicubic');
end
save([mat_name,num2str(k),'.mat'],'displacements');
