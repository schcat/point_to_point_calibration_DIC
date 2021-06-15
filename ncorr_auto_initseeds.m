function [seedinfo,threaddiagram,outstate] = ncorr_auto_initseeds(reference,current,roi,radius,spacing,cutoff_diffnorm,cutoff_iteration,total_threads,enabled_stepanalysis,subsettrunc,num_img,total_imgs)
res_camera_w = 384;

    % Data ---------------------------------------------------------------%     
    % Initialize outputs
    seedinfo = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{}); % paramvector = [x y u v du/dx du/dy dv/dx dv/dy corrcoef]
    threaddiagram = [];
    outstate = out.cancelled; 
    
    % Initialize buffers
    % Initialize reduced ref
    ref_reduced = reference.reduce(spacing);
    % Initialize reduced ROI
    roi_reduced = roi.reduce(spacing);
    % Initialize seedinfo buffer
    seedinfo_prelim = struct('paramvector',{},'num_region',{},'num_thread',{},'computepoints',{});
    % Initialize thread diagram buffer
    threaddiagram_buffer = -ones(size(roi_reduced.mask));
      
    % Number of images successfully seeded
    num_imgs_success = inf;
    pos_region = [round(reference.width/2)-1 round(reference.height/2)-1];
    num_region = 0;
  
    % Get region mask
    regionmask = roi_reduced.get_regionmask(num_region);   
    % Initialize threaddiagram_prelim and threaddiagram preview - these
    % are modified in-place
    preview_threaddiagram = zeros(size(roi_reduced.mask));     
    % Update thread diagram and preview - must do this after setting the data
    %get_threaddiagram([]);

    
    % Get position - and make sure no two seeds overlap
    pos_seed = [];
    xx = 0;
    yy = 0;
    meanxx = sum(roi_reduced.mask,1);
    for i=1:1:size(roi_reduced.mask,2)
        xx = xx + i * meanxx(i);
    end
    meanyy = sum(roi_reduced.mask,2);
    for i=1:1:size(roi_reduced.mask,1)
        yy = yy + i * meanyy(i);
    end
    m=sum(sum(roi_reduced.mask));
%    meanx=round(xx/m)*(spacing+1);
%    meany=round(yy/m)*(spacing+1);
    
    meanx=round(xx/m)*(spacing+1);
    meany=round(yy/m)*(spacing+1);
    
    % Convert to zero based indexing and round
    pos_buffer = [round(reference.width/2) round(reference.height/2)];
    pos_seed = vertcat(pos_seed,pos_buffer); %#ok<AGROW>

        
    % Convert pos_seed to regular coordinates since they are WRT to
    % the reduced ROI at this point
%    pos_seed = pos_seed*(spacing+1);
    pos_seed = [meanx,meany];
    
                                                                   
    [seedinfo_analy,convergence_prelim,outstate_seeds] = ncorr_alg_seedanalysis(reference, ...
                                                                                     current, ...
                                                                                     roi, ...
                                                                                     num_region, ...
                                                                                     pos_seed, ...
                                                                                     radius, ...
                                                                                     cutoff_diffnorm, ...
                                                                                     cutoff_iteration, ...
                                                                                     enabled_stepanalysis, ...
                                                                                     subsettrunc, ...
                                                                                     num_img, ...
                                                                                     total_imgs);            
                                                                                       
    
    % Take minimum of num_imgs_success and
    % seedinfo_buffer. Buffer can be more or less than
    % num_imgs_success
    num_imgs_success = min(num_imgs_success,size(seedinfo_analy,3));

    % Clear out other images in prelim and buffer
%    if (~isempty(seedinfo_prelim))
%        seedinfo_prelim = seedinfo_prelim(:,:,1:num_imgs_success);
%    end
    seedinfo_analy = seedinfo_analy(:,:,1:num_imgs_success);

    % Append seedinfo_buffer - append along 2nd dimension
%    seedinfo_prelim = horzcat(seedinfo_prelim,seedinfo_analy); %#ok<AGROW>

    % Merge threaddiagram from previous iteration
 %   threaddiagram(threaddiagram_buffer ~= -1) = threaddiagram_buffer(threaddiagram_buffer ~= -1); 
 %   pad=0;
    if reference.width == 798
        pad = 15;
    elseif reference.width == 1024
        pad = 20;
    elseif reference.width == 1000
        pad = 100;
	pad_n = pad/(spacing+1);
        width = (reference.width-2*pad)/(spacing+1);
        height = (reference.height-2*pad)/(spacing+1);
        width_expand = (reference.width)/(spacing+1);
        threaddiagram = [-1*ones(pad_n,width_expand);-1*ones(height,pad_n),zeros(height,width),-1*ones(height,pad_n);-1*ones(pad_n,width_expand)];
    elseif reference.width == 750
        pad = 15;
        width = (reference.width)/(spacing+1)-2*pad;
        height = (reference.height)/(spacing+1)-2*pad;
        width_expand = width + 2*pad;
        threaddiagram = [-1*ones(pad,width_expand);-1*ones(height,pad),zeros(height,width),-1*ones(height,pad);-1*ones(pad,width_expand)];
    else
        threaddiagram = roi.mask(1:spacing+1:size(roi.mask,1),1:spacing+1:size(roi.mask,2))-1;
    end
    
    for i = 0:size(seedinfo_analy,3)-1
        for j = 0:size(seedinfo_analy,1)-1
         %   seedinfo_analy(j+1,1,i+1).computepoints = width*height;
            seedinfo_analy(j+1,1,i+1).computepoints = sum(sum(roi_reduced.mask));
        end
    end
        
    % Set outputs
    for i = 0:size(seedinfo_analy,1)-1
        for j = 0:size(seedinfo_analy,2)-1
            for k = 0:size(seedinfo_analy,3)-1
                seedinfo(i+1,j+1,k+1) = seedinfo_analy(i+1,j+1,k+1);
            end
        end
    end
    outstate = out.success;
 
end
