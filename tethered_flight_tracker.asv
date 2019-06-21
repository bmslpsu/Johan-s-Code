function tethered_flight_tracker(mov_folder,mov_nr)

    close all;
    
    %----------------------------------------------------------------------
    % High-speed video tracker for tethered flight using 3D hull
    % reconstruction, minimum bounding box fitting and fast ICP fitting.
    % By Johan Melis, 12/21/2017
    %----------------------------------------------------------------------
    
    addpath('C:\Users\Flyami\Documents\tethered_flight_tracker\toolbox_graph');
    addpath('C:\Users\Flyami\Documents\tethered_flight_tracker\minboundbox');
    addpath('C:\Users\Flyami\Documents\tethered_flight_tracker\fast_icp_tolga_birdal');
    
    % Load background images:
    
    cd([mov_folder '/background']);
    
    bckg = cell(3,1);

    bckg{1} = 1-im2double(imread('background_cam_1.tif')); %inverts the colors of the image
    bckg{2} = 1-im2double(imread('background_cam_2.tif')); %this probably has to do with segmentation or binirization later on
    bckg{3} = 1-im2double(imread('background_cam_3.tif'));%
    %turns the background black
    % Ask user to draw tether and body masks:
    %for each image we have an image, body and tether mask (WS)
    masks = user_input_masks(mov_folder, mov_nr);
    
    % Ask user to construct body and wing models:
    %wing and wing hinges location in 3D, this is only valid for the rigid
    %tether (WS)
    %also its manual....f....
    model = user_input_model(mov_folder, mov_nr);
    
    % Load focal grid:
    
    try
        % Try to load the focal grid matrix:
        %but what the f1-1ck does a focal grid do?(WS)
        cd([mov_folder '/calibration']);
        disp('Loading focal grid, the hell does it do?')
        f_grid = load('Focal_grid.mat');
        
    catch
        
        % Construct focal grid
        find_focal_grid([mov_folder '/calibration'],model.wing_hinge_center);
        cd([mov_folder '/calibration']);
        f_grid = load('Focal_grid.mat');
        
    end
    
    % Create a solutions folder if it doesn't exist already:
    try
        cd([mov_folder '/mov_' int2str(mov_nr) '/solutions']);
    catch
        cd([mov_folder '/mov_' int2str(mov_nr)]);
        mkdir('solutions');
    end
    
    % Start saving batches until
    batch_size = 50;
    %flips the array using fliplr
    pre_trigger_batches = fliplr([(16375-batch_size+1):(-batch_size):8188 8188; 16375:(-batch_size):8188]);
    post_trigger_batches = [0:batch_size:8187; (batch_size-1):batch_size:8187 8187];
    
    % drop the first and last batch:
    pre_trigger_batches(:,1) = [];
    post_trigger_batches(:,end) = [];
    
    N_pre_trigger_batches = size(pre_trigger_batches,2);
    N_post_trigger_batches = size(post_trigger_batches,2);
    
    disp(['Number of batches: ' int2str(N_pre_trigger_batches+N_post_trigger_batches)]);
    
    for i = 1:(N_pre_trigger_batches+N_post_trigger_batches)
        if i <=N_pre_trigger_batches
            k = i;
            tic
            extract_wing_motion(mov_folder,mov_nr,bckg,model,f_grid,masks,pre_trigger_batches(1,k),pre_trigger_batches(2,k),i);
            toc
        elseif i > N_pre_trigger_batches
            k = i-N_pre_trigger_batches;
            tic
            extract_wing_motion(mov_folder,mov_nr,bckg,model,f_grid,masks,post_trigger_batches(1,k),post_trigger_batches(2,k),i);
            toc
        end
    end

    disp('finished');

end

function extract_wing_motion(mov_folder,mov_nr,bckg,model,f_grid,masks,start_frame,end_frame,batch_nr)

    disp(['Analyzing batch ' int2str(batch_nr)]);
    
    % Load images:
    disp('load images');
    frames = load_frames(mov_folder,mov_nr,masks,bckg,start_frame,end_frame);
    
    % Extract wing contours from the images:
    %these are used to create the voxels later on(WS)
    disp('extract image contours');
    w_cont = extract_wing_contours(frames,masks); 
    
    % Project the images to voxels and select left and right wing voxels:
    disp('select wing voxels');
    wing_voxels = select_wing_voxels(frames,f_grid,w_cont,masks);
    
    % Obtain the wing orientation using pca with dropout:
    disp('obtain wing orientation');
    sol = wing_fitting(wing_voxels,model);
    
    % Save results:
    disp('save batch results');
    cd([mov_folder '/mov_' int2str(mov_nr) '/solutions']);
    file_name = ['sol_batch_' int2str(batch_nr) '.mat'];
    save(file_name,'sol');
    
%     figure()
%     hold on
%     subplot(4,1,1); hold on
%     plot(start_frame:end_frame,sol.q_wing_L1(1,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_L2(1,:),'Color','b')
%     hold off
%     subplot(4,1,2); hold on
%     plot(start_frame:end_frame,sol.q_wing_L1(2,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_L2(2,:),'Color','b')
%     hold off
%     subplot(4,1,3); hold on
%     plot(start_frame:end_frame,sol.q_wing_L1(3,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_L2(3,:),'Color','b')
%     hold off
%     subplot(4,1,4); hold on
%     plot(start_frame:end_frame,sol.q_wing_L1(4,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_L2(4,:),'Color','b')
%     hold off
%     hold off
%     
%     figure()
%     hold on
%     subplot(4,1,1); hold on
%     plot(start_frame:end_frame,sol.q_wing_R1(1,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_R2(1,:),'Color','b')
%     hold off
%     subplot(4,1,2); hold on
%     plot(start_frame:end_frame,sol.q_wing_R1(2,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_R2(2,:),'Color','b')
%     hold off
%     subplot(4,1,3); hold on
%     plot(start_frame:end_frame,sol.q_wing_R1(3,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_R2(3,:),'Color','b')
%     hold off
%     subplot(4,1,4); hold on
%     plot(start_frame:end_frame,sol.q_wing_R1(4,:),'Color','r')
%     plot(start_frame:end_frame,sol.q_wing_R2(4,:),'Color','b')
%     hold off
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); hold on
%     plot(start_frame:end_frame,radtodeg(sol.theta_L),'Color','b')
%     hold off
%     subplot(3,1,2); hold on
%     plot(start_frame:end_frame,radtodeg(sol.eta_L1),'Color','r')
%     plot(start_frame:end_frame,radtodeg(sol.eta_L2),'Color','b')
%     hold off
%     subplot(3,1,3); hold on
%     plot(start_frame:end_frame,radtodeg(sol.phi_L),'Color','b')
%     hold off
%     hold off
%     
%     figure()
%     hold on
%     subplot(3,1,1); hold on
%     plot(start_frame:end_frame,radtodeg(sol.theta_R),'Color','b')
%     hold off
%     subplot(3,1,2); hold on
%     plot(start_frame:end_frame,radtodeg(sol.eta_R1),'Color','r')
%     plot(start_frame:end_frame,radtodeg(sol.eta_R2),'Color','b')
%     hold off
%     subplot(3,1,3); hold on
%     plot(start_frame:end_frame,radtodeg(sol.phi_R),'Color','b')
%     hold off
%     hold off
    
end

function wing_sol = wing_fitting(wing_voxels,model)
%takes voxels and fits them. won't be needed for my modification(WS)
    N_frames = length(wing_voxels.left);
    
    R_strk = model.R_strk;
    q_strk = model.q_strk;
    R_body = model.R_body;
    q_body = model.q_body;
    
    joint_L = model.wing_hinge_L;
    joint_R = model.wing_hinge_R;
    
    wing_sol = {};
    
    wing_sol.R_wing_L1 = zeros(3,3,N_frames);
    wing_sol.R_wing_L2 = zeros(3,3,N_frames);
    wing_sol.R_wing_R1 = zeros(3,3,N_frames);
    wing_sol.R_wing_R2 = zeros(3,3,N_frames);
    
    wing_sol.q_wing_L1 = zeros(4,N_frames);
    wing_sol.q_wing_L2 = zeros(4,N_frames);
    wing_sol.q_wing_R1 = zeros(4,N_frames);
    wing_sol.q_wing_R2 = zeros(4,N_frames);
    
    wing_sol.R_wing_strk_L1 = zeros(3,3,N_frames);
    wing_sol.R_wing_strk_L2 = zeros(3,3,N_frames);
    wing_sol.R_wing_strk_R1 = zeros(3,3,N_frames);
    wing_sol.R_wing_strk_R2 = zeros(3,3,N_frames);
    
    wing_sol.q_wing_strk_L1 = zeros(4,N_frames);
    wing_sol.q_wing_strk_L2 = zeros(4,N_frames);
    wing_sol.q_wing_strk_R1 = zeros(4,N_frames);
    wing_sol.q_wing_strk_R2 = zeros(4,N_frames);
    
    wing_sol.theta_L = zeros(1,N_frames);
    wing_sol.eta_L1  = zeros(1,N_frames);
    wing_sol.eta_L2  = zeros(1,N_frames);
    wing_sol.phi_L   = zeros(1,N_frames);
    
    wing_sol.theta_R = zeros(1,N_frames);
    wing_sol.eta_R1  = zeros(1,N_frames);
    wing_sol.eta_R2  = zeros(1,N_frames);
    wing_sol.phi_R   = zeros(1,N_frames);
    
    for i = 1:N_frames

        % Find rotation matrices and quaternions using minimum bounding box
        if isempty(wing_voxels.left{i}) == 0
            try 
            wing_vox_L = remove_outlying_voxels(full(wing_voxels.left{i}(1:3,:)),joint_L);
            [~,corners_L,~,~,~] = minboundbox(wing_vox_L(1,:),wing_vox_L(2,:),wing_vox_L(3,:),'v',1);
            wing_ref_L = find_wing_reference_frame(wing_vox_L,corners_L',joint_L,q_body,q_strk,0);
            catch
                wing_ref_L.R_wing_1 = nan(3,3);
                wing_ref_L.R_wing_2 = nan(3,3);
                wing_ref_L.q_wing_1 = nan(4,1);
                wing_ref_L.q_wing_2 = nan(4,1);
                wing_ref_L.R_wing_strk_1 = nan(3,3);
                wing_ref_L.R_wing_strk_2 = nan(3,3);
                wing_ref_L.q_wing_strk_1 = nan(4,1);
                wing_ref_L.q_wing_strk_2 = nan(4,1);
                wing_ref_L.theta = nan;
                wing_ref_L.eta_1 = nan;
                wing_ref_L.eta_2 = nan;
                wing_ref_L.phi = nan;
            end
        else
            wing_ref_L.R_wing_1 = nan(3,3);
            wing_ref_L.R_wing_2 = nan(3,3);
            wing_ref_L.q_wing_1 = nan(4,1);
            wing_ref_L.q_wing_2 = nan(4,1);
            wing_ref_L.R_wing_strk_1 = nan(3,3);
            wing_ref_L.R_wing_strk_2 = nan(3,3);
            wing_ref_L.q_wing_strk_1 = nan(4,1);
            wing_ref_L.q_wing_strk_2 = nan(4,1);
            wing_ref_L.theta = nan;
            wing_ref_L.eta_1 = nan;
            wing_ref_L.eta_2 = nan;
            wing_ref_L.phi = nan;
        end
        
        if isempty(wing_voxels.right{i}) == 0
            try
                wing_vox_R = remove_outlying_voxels(full(wing_voxels.right{i}(1:3,:)),joint_R);
                [~,corners_R,~,~,~] = minboundbox(wing_vox_R(1,:),wing_vox_R(2,:),wing_vox_R(3,:),'v',1);
                wing_ref_R = find_wing_reference_frame(wing_vox_R,corners_R',joint_R,q_body,q_strk,1);
            catch
                wing_ref_R.R_wing_1 = nan(3,3);
                wing_ref_R.R_wing_2 = nan(3,3);
                wing_ref_R.q_wing_1 = nan(4,1);
                wing_ref_R.q_wing_2 = nan(4,1);
                wing_ref_R.R_wing_strk_1 = nan(3,3);
                wing_ref_R.R_wing_strk_2 = nan(3,3);
                wing_ref_R.q_wing_strk_1 = nan(4,1);
                wing_ref_R.q_wing_strk_2 = nan(4,1);
                wing_ref_R.theta = nan;
                wing_ref_R.eta_1 = nan;
                wing_ref_R.eta_2 = nan;
                wing_ref_R.phi = nan;
            end
        else
            wing_ref_R.R_wing_1 = nan(3,3);
            wing_ref_R.R_wing_2 = nan(3,3);
            wing_ref_R.q_wing_1 = nan(4,1);
            wing_ref_R.q_wing_2 = nan(4,1);
            wing_ref_R.R_wing_strk_1 = nan(3,3);
            wing_ref_R.R_wing_strk_2 = nan(3,3);
            wing_ref_R.q_wing_strk_1 = nan(4,1);
            wing_ref_R.q_wing_strk_2 = nan(4,1);
            wing_ref_R.theta = nan;
            wing_ref_R.eta_1 = nan;
            wing_ref_R.eta_2 = nan;
            wing_ref_R.phi = nan;
        end

        % Save solutions:
        
        if isempty(wing_ref_L) == 0

            wing_sol.R_wing_L1(:,:,i) = wing_ref_L.R_wing_1;
            wing_sol.R_wing_L2(:,:,i) = wing_ref_L.R_wing_2;

            wing_sol.q_wing_L1(:,i) = wing_ref_L.q_wing_1;
            wing_sol.q_wing_L2(:,i) = wing_ref_L.q_wing_2;

            wing_sol.R_wing_strk_L1(:,:,i) = wing_ref_L.R_wing_strk_1;
            wing_sol.R_wing_strk_L2(:,:,i) = wing_ref_L.R_wing_strk_2;

            wing_sol.q_wing_strk_L1(:,i) = wing_ref_L.q_wing_strk_1;
            wing_sol.q_wing_strk_L2(:,i) = wing_ref_L.q_wing_strk_2;

            wing_sol.theta_L(:,i) = wing_ref_L.theta;
            wing_sol.eta_L1(:,i)  = wing_ref_L.eta_1;
            wing_sol.eta_L2(:,i)  = wing_ref_L.eta_2;
            wing_sol.phi_L(:,i)   = wing_ref_L.phi;
        
        else
            
            wing_sol.R_wing_L1(:,:,i) = nan(3,3,1);
            wing_sol.R_wing_L2(:,:,i) = nan(3,3,1);

            wing_sol.q_wing_L1(:,i) = nan(4,1);
            wing_sol.q_wing_L2(:,i) = nan(4,1);

            wing_sol.R_wing_strk_L1(:,:,i) = nan(3,3,1);
            wing_sol.R_wing_strk_L2(:,:,i) = nan(3,3,1);

            wing_sol.q_wing_strk_L1(:,i) = nan(4,1);
            wing_sol.q_wing_strk_L2(:,i) = nan(4,1);

            wing_sol.theta_L(:,i) = nan;
            wing_sol.eta_L1(:,i)  = nan;
            wing_sol.eta_L2(:,i)  = nan;
            wing_sol.phi_L(:,i)   = nan;
            
        end
        
        if isempty(wing_ref_R) == 0
        
            wing_sol.R_wing_R1(:,:,i) = wing_ref_R.R_wing_1;
            wing_sol.R_wing_R2(:,:,i) = wing_ref_R.R_wing_2;

            wing_sol.q_wing_R1(:,i) = wing_ref_R.q_wing_1;
            wing_sol.q_wing_R2(:,i) = wing_ref_R.q_wing_2;

            wing_sol.R_wing_strk_R1(:,:,i) = wing_ref_R.R_wing_strk_1;
            wing_sol.R_wing_strk_R2(:,:,i) = wing_ref_R.R_wing_strk_2;

            wing_sol.q_wing_strk_R1(:,i) = wing_ref_R.q_wing_strk_1;
            wing_sol.q_wing_strk_R2(:,i) = wing_ref_R.q_wing_strk_2;

            wing_sol.theta_R(i) = wing_ref_R.theta;
            wing_sol.eta_R1(i)  = wing_ref_R.eta_1;
            wing_sol.eta_R2(i)  = wing_ref_R.eta_2;
            wing_sol.phi_R(i)   = wing_ref_R.phi;
        
        else
            
            wing_sol.R_wing_R1(:,:,i) = nan(3,3);
            wing_sol.R_wing_R2(:,:,i) = nan(3,3);

            wing_sol.q_wing_R1(:,i) = nan(4,1);
            wing_sol.q_wing_R2(:,i) = nan(4,1);

            wing_sol.R_wing_strk_R1(:,:,i) = nan(3,3);
            wing_sol.R_wing_strk_R2(:,:,i) = nan(3,3);

            wing_sol.q_wing_strk_R1(:,i) = nan(4,1);
            wing_sol.q_wing_strk_R2(:,i) = nan(4,1);

            wing_sol.theta_R(i) = nan;
            wing_sol.eta_R1(i)  = nan;
            wing_sol.eta_R2(i)  = nan;
            wing_sol.phi_R(i)   = nan;
            
        end
        
    end

end

function wing_ref = find_wing_reference_frame(vox,corners,joint,q_body,q_strk,left_or_right)

    % Calculate center bounding box:
    center = mean(corners,2);
        
    % Calculate distance from joint:
    corner_vec = bsxfun(@minus,corners,joint);
    dist_corner = sqrt(sum(corner_vec.^2,1));
        
    % Sort based on joint distance:
    [~,dist_corner_ind] = sort(dist_corner);
        
    %corners_base = corners(:,dist_corner_ind(1:4));
    corners_tip = corners(:,dist_corner_ind(5:8));
        
    center_tip = mean(corners_tip,2);
    
    % Select corner pairs based on the distance between the pairs
    
    tip_corner_dist = [ sqrt(sum(bsxfun(@minus,corners_tip(:,2:4),corners_tip(:,1)).^2,1)) ...
                        sqrt(sum(bsxfun(@minus,corners_tip(:,3:4),corners_tip(:,2)).^2,1)) ...
                        sqrt(sum(bsxfun(@minus,corners_tip(:,4),corners_tip(:,3)).^2,1)) ];
                        
    tip_corner_pairs = [ 2 3 4 3 4 4; 1 1 1 2 2 3 ];
                    
    [~,tip_corner_ind] = sort(tip_corner_dist,'descend');
    
    corner_diag_points_1 = [corners_tip(:,tip_corner_pairs(1,tip_corner_ind(1))) corners_tip(:,tip_corner_pairs(2,tip_corner_ind(1)))];
    corner_diag_points_2 = [corners_tip(:,tip_corner_pairs(1,tip_corner_ind(2))) corners_tip(:,tip_corner_pairs(2,tip_corner_ind(2)))];
    
    % Y-axis is formed by the vector between the joint and the mean
    % position of corners_tip:
        
    if left_or_right == 0
        Y_wing = center_tip-joint;
        Y_wing = Y_wing/norm(Y_wing);
        
        corner_vec_points_1 = corner_diag_points_1(:,1)-corner_diag_points_1(:,2);
        corner_vec_points_2 = corner_diag_points_2(:,1)-corner_diag_points_2(:,2);
        
        X_wing_1 = corner_vec_points_1/norm(corner_vec_points_1);
        X_wing_2 = corner_vec_points_2/norm(corner_vec_points_2);
        
        Z_wing_1 = cross(X_wing_1,Y_wing)/(norm(X_wing_1)*norm(Y_wing));
        Z_wing_2 = cross(X_wing_2,Y_wing)/(norm(X_wing_2)*norm(Y_wing));
        
    elseif left_or_right == 1
        Y_wing = -(center_tip-joint);
        Y_wing = Y_wing/norm(Y_wing);            
        
        corner_vec_points_1 = corner_diag_points_1(:,1)-corner_diag_points_1(:,2);
        corner_vec_points_2 = corner_diag_points_2(:,1)-corner_diag_points_2(:,2);
        
        X_wing_1 = corner_vec_points_1/norm(corner_vec_points_1);
        X_wing_2 = corner_vec_points_2/norm(corner_vec_points_2);
        
        Z_wing_1 = cross(X_wing_1,Y_wing)/(norm(X_wing_1)*norm(Y_wing));
        Z_wing_2 = cross(X_wing_2,Y_wing)/(norm(X_wing_2)*norm(Y_wing));
        
    end
    
    % Calculate the number of voxels within 0.05 mm of the diagonal plane
    % to determine the most likely wing orientation
    
%     vox_diag_1 = drop_points(vox,[X_wing_1 Y_wing Z_wing_1],center,0.1,3);
%     vox_diag_2 = drop_points(vox,[X_wing_2 Y_wing Z_wing_2],center,0.1,3);
%     
%     N_vox = [size(vox_diag_1,2) size(vox_diag_2,2)];
%     
%     [~, vox_sort_ind] = sort(N_vox,'descend');
    
    %if vox_sort_ind(1) == 1
        corner_vec_strk_1 = quat2mat(quat_multiply(q_strk,q_body))'*X_wing_1;
        if corner_vec_strk_1(3) < 0
            corner_vec_strk_1 = -corner_vec_strk_1;
        end
        corner_vec_strk_2 = quat2mat(quat_multiply(q_strk,q_body))'*X_wing_2;
        if corner_vec_strk_2(3) < 0
            corner_vec_strk_2 = -corner_vec_strk_2;
        end
%     elseif vox_sort_ind(1) == 2
%         corner_vec_strk_1 = quat2mat(quat_multiply(q_strk,q_body))'*X_wing_2;
%         if corner_vec_strk_1(3) < 0
%             corner_vec_strk_1 = -corner_vec_strk_1;
%         end
%         corner_vec_strk_2 = quat2mat(quat_multiply(q_strk,q_body))'*X_wing_1;
%         if corner_vec_strk_2(3) < 0
%             corner_vec_strk_2 = -corner_vec_strk_2;
%         end
%     end
    
    
    % Transform the wing's Y-axis to the strokeplane reference frame:
    Y_wing_strk = quat2mat(quat_multiply(q_strk,q_body))'*Y_wing;
    
    % Calculate the wing kinematic angles:
    
    if left_or_right == 0
    
        phi = atan2(-Y_wing_strk(1),Y_wing_strk(2));
        theta = asin(Y_wing_strk(3));
        Z_axis = [sin(phi)*sin(theta); cos(phi)*sin(theta); cos(theta)];
        eta_1_t = asin(norm(cross(corner_vec_strk_1,Z_axis))/(norm(corner_vec_strk_1)*norm(Z_axis)));
        eta_2_t = asin(norm(cross(corner_vec_strk_2,Z_axis))/(norm(corner_vec_strk_2)*norm(Z_axis)));
        
        eta_1 = min([eta_1_t eta_2_t]);
        eta_2 = max([eta_1_t eta_2_t]);
    
    elseif left_or_right == 1
        
        phi = atan2(Y_wing_strk(1),Y_wing_strk(2));
        theta = asin(-Y_wing_strk(3));
        Z_axis = [sin(phi)*sin(theta); cos(phi)*sin(theta); cos(theta)];
        eta_1_t = asin(norm(cross(corner_vec_strk_1,Z_axis))/(norm(corner_vec_strk_1)*norm(Z_axis)));
        eta_2_t = asin(norm(cross(corner_vec_strk_2,Z_axis))/(norm(corner_vec_strk_2)*norm(Z_axis)));
        
        eta_1 = min([eta_1_t eta_2_t]);
        eta_2 = max([eta_1_t eta_2_t]);
        
    end
    
    % Construct rotation matrices and quaternions:
    
    if left_or_right == 0

        R_wing_strk_1 = [ cos(eta_1)*cos(-theta) sin(eta_1)*sin(phi)-cos(eta_1)*cos(phi)*sin(-theta) cos(phi)*sin(eta_1)+cos(eta_1)*sin(-theta)*sin(phi); ...
                          sin(-theta) cos(-theta)*cos(phi) -cos(-theta)*sin(phi); ...
                          -cos(theta)*sin(eta_1) cos(eta_1)*sin(phi)+cos(phi)*sin(eta_1)*sin(-theta) cos(eta_1)*cos(phi)-sin(eta_1)*sin(-theta)*sin(phi)];
                      
        R_wing_strk_2 = [ cos(eta_2)*cos(-theta) sin(eta_2)*sin(phi)-cos(eta_2)*cos(phi)*sin(-theta) cos(phi)*sin(eta_2)+cos(eta_2)*sin(-theta)*sin(phi); ...
                          sin(-theta) cos(-theta)*cos(phi) -cos(-theta)*sin(phi); ...
                          -cos(-theta)*sin(eta_2) cos(eta_2)*sin(phi)+cos(phi)*sin(eta_2)*sin(-theta) cos(eta_2)*cos(phi)-sin(eta_2)*sin(-theta)*sin(phi)];

        % Rotation matrix to rotate the wing in vertical pitch position
        R_wing_vertical = [ 0 0 -1; 0 1 0; 1 0 0];
                  
    elseif left_or_right == 1

        R_wing_strk_1 = [ cos(eta_1)*cos(-theta) sin(eta_1)*sin(phi)-cos(eta_1)*cos(phi)*sin(-theta) cos(phi)*sin(eta_1)+cos(eta_1)*sin(-theta)*sin(phi); ...
                          sin(-theta) cos(-theta)*cos(phi) -cos(-theta)*sin(phi); ...
                          -cos(theta)*sin(eta_1) cos(eta_1)*sin(phi)+cos(phi)*sin(eta_1)*sin(-theta) cos(eta_1)*cos(phi)-sin(eta_1)*sin(-theta)*sin(phi)];
                      
        R_wing_strk_2 = [ cos(eta_2)*cos(-theta) sin(eta_2)*sin(phi)-cos(eta_2)*cos(phi)*sin(-theta) cos(phi)*sin(eta_2)+cos(eta_2)*sin(-theta)*sin(phi); ...
                          sin(-theta) cos(-theta)*cos(phi) -cos(-theta)*sin(phi); ...
                          -cos(-theta)*sin(eta_2) cos(eta_2)*sin(phi)+cos(phi)*sin(eta_2)*sin(-theta) cos(eta_2)*cos(phi)-sin(eta_2)*sin(-theta)*sin(phi)];
                      
        % Rotation matrix to rotate the wing in vertical pitch position
        R_wing_vertical = [ 0 0 1; 0 1 0; -1 0 0];
        
    end
    
    q_wing_strk_1 = quat_multiply(mat2quat(R_wing_strk_1),mat2quat(R_wing_vertical));
    q_wing_strk_2 = quat_multiply(mat2quat(R_wing_strk_2),mat2quat(R_wing_vertical));
    
    % Caclculate world to wing transformation matrix:
    q_wing_1 = quat_multiply(q_wing_strk_1,quat_multiply(q_strk,q_body));
    q_wing_2 = quat_multiply(q_wing_strk_2,quat_multiply(q_strk,q_body));
    R_wing_1 = quat2mat(q_wing_1);
    R_wing_2 = quat2mat(q_wing_2);
    
    % Return wing_ref structure:
    
    wing_ref = {};
    
    wing_ref.theta = theta;
    wing_ref.eta_1 = eta_1;
    wing_ref.eta_2 = eta_2;
    wing_ref.phi   = phi;
    
    wing_ref.R_wing_strk_1 = R_wing_strk_1;
    wing_ref.R_wing_strk_2 = R_wing_strk_2;
    wing_ref.q_wing_strk_1 = q_wing_strk_1; 
    wing_ref.q_wing_strk_2 = q_wing_strk_2;
    
    wing_ref.R_wing_1 = R_wing_1;
    wing_ref.R_wing_2 = R_wing_2;
    wing_ref.q_wing_1 = q_wing_1;
    wing_ref.q_wing_2 = q_wing_2;
                  
end

function wing_voxels = select_wing_voxels(frames,f_grid,w_cont,masks)
%where
%w_cont is the contour of each wing 
    N_frames = size(frames{1},3);
    
    % Reshape all frames to vectors:
    
    im1_vec = zeros(256^2,N_frames);
    im2_vec = zeros(256^2,N_frames);
    im3_vec = zeros(256^2,N_frames);
    im1_wing_vec = zeros(256^2,N_frames);
    im2_wing_vec = zeros(256^2,N_frames);
    im3_wing_vec_L = zeros(256^2,N_frames);
    im3_wing_vec_R = zeros(256^2,N_frames);
    
    for i = 1:N_frames
        %image is reshaped to vector and written onto the vector defined
        %above(WS)
        %frame is a cell of dim 3. each dimension has 50 frames of
        %images(WS)
        im1_vec(:,i) = reshape(frames{1}(:,:,i),256^2,1);
        im2_vec(:,i) = reshape(frames{2}(:,:,i),256^2,1);
        im3_vec(:,i) = reshape(frames{3}(:,:,i),256^2,1);
        
        im1_wing_vec(:,i) = reshape(w_cont{1}(:,:,i),256^2,1);
        im2_wing_vec(:,i) = reshape(w_cont{2}(:,:,i),256^2,1);
        im3_wing_vec_L(:,i) = reshape(masks{3}.im_left_side.*w_cont{3}(:,:,i),256^2,1);
        im3_wing_vec_R(:,i) = reshape((1-masks{3}.im_left_side).*w_cont{3}(:,:,i),256^2,1);
        
    end
    
    % Multiply the vectorized images with the voxel projection matrices:
    im1_proj = f_grid.vox_1_mat*sparse(im1_vec); %dim of this is 257^3, no idea why(WS)
    im2_proj = f_grid.vox_2_mat*sparse(im2_vec);
    im3_proj = f_grid.vox_3_mat*sparse(im3_vec);
    
    % WS: this is a bit of code I added to do some debugging stuff
    %The lack of alcohol in the lab bothers me
    %
    figure
    im_11=reshape(full(im1_proj(:,1)),[257,257,257]);
    idx = find(im_11);
    [X, Y, Z] = ind2sub(size(im_11), idx);
    pointsize = 20;
    scatter3(X(:), Y(:), Z(:), pointsize, im_11(idx));
    colormap(gray(256));
    
    figure
    im_22=reshape(full(im2_proj(:,1)),[257,257,257]);
    idx = find(im_22);
    [X, Y, Z] = ind2sub(size(im_22), idx);
    pointsize = 20;
    scatter3(X(:), Y(:), Z(:), pointsize, im_22(idx));
    colormap(gray(256));
    
    figure
    im_33=reshape(full(im3_proj(:,1)),[257,257,257]);
    idx = find(im_33);
    [X, Y, Z] = ind2sub(size(im_33), idx);
    pointsize = 20;
    scatter3(X(:), Y(:), Z(:), pointsize, im_33(idx));
    colormap(gray(256));
    
    im_3d=im1_proj(:,1).*im2_proj(:,1).*im3_proj(:,1);
    im_3d(im_3d<0.15)=0;
    im_33d=reshape(full(im_3d),[257,257,257]);
    idx = find(im_33d);
    [X, Y, Z] = ind2sub(size(im_33d), idx);
    pointsize = 20;
    figure
    scatter3(X(:), Y(:), Z(:), pointsize, im_33d(idx));
    colormap(gray(256));
    %Wael Modification ends here. still no alcohol
    im1_wing_proj = f_grid.vox_1_mat*sparse(im1_wing_vec);
    im2_wing_proj = f_grid.vox_2_mat*sparse(im2_wing_vec);
    im3_wing_proj_L = f_grid.vox_3_mat*sparse(im3_wing_vec_L);
    im3_wing_proj_R = f_grid.vox_3_mat*sparse(im3_wing_vec_R);
    
    % Select wing voxels left and right wing:
    
    wing_voxels = {};
    
    wing_voxels.left = cell(N_frames,1);
    wing_voxels.right = cell(N_frames,1);
    
    for i = 1:N_frames
        
        vox_select_L = ((im1_proj(:,i).*im2_proj(:,i).*im3_proj(:,i))>0)&((im1_wing_proj(:,i).*im3_wing_proj_L(:,i))>0 | (im2_wing_proj(:,i).*im3_wing_proj_L(:,i))>0);
        vox_select_R = ((im1_proj(:,i).*im2_proj(:,i).*im3_proj(:,i))>0)&((im1_wing_proj(:,i).*im3_wing_proj_R(:,i))>0 | (im2_wing_proj(:,i).*im3_wing_proj_R(:,i))>0);
        
        wing_voxels.left{i} = [f_grid.P_vox(:,vox_select_L); im1_proj(vox_select_L,i)'; im2_proj(vox_select_L,i)'; im3_proj(vox_select_L,i)'];
        wing_voxels.right{i} = [f_grid.P_vox(:,vox_select_R); im1_proj(vox_select_R,i)'; im2_proj(vox_select_R,i)'; im3_proj(vox_select_R,i)'];
        
    end

end

function vox_select = select_diagonal_voxels(vox,R)

    % Transform the voxels to the wing reference frame and only select the
    % voxels in or close to the XY plane in the wing reference frame:
    
    vox_wing_ref = R'*vox;
    
    vox_wing_std = std(vox_wing_ref,1,2);
    
    vox_select = vox(:,vox_wing_std(3,:)<=0.75);

end

function vox_select = remove_outlying_voxels(vox,joint)

    % Remove all voxels in a radius of 1.0 mm around the wing hinge
    
    vox_joint_dist = sqrt(sum(bsxfun(@minus,vox,joint).^2,1));
    
    vox_select = vox(:,vox_joint_dist>1.0);
    
    % Remove all voxels with a z score distance higher than 3
    
    z_score_vox = zscore(vox_select,1,2);
    
    z_dist_vox = sqrt(sum(z_score_vox.^2,1));
    
    vox_select = vox_select(:,z_dist_vox<=3);

end

function w_cont = extract_wing_contours(frames,masks)

    N_frames = size(frames{1},3);
    
    w_cont = cell(3,1);
    
    for i = 1:3
        
        w_cont{i} = zeros(256,256,N_frames);
        
        for j = 1:N_frames
            
            im_t = frames{i}(:,:,j); %im_t is the image j from camera i
            im_wing_blur = imgaussfilt(((1-masks{i}.body_mask).*im_t),1);
            %this function only works with binary images
            [B_wing,~] = bwboundaries(((im_wing_blur>=0.2).*im_wing_blur),'noholes');
            
            wing_cont_im = zeros(256,256);

            for k = 1:length(B_wing)
                if size(B_wing{k},1) > 60 % not really sure how he chose 60, maybe to filter out noise and other contours that arent part of the fly
                    bw = B_wing{k};
                    for m = 1:size(bw,1)
                        wing_cont_im(bw(m,1),bw(m,2)) = 1; %plots contour on the image
                    end
                end
            end
            %fills holes in an image which are the contours
            wing_mask_im = imfill(wing_cont_im,'holes');
            w_cont{i}(:,:,j) = bwmorph(wing_mask_im,'thicken',1); %applies a morphological operation to a BW image, thicken connects unconnected bodies together
            % check the documentation for a list of those
        end
    end

end

function frames = load_frames(mov_folder,mov_nr,masks,bckg,start_frame,end_frame)

    % Batch load frames from the 3 cameras:

    N_frames = end_frame-start_frame+1;

    frames = cell(3,1);
    
    for i = 1:3 %3 is the number of cameras 

        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_' int2str(i)]);

        frames{i} = zeros(256,256,N_frames); % size of the image recorded

        for j = start_frame:end_frame

            img_t = (1-im2double(imread(['frame_' int2str(j) '.bmp'])))-bckg{1}; %here bckg is 1-originalImage
            img_t(img_t<0.2) = 0;%this seems to be soem sort of binirization
            frames{i}(:,:,j-start_frame+1) = (masks{i}.img_mask.*(1-masks{i}.tether_mask)).*img_t;

        end
    
    end

end

function [point_1,point_2,point_3] = get_image_centers(mov_folder, mov_nr)

    cd([mov_folder '/mov_' int2str(mov_nr) '/cam_1']);
    img_1_t = imread('frame_0.bmp');
    cd([mov_folder '/mov_' int2str(mov_nr) '/cam_2']);
    img_2_t = imread('frame_0.bmp');
    cd([mov_folder '/mov_' int2str(mov_nr) '/cam_3']);
    img_3_t = imread('frame_0.bmp');
    
    disp('Locate body center of mass cam 1, double click to proceed');
    figure, imshow(img_1_t);
    h1 = impoint(gca,[]);
    point_1 = wait(h1);
    close all;
    
    disp('Locate body center of mass cam 2, double click to proceed');
    figure, imshow(img_2_t);
    h2 = impoint(gca,[]);
    point_2 = wait(h2);
    close all;
    
    disp('Locate body center of mass cam 3, double click to proceed');
    figure, imshow(img_3_t);
    h3 = impoint(gca,[]);
    point_3 = wait(h3);
    close all;

end

function model = user_input_model(mov_folder, mov_nr)

    % Get body center, wing hinge centers and wing lengths:
    %this is manual digitization of the fly model which includes wing
    %hinge, abdomen hinge, wing tip.....
    
    try 
        cd([mov_folder '/fly_model']);
        load('fly_model.mat');
    catch
        % Ask user to iterate trough frames until a frame comes up with a
        % good view of the wings fully extended:
        show_frames = 1;
        
        frame_nr = 0;
        %reads teh first images of the 3 cameras (WS)
        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_1']);
        img_1_t = imread(['frame_' int2str(frame_nr) '.bmp']);
        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_2']);
        img_2_t = imread(['frame_' int2str(frame_nr) '.bmp']);
        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_3']);
        img_3_t = imread(['frame_' int2str(frame_nr) '.bmp']);
        
        figure()
        hold on
        subplot(1,3,1); hold on
        h_im1 = imshow(img_1_t);
        hold off
        subplot(1,3,2); hold on
        h_im2 = imshow(img_2_t);
        hold off
        subplot(1,3,3); hold on
        h_im3 = imshow(img_3_t);
        hold off
        hold off
        %cycles through the frames until you get a set of images you
        %like(WS)
        while show_frames == 1
            
            prompt = 'Press enter to proceed to the next frame, press y+enter to select current frame.';
            str = input(prompt,'s');
            if isempty(str)==1
                frame_nr = frame_nr+1;
                
                cd([mov_folder '/mov_' int2str(mov_nr) '/cam_1']);
                img_1_t = imread(['frame_' int2str(frame_nr) '.bmp']);
                cd([mov_folder '/mov_' int2str(mov_nr) '/cam_2']);
                img_2_t = imread(['frame_' int2str(frame_nr) '.bmp']);
                cd([mov_folder '/mov_' int2str(mov_nr) '/cam_3']);
                img_3_t = imread(['frame_' int2str(frame_nr) '.bmp']);
                
                set(h_im1, 'CData', img_1_t);
                set(h_im2, 'CData', img_2_t);
                set(h_im3, 'CData', img_3_t);
                
                drawnow;
                
            elseif strcmp(str,'y')==1
                show_frames = 0;
                close all;
            else
                % do nothing
            end
        end
        
        disp(['frame nr ' int2str(frame_nr) ' has been selected.']);
        
        % Display the 3 different views of the selected frame in 3 figures
        % and ask the user to select the wing, neck and
        % abdomen hinge points:
        
        disp('Select the following points in the Camera 1 view');
        disp('Neck hinge (yellow)');
        disp('Abdomen hinge (green)');
        disp('Right wing hinge (cyan)');
        disp('Right wing tip (blue)');
        disp('Double click on the blue point once all points have been positioned correctly');
        
        figure('name','Camera 1')%select nech and wing hinge positions in each frame
        ax_1 = gca;
        hold on
        imshow(img_1_t)
        hold off
        h1_neck = impoint(ax_1);
        setColor(h1_neck,'y');
        h1_abdm = impoint(ax_1);
        setColor(h1_abdm,'g');
        h1_hinge_R = impoint(ax_1);
        setColor(h1_hinge_R,'c');
        h1_tip_R = impoint(ax_1);
        setColor(h1_tip_R,'b');
        wait(h1_tip_R);
        
        hinge_point_neck_1 = getPosition(h1_neck);
        hinge_point_abd_1 = getPosition(h1_abdm);
        hinge_point_wing_R_1 = getPosition(h1_hinge_R);
        wingtip_point_R_1 = getPosition(h1_tip_R);
        
        disp('Select the following points in the Camera 2 view');
        disp('Neck hinge (yellow)');
        disp('Abdomen hinge (green)');
        disp('Left wing hinge (mangenta)');
        disp('Left wing tip (red)');
        disp('Double click on the red point once all points have been positioned correctly');
        
        figure('name','Camera 2')
        ax_2 = gca;
        hold on
        imshow(img_2_t);
        hold off
        h2_neck = impoint(ax_2);
        setColor(h2_neck,'y');
        h2_abdm = impoint(ax_2);
        setColor(h2_abdm,'g');
        h2_hinge_L = impoint(ax_2);
        setColor(h2_hinge_L,'m');
        h2_tip_L = impoint(ax_2);
        setColor(h2_tip_L,'r');
        wait(h2_tip_L);
        
        hinge_point_neck_2 = getPosition(h2_neck);
        hinge_point_abd_2 = getPosition(h2_abdm);
        hinge_point_wing_L_2 = getPosition(h2_hinge_L);
        wingtip_point_L_2 = getPosition(h2_tip_L);
        
        disp('Select the following points in the Camera 3 view');
        disp('Neck hinge (yellow)');
        disp('Abdomen hinge (green)');
        disp('Left wing hinge (mangenta)');
        disp('Left wing tip (red)');
        disp('Right wing hinge (cyan)');
        disp('Right wing tip (blue)');
        disp('Double click on the blue point once all points have been positioned correctly');
        
        figure('name','Camera 3')
        ax_3 = gca;
        hold on
        imshow(img_3_t);
        hold off
        h3_neck = impoint(ax_3);
        setColor(h3_neck,'y');
        h3_abdm = impoint(ax_3);
        setColor(h3_abdm,'g');
        h3_hinge_L = impoint(ax_3);
        setColor(h3_hinge_L,'m');
        h3_tip_L = impoint(ax_3);
        setColor(h3_tip_L,'r');
        h3_hinge_R = impoint(ax_3);
        setColor(h3_hinge_R,'c');
        h3_tip_R = impoint(ax_3);
        setColor(h3_tip_R,'b');
        wait(h3_tip_R);
        
        hinge_point_neck_3 = getPosition(h3_neck);
        hinge_point_abd_3 = getPosition(h3_abdm);
        hinge_point_wing_L_3 = getPosition(h3_hinge_L);
        wingtip_point_L_3 = getPosition(h3_tip_L);
        hinge_point_wing_R_3 = getPosition(h3_hinge_R);
        wingtip_point_R_3 = getPosition(h3_tip_R);
        
        close all;
        
        % Calculate 3D locations by averaging the point locations:
        
        cd([mov_folder '/calibration']);
    
        load('cam_calib.mat');%loads the DLT parameters (WS)
        %385 and 373 are the offsets of the image. original image was
        %larger but was then cropped. I assume...(WS)
        neck_loc_1 = camera_to_world_projection(calib_par_cam_1,[385+hinge_point_neck_1(1); 373+hinge_point_neck_1(2); 1]);
        neck_loc_2 = camera_to_world_projection(calib_par_cam_2,[385+hinge_point_neck_2(1); 373+hinge_point_neck_2(2); 1]);
        neck_loc_3 = camera_to_world_projection(calib_par_cam_3,[385+hinge_point_neck_3(1); 373+hinge_point_neck_3(2); 1]);
        %convert 2D from 3 images to 3D. However points were selected by
        %user
        neck_loc = [(neck_loc_2(1)+neck_loc_3(1))/2; (neck_loc_1(2)+neck_loc_3(2))/2; (neck_loc_1(3)+neck_loc_2(3))/2];
        
        abdm_loc_1 = camera_to_world_projection(calib_par_cam_1,[385+hinge_point_abd_1(1); 373+hinge_point_abd_1(2); 1]);
        abdm_loc_2 = camera_to_world_projection(calib_par_cam_2,[385+hinge_point_abd_2(1); 373+hinge_point_abd_2(2); 1]);
        abdm_loc_3 = camera_to_world_projection(calib_par_cam_3,[385+hinge_point_abd_3(1); 373+hinge_point_abd_3(2); 1]);
        
        abdm_loc = [(abdm_loc_2(1)+abdm_loc_3(1))/2; (abdm_loc_1(2)+abdm_loc_3(2))/2; (abdm_loc_1(3)+abdm_loc_2(3))/2];
        
        wing_hinge_loc_L_2 = camera_to_world_projection(calib_par_cam_2,[385+hinge_point_wing_L_2(1); 373+hinge_point_wing_L_2(2); 1]);
        wing_hinge_loc_L_3 = camera_to_world_projection(calib_par_cam_3,[385+hinge_point_wing_L_3(1); 373+hinge_point_wing_L_3(2); 1]);
        
        wing_hinge_loc_L = [(wing_hinge_loc_L_2(1)+wing_hinge_loc_L_3(1))/2; wing_hinge_loc_L_3(2); wing_hinge_loc_L_2(3)];
        
        wingtip_point_L_2 = camera_to_world_projection(calib_par_cam_2,[385+wingtip_point_L_2(1); 373+wingtip_point_L_2(2); 1]);
        wingtip_point_L_3 = camera_to_world_projection(calib_par_cam_3,[385+wingtip_point_L_3(1); 373+wingtip_point_L_3(2); 1]);
        %385 and 373 are added since calibration was performed on a larger
        wing_tip_loc_L = [(wingtip_point_L_2(1)+wingtip_point_L_3(1))/2; wingtip_point_L_3(2); wingtip_point_L_2(3)];
        
        wing_hinge_loc_R_1 = camera_to_world_projection(calib_par_cam_1,[385+hinge_point_wing_R_1(1); 373+hinge_point_wing_R_1(2); 1]);
        wing_hinge_loc_R_3 = camera_to_world_projection(calib_par_cam_3,[385+hinge_point_wing_R_3(1); 373+hinge_point_wing_R_3(2); 1]);
        
        wing_hinge_loc_R = [wing_hinge_loc_R_3(1); (wing_hinge_loc_R_1(2)+wing_hinge_loc_R_3(2))/2; wing_hinge_loc_R_1(3)];
        
        wingtip_point_R_1 = camera_to_world_projection(calib_par_cam_1,[385+wingtip_point_R_1(1); 373+wingtip_point_R_1(2); 1]);
        wingtip_point_R_3 = camera_to_world_projection(calib_par_cam_3,[385+wingtip_point_R_3(1); 373+wingtip_point_R_3(2); 1]);
        
        wing_tip_loc_R = [wingtip_point_R_3(1); (wingtip_point_R_1(2)+wingtip_point_R_3(2))/2; wingtip_point_R_1(3)];
        
        % Calculate wing lengths:
        
        wing_length_L = sqrt((wing_tip_loc_L(1)-wing_hinge_loc_L(1))^2+(wing_tip_loc_L(2)-wing_hinge_loc_L(2))^2+(wing_tip_loc_L(3)-wing_hinge_loc_L(3))^2);
        wing_length_R = sqrt((wing_tip_loc_R(1)-wing_hinge_loc_R(1))^2+(wing_tip_loc_R(2)-wing_hinge_loc_R(2))^2+(wing_tip_loc_R(3)-wing_hinge_loc_R(3))^2);
        
        % Calculate strokeplane orientation:
        
        wing_hinge_center = mean([wing_hinge_loc_L wing_hinge_loc_R],2);
        
        % Construct body reference frame (x = longitudinal axis, y = line
        % between the two hinges, z = normal to x & y (dorsal side is
        % positive)):
        X_body = neck_loc-abdm_loc;
        X_body = X_body/norm(X_body);
        Y_body = wing_hinge_loc_L-wing_hinge_center;
        Y_body = Y_body/norm(Y_body);
        Z_body = cross(X_body,Y_body);
        Z_body = Z_body/norm(Z_body);
        
        % Calculate the Tait-Bryan angles which rotate from world to body
        % reference frame:
        %q-body is the quaternion rotation
        [R_body, q_body] = find_rotation([X_body Y_body Z_body],eye(3));
        
%         %alpha = asin(X_body(2)/sqrt(1-X_body(3)^2)); % yaw rotation
%         alpha = atan2(X_body(2),X_body(1));          % yaw rotation
%         beta  = asin(-X_body(3));                    % pitch rotation
%         
%         if cos(beta)>0.5
%             gamma = atan2(-Z_body(2),(-Z_body(1)/cos(beta)));
%         elseif cos(beta)<=0.5
%             gamma = atan2(-Z_body(2),(Z_body(3)/sin(beta)));
%         end
%         
%         R_body = [cos(alpha)*cos(beta) cos(alpha)*sin(beta)*sin(gamma)-cos(gamma)*sin(alpha) sin(alpha)*sin(gamma)+cos(alpha)*cos(gamma)*sin(beta); ...
%                   cos(beta)*sin(alpha) cos(alpha)*cos(gamma)+sin(alpha)*sin(beta)*sin(gamma) cos(gamma)*sin(alpha)*sin(beta)-cos(alpha)*sin(gamma); ...
%                   -sin(beta) cos(beta)*sin(gamma) cos(beta)*cos(gamma)];
%         
%         q_body = mat2quat(R_body);
        
        % Calculate strokeplane reference frame by rotating the body
        % reference frame 45 degrees w.r.t. the body-Y-axis.
        
        beta_strk = -55/180*pi;

        R_strk = [ cos(beta_strk)  0 -sin(beta_strk); ...
                   0               1 0             ; ...
                   sin(beta_strk)  0 cos(beta_strk) ];
        
        q_strk = mat2quat(R_strk);
        
        % Scale the fly wings to wing lengths:
        
        % Save the model:
        model = {};
        model.neck_hinge = neck_loc;
        model.abdm_hinge = abdm_loc;
        model.wing_hinge_L = wing_hinge_loc_L;
        model.wing_hinge_R = wing_hinge_loc_R;
        model.wing_hinge_center = wing_hinge_center;
        model.wing_tip_L = wing_tip_loc_L;
        model.wing_tip_R = wing_tip_loc_R;
        model.wing_length_L = wing_length_L;
        model.wing_length_R = wing_length_R;
        model.R_body = R_body;
        model.q_body = q_body;
        model.R_strk = R_strk;
        model.q_strk = q_strk;
        
        try
            cd([mov_folder '/fly_model']);
            save('fly_model.mat','model');
        catch
            cd(mov_folder);
            mkdir('fly_model');
            save('fly_model.mat','model');
        end
        
    end

end

function masks = user_input_masks(mov_folder, mov_nr)
%this function finds the mask of the videos (for the tether and body)
%based on the code, the user can either load a predefined mask or create a
%new one
    try
        cd([mov_folder '/fly_model']);
        load('img_masks.mat');
    catch
        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_1']);
        img_1_t = imread('frame_0.bmp');
        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_2']);
        img_2_t = imread('frame_0.bmp');
        cd([mov_folder '/mov_' int2str(mov_nr) '/cam_3']);
        img_3_t = imread('frame_0.bmp');
        
%this portion is used to select the mask around the tether and determines
%the image ROI
        disp('select ROI image 1, close ROI by double clicking on start point');
        [~, ~, mask_img_1, ~, ~] = roipoly(img_1_t); %reurns the ROI as a white patch. Functions allows us to select a polygon ROI
        disp('select mask region around the tether, close ROI by double clicking on start point');
        [~, ~, mask_tether_1, ~, ~] = roipoly(im2double(img_1_t).*mask_img_1);
        disp('select ROI image 2, close ROI by double clicking on start point');
        [~, ~, mask_img_2, ~, ~] = roipoly(img_2_t);
        disp('select mask region around the tether, close ROI by double clicking on start point');
        [~, ~, mask_tether_2, ~, ~] = roipoly(im2double(img_2_t).*mask_img_2);
        disp('select ROI image 3, close ROI by double clicking on start point');
        [~, ~, mask_img_3, ~, ~] = roipoly(img_3_t);
        disp('select mask region around the tether, close ROI by double clicking on start point');
        [~, ~, mask_tether_3, ~, ~] = roipoly(im2double(img_3_t).*mask_img_3);
%this part selects the mask around the body for the three different views
        disp('select body mask image 1');
        [~, ~, mask_body_1, ~, ~] = roipoly(img_1_t);
        disp('select body mask image 2');
        [~, ~, mask_body_2, ~, ~] = roipoly(img_2_t);
        disp('select body mask image 3');
        [~, ~, mask_body_3, ~, ~] = roipoly(img_3_t);

        disp('select left image w.r.t. symmetry axis')
        [~, ~, im_left, ~, ~] = roipoly(img_3_t);
%moves into the folder that has the masks and saves them there
        cd(mov_folder);
        mkdir('fly_model');
        cd([mov_folder '/fly_model']);
        save('img_masks.mat','mask_img_1','mask_img_2','mask_img_3','mask_tether_1','mask_tether_2','mask_tether_3','mask_body_1','mask_body_2','mask_body_3','im_left');

    end
    
    masks = cell(3,1);
    masks{1}.img_mask       = mask_img_1;
    masks{1}.body_mask      = mask_body_1;
    masks{1}.tether_mask    = mask_tether_1;
    masks{2}.img_mask       = mask_img_2;
    masks{2}.body_mask      = mask_body_2;
    masks{2}.tether_mask    = mask_tether_2;
    masks{3}.img_mask       = mask_img_3;
    masks{3}.body_mask      = mask_body_3;
    masks{3}.tether_mask    = mask_tether_3;
    masks{3}.im_left_side   = im_left;

end

function [vox_1, vox_2, e1, e2, c1, c2] = iterative_svd_w_dropout(voxels,n_iter)
    tic
    for i = 1:n_iter
        if i == 1
            [e1, c1] = fit_plane(voxels);
            [e2, c2] = fit_plane(voxels);
            vox_1 = drop_points(voxels,e1,c1,1.0,2);
            vox_2 = drop_points(voxels,e2,c2,1.0,3);
        else
            [e1, c1] = fit_plane(vox_1);
            [e2, c2] = fit_plane(vox_2);
            vox_1 = drop_points(vox_1,e1,c1,2.0,3);
            vox_2 = drop_points(vox_2,e2,c2,2.0,3);
        end
    end
    toc
end

function vox_new = drop_points(vox,e,c,tresh,ax)

    A = bsxfun(@minus,vox,c);
    e_dist_1 = (A'*e(:,1));
    e_dist_2 = (A'*e(:,2));
    e_dist_3 = (A'*e(:,3));
    
    if ax == 1
        vox_new = vox(:,abs(e_dist_1)<=tresh);
    elseif ax == 2
        vox_new = vox(:,abs(e_dist_2)<=tresh);
    elseif ax == 3
        vox_new = vox(:,abs(e_dist_3)<=tresh);
    end

end

% function vox_new = drop_points(vox,e,c,std_tresh,ax)
% 
%     A = bsxfun(@minus,vox,c);
%     e_dist_1 = (A'*e(:,1));
%     e_dist_2 = (A'*e(:,2));
%     e_dist_3 = (A'*e(:,3));
%     if ax == 2 && size(vox,2) > 100
%         %e_slct = (e_dist_2.^2).*sqrt(e_dist_2.^2+e_dist_3.^2);
%         %e_slct = e_dist_2.^2;
%         e_slct = sqrt(e_dist_1.^2+e_dist_3.^2);
%         vox_new = vox(:,e_slct>(std_tresh*std(e_slct)));
%     elseif ax == 3 && size(vox,2) > 100
%         %e_slct = (e_dist_3.^2).*sqrt(e_dist_2.^2+e_dist_3.^2);
%         %e_slct = e_dist_3.^2;
%         e_slct = sqrt(e_dist_1.^2+e_dist_2.^2);
%         vox_new = vox(:,e_slct>(std_tresh*std(e_slct)));
%     else
%         vox_new = vox;
%         disp('100 points');
%     end
% 
% end

function [e, c] = fit_plane(vox)

    c = mean(vox,2);
    A = bsxfun(@minus,vox,c);
    [U,~,~] = svds(A);
    
    e1 = U(:,1);
    e2 = U(:,2);
    e3 = U(:,3);
    
    %e1 = e1/norm(e1);
    %e2 = e2/norm(e2);
    %e3 = e3/norm(e3);
    
    e = [e1 e2 e3];

end

function [R, q] = find_rotation(E_W,E_B)
    
    % Find the rotation matrix and quaternion using SVD, matrix E_W is the
    % orientation of the body reference frame in the world reference frame
    % whilst E_B is the orientation of tye body reference frame in the body
    % reference frame (eye(3)):
    
    B = zeros(3,3);
    
    for i = 1:3
       B = B+E_W(:,i)*E_B(:,i)' ;
    end
    %check wahba's problem amd its solution by singular value decomposition
    %(WS)
   
    [U,S,V] = svd(B);
    
    M = diag([1 1 det(U)*det(V)]);
    
    R = U*M*V';
    
    q = mat2quat(R);

end

function R_w = return_rot_mat(calib_param)

    theta = sqrt(calib_param(4)^2+calib_param(5)^2+calib_param(6)^2);
    omega = [ 0 -calib_param(6) calib_param(5); ...
              calib_param(6) 0 -calib_param(4); ...
              -calib_param(5) calib_param(4) 0 ];
    R_w = eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);

end

function X_w = camera_to_world_projection(calib_param,X_uv)

    C = [calib_param(1) calib_param(3) 0 0; ...
         0              calib_param(2) 0 0; ...
         0              0              0 1];
     
    theta = sqrt(calib_param(4)^2+calib_param(5)^2+calib_param(6)^2);
    omega = [ 0 -calib_param(6) calib_param(5); ...
              calib_param(6) 0 -calib_param(4); ...
              -calib_param(5) calib_param(4) 0 ];
    R = eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);
    
    t = [calib_param(7); calib_param(8); calib_param(9)];
    
    K_inv = inv([R t; 0 0 0 1]);
    
    C_inv = pinv(C);
    
    X_i = X_uv;
    
    X_w = K_inv*C_inv*X_i;

end

function X_uv = world_to_camera_projection(calib_param,X_w)

    C = [calib_param(1) calib_param(3) 0 0; ...
         0              calib_param(2) 0 0; ...
         0              0              0 1];
     
    theta = sqrt(calib_param(4)^2+calib_param(5)^2+calib_param(6)^2);
    omega = [ 0 -calib_param(6) calib_param(5); ...
              calib_param(6) 0 -calib_param(4); ...
              -calib_param(5) calib_param(4) 0 ];
    R = eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);
    
    t = [calib_param(7); calib_param(8); calib_param(9)];
    
    K = [R t; 0 0 0 1];
    
    X_i = C*K*X_w;
    
    X_uv = X_i;

end

function q = quat_subtraction(q1,q2)

    % Perform a quaternion "subtraction" by multiplying q1 with the complex
    % conjugate of q2:
    
    q2_conj = [q2(1); -q2(2); -q2(3); -q2(4)];
    
    q = quat_multiply(q1,q2_conj);

end

function q = quat_multiply(q1,q2)

    % Perform a quaternion multiplication:
    
    Q1 = [q1(1) -q1(2) -q1(3) -q1(4); ...
          q1(2) q1(1) -q1(4) q1(3); ...
          q1(3) q1(4) q1(1) -q1(2); ...
          q1(4) -q1(3) q1(2) q1(1)];
    
    q = Q1*q2;
    
    q = q/norm(q);

end

function R = quat2mat(q)

    q = q/norm(q);
    
    R = [ 2*q(1)^2-1+2*q(2)^2 2*q(2)*q(3)+2*q(1)*q(4) 2*q(2)*q(4)-2*q(1)*q(3); ...
          2*q(2)*q(3)-2*q(1)*q(4) 2*q(1)^2-1+2*q(3)^2 2*q(3)*q(4)+2*q(1)*q(2); ...
          2*q(2)*q(4)+2*q(1)*q(3) 2*q(3)*q(4)-2*q(1)*q(2) 2*q(1)^2-1+2*q(4)^2 ];

end

function q = mat2quat(R)

    q0 = 0.5*sqrt(R(1,1)+R(2,2)+R(3,3)+1);
    q1 = (R(2,3)-R(3,2))/(4*q0);
    q2 = (R(3,1)-R(1,3))/(4*q0);
    q3 = (R(1,2)-R(2,1))/(4*q0);
    
    % Normalize the quaternion:
    
    q = [q0; q1; q2; q3]/norm([q0; q1; q2; q3]);

%     % Convert a rotation matrix to a quaternion:
%     
%     if trace(R) > 0
%         S = sqrt(R(1,1)+R(2,2)+R(3,3)+1);
%         q0 = S/2;
%         q1 = (R(2,3)-R(3,2))/(2*S);
%         q2 = (R(3,1)-R(1,3))/(2*S);
%         q3 = (R(1,2)-R(2,1))/(2*S);
%     elseif (R(1,1)>R(2,2))&&(R(1,1)>R(3,3))
%         S = sqrt(R(1,1)-R(2,2)-R(3,3)+1);
%         q0 = (R(2,3)-R(3,2))/(2*S);
%         q1 = S/2;
%         q2 = (R(1,2)+R(2,1))/(2*S);
%         q3 = (R(1,3)+R(3,1))/(2*S);
%     elseif (R(2,2)>R(3,3))
%         S = sqrt(-R(1,1)+R(2,2)-R(3,3)+1);
%         q0 = (R(3,2)-R(1,3))/(2*S);
%         q1 = (R(1,2)+R(2,1))/(2*S);
%         q2 = S/2;
%         q3 = (R(2,3)+R(3,2))/(2*S);
%     else
%         S = sqrt(-R(1,1)-R(2,2)+R(3,3)+1);
%         q0 = (R(1,2)-R(2,1))/(2*S);
%         q1 = (R(1,3)+R(3,1))/(2*S);
%         q2 = (R(2,3)+R(3,2))/(2*S);
%         q3 = S/2;
%     end
% 
%     % Normalize the quaternion:
%     
%     q = [q0; q1; q2; q3]/norm([q0; q1; q2; q3]);
end