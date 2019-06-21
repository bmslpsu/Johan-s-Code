function find_focal_grid(calibration_loc,wing_hinge_center)
%consider replacing this with image center (WS)
% can probably remove this input to the function(WS)
    close all;

    % Reconstruct the focal grid using the calibration file:
    
    cd(calibration_loc);
    
    load('cam_calib.txt'); %loads the calib_par_cam_(1,2,3)
    calib_par_cam_1=cam_calib(:,1);
    calib_par_cam_2=cam_calib(:,2);
    calib_par_cam_3=cam_calib(:,3);
    % Project the uv grid of 256 by 256 pixels to the world reference frame
    %it seems the original image is either cropped or calibration was
    %performed on a larger image (WS)
    [U_grid, V_grid] = meshgrid(385:640,373:628);
    %for my apploication can use 1:end for length and width of image(WS)
    
    X_uv = [reshape(U_grid,1,256^2); reshape(V_grid,1,256^2); ones(1,256^2)]; %converts into a 3xsomething matrix containing the coordinates of each point in the image
    
    XW_cam_1 = camera_to_world_projection(calib_par_cam_1,X_uv); 
    XW_cam_2 = camera_to_world_projection(calib_par_cam_2,X_uv);
    XW_cam_3 = camera_to_world_projection(calib_par_cam_3,X_uv);
%some rotations, need clarifications (WS)
    RW_1 = return_rot_mat(calib_par_cam_1); %rotation from camera to world reference frame
    RW_2 = return_rot_mat(calib_par_cam_2);
    RW_3 = return_rot_mat(calib_par_cam_3);
    
    e1 = RW_1'*[0; 0; -1]; %3 axis of the world reference frame
    e2 = RW_2'*[0; 0; -1];%vector along which i am extruding the image
    e3 = RW_3'*[0; 0; -1];
    
%     TW_cam_1 = camera_to_world_projection(calib_par_cam_1,[385+im_cent_cam_1(1); 373+im_cent_cam_1(2); 1]);
%     TW_cam_2 = camera_to_world_projection(calib_par_cam_2,[385+im_cent_cam_2(1); 373+im_cent_cam_2(2); 1]);
%     TW_cam_3 = camera_to_world_projection(calib_par_cam_3,[385+im_cent_cam_3(1); 373+im_cent_cam_3(2); 1]);
    
    P1 = XW_cam_1(1:3,:);
    P2 = XW_cam_2(1:3,:);
    P3 = XW_cam_3(1:3,:);
    %P_vox coordinates of the voxels 
    %P_vox_neighbors coordinates :helps disregard certain voxels
    [vox_1_mat, vox_2_mat, vox_3_mat, P_vox, P_vox_neighbor] = calculate_focal_grid(e1,e2,e3,P1,P2,P3,wing_hinge_center);
    %vox matricies, need a reference(WS)
    whos vox_1_mat
    whos vox_2_mat
    whos vox_3_mat
    
    % Save the sparse matrices:
    
    cd(calibration_loc);
    
    save('Focal_grid.mat','vox_1_mat','vox_2_mat','vox_3_mat', 'P_vox', 'P_vox_neighbor');
    
end

function [vox_1_mat, vox_2_mat, vox_3_mat, P_vox, P_vox_neighbor] = calculate_focal_grid(e1,e2,e3,P1,P2,P3,wing_hinge_center)
    
    voxel_size = 0.040; 

    % Step 1: project points P1, P2 and P3 on the XZ, YZ and XY plane
    % respectively:
    
    P1_YZ = P1-e1*(P1(1,:)./e1(1));
    P2_XZ = P2-e2*(P2(2,:)./e2(2));
    P3_XY = P3-e3*(P3(3,:)./e3(3));
    
    %D1_YZ = delaunayn([P1_YZ(2,:)' P1_YZ(3,:)']);
    %D2_XZ = delaunayn([P2_XZ(1,:)' P2_XZ(3,:)']);
    %D3_XY = delaunayn([P3_XY(1,:)' P3_XY(2,:)']);
    
    % Step 2: create a grid of 513^3 voxels of 40 by 40 by 40 micrometers
    %wing hinge center is the center of the line formed by the left and right wing hinge
    X_grid_trans = wing_hinge_center(1);
    Y_grid_trans = wing_hinge_center(2);
    Z_grid_trans = wing_hinge_center(3);
    
%     X_grid_trans = 0;
%     Y_grid_trans = 0;
%     Z_grid_trans = 0;
    
    [X_grid, Y_grid, Z_grid] = meshgrid((-128*voxel_size+X_grid_trans):voxel_size:(128*voxel_size+X_grid_trans), ...
                                        (-128*voxel_size+Y_grid_trans):voxel_size:(128*voxel_size+Y_grid_trans), ...
                                        (-128*voxel_size+Z_grid_trans):voxel_size:(128*voxel_size+Z_grid_trans));
    
    % Step 3: iterate through all voxels in the grid and project along e1,
    % e2 and e3 on the XZ, YZ and XY plane respectively. Find the closest
    % corresponding P1_XZ, P2_YZ and P3_XY triplets and add a one at these
    % point indices in the sparse image to voxel transformation matrix.
    
    P_vox = [reshape(X_grid,1,257^3); reshape(Y_grid,1,257^3); reshape(Z_grid,1,257^3)]; % consider replacing 257 with size of each array(WS)
    
    P_vox_YZ = P_vox-e1*(P_vox(1,:)./e1(1));
    P_vox_XZ = P_vox-e2*(P_vox(2,:)./e2(2));
    P_vox_XY = P_vox-e3*(P_vox(3,:)./e3(3));
    
    ind = 1:(257^3);
    
    P_vox_neighbor = [ind; ind-1; ind+1; ind-257; ind+257; ind-257^2; ind+257^2]; % no idea what this for as of yet, seems to have the indicies of the voxels
    %with some variations within the rows
    
    % Create a grid:
    dy_1 = (P1_YZ(2,256^2)-P1_YZ(2,256))/255;
    dz_1 = (P1_YZ(3,1)-P1_YZ(3,256))/255;
    dz_dy_1 = (P1_YZ(3,256^2)-P1_YZ(3,256))/(P1_YZ(2,256^2)-P1_YZ(2,256));
    dy_dz_1 = (P1_YZ(2,1)-P1_YZ(2,256))/(P1_YZ(3,1)-P1_YZ(3,256));
    
    dx_2 = (P2_XZ(1,256)-P2_XZ(1,256^2))/255;
    dz_2 = (P2_XZ(3,1)-P2_XZ(3,256))/255;
    dz_dx_2 = (P2_XZ(3,256)-P2_XZ(3,256^2))/(P2_XZ(1,256)-P2_XZ(1,256^2));
    dx_dz_2 = (P2_XZ(1,1)-P2_XZ(1,256))/(P2_XZ(3,1)-P2_XZ(3,256));
    
    dx_3 = (P3_XY(1,1)-P3_XY(1,256))/255;
    dy_3 = (P3_XY(2,256)-P3_XY(2,256^2))/255;
    dy_dx_3 = (P3_XY(2,1)-P3_XY(2,256))/(P3_XY(1,1)-P3_XY(1,256));
    dx_dy_3 = (P3_XY(1,256)-P3_XY(1,256^2))/(P3_XY(2,256)-P3_XY(2,256^2));
    
%     dx_A3 = (P3_XY(1,256)-P3_XY(1,256^2))/510;
%     dy_A3 = (P3_XY(2,256)-P3_XY(2,256^2))/510;
%     dx_B3 = (P3_XY(1,1)-P3_XY(1,256))/510;
%     dy_B3 = (P3_XY(2,1)-P3_XY(2,256))/510;
%     dx_C3 = -(P3_XY(1,256)-P3_XY(1,256^2))/510;
%     dy_C3 = -(P3_XY(2,256)-P3_XY(2,256^2))/510;
%     dx_D3 = -(P3_XY(1,1)-P3_XY(1,256))/510;
%     dy_D3 = -(P3_XY(2,1)-P3_XY(2,256))/510;
    
    P1_grid = [reshape(P1_YZ(2,:),256^2,1) reshape(P1_YZ(3,:),256^2,1) ...
               reshape(P1_YZ(2,:),256^2,1)+dy_dz_1*(dz_1/2) reshape(P1_YZ(3,:),256^2,1)+dz_1/2 ...
               reshape(P1_YZ(2,:),256^2,1)+dy_1/2 reshape(P1_YZ(3,:),256^2,1)+dz_dy_1*(dy_1/2) ...
               reshape(P1_YZ(2,:),256^2,1)-dy_dz_1*(dz_1/2) reshape(P1_YZ(3,:),256^2,1)-dz_1/2 ...
               reshape(P1_YZ(2,:),256^2,1)-dy_1/2 reshape(P1_YZ(3,:),256^2,1)-dz_dy_1*(dy_1/2)];
    P2_grid = [reshape(P2_XZ(1,:),256^2,1) reshape(P2_XZ(3,:),256^2,1) ...
               reshape(P2_XZ(1,:),256^2,1)+dx_dz_2*(dz_2/2) reshape(P2_XZ(3,:),256^2,1)+dz_2/2 ...
               reshape(P2_XZ(1,:),256^2,1)+dx_2/2 reshape(P2_XZ(3,:),256^2,1)+dz_dx_2*(dx_2/2) ...
               reshape(P2_XZ(1,:),256^2,1)-dx_dz_2*(dz_2/2) reshape(P2_XZ(3,:),256^2,1)-dz_2/2 ...
               reshape(P2_XZ(1,:),256^2,1)-dx_2/2 reshape(P2_XZ(3,:),256^2,1)-dz_dx_2*(dx_2/2)];
    P3_grid = [reshape(P3_XY(1,:),256^2,1) reshape(P3_XY(2,:),256^2,1) ...
               reshape(P3_XY(1,:),256^2,1)+dx_dy_3*(dy_3/2) reshape(P3_XY(2,:),256^2,1)+dy_3/2 ...
               reshape(P3_XY(1,:),256^2,1)+dx_3/2 reshape(P3_XY(2,:),256^2,1)+dy_dx_3*(dx_3/2) ...
               reshape(P3_XY(1,:),256^2,1)-dx_dy_3*(dy_3/2) reshape(P3_XY(2,:),256^2,1)-dy_3/2 ...
               reshape(P3_XY(1,:),256^2,1)-dx_3/2 reshape(P3_XY(2,:),256^2,1)-dy_dx_3*(dx_3/2)];
    
    P_YZ_y = [min(P1_grid(:,1)) max(P1_grid(:,1))];
    P_YZ_z = [min(P1_grid(:,2)) max(P1_grid(:,2))];
    
    TA1 = P1_grid(1,4)+dz_dy_1*(P_YZ_y-P1_grid(1,3));
    TC1 = P1_grid(1,8)+dz_dy_1*(P_YZ_y-P1_grid(1,7));
    TB1 = P1_grid(1,5)+dy_dz_1*(P_YZ_z-P1_grid(1,6));
    TD1 = P1_grid(1,9)+dy_dz_1*(P_YZ_z-P1_grid(1,10));
    
    P_XZ_x = [min(P2_grid(:,1)) max(P2_grid(:,1))];
    P_XZ_z = [min(P2_grid(:,2)) max(P2_grid(:,2))];
    
    TA2 = P2_grid(1,4)+dz_dx_2*(P_XZ_x-P2_grid(1,3));
    TC2 = P2_grid(1,8)+dz_dx_2*(P_XZ_x-P2_grid(1,7));
    TB2 = P2_grid(1,5)+dx_dz_2*(P_XZ_z-P2_grid(1,6));
    TD2 = P2_grid(1,9)+dx_dz_2*(P_XZ_z-P2_grid(1,10));
    
    P_XY_x = [min(P3_grid(:,1)) max(P3_grid(:,1))];
    P_XY_y = [min(P3_grid(:,2)) max(P3_grid(:,2))];
    
    TA3 = P3_grid(1,4)+dy_dx_3*(P_XY_x-P3_grid(1,3));
    TC3 = P3_grid(1,8)+dy_dx_3*(P_XY_x-P3_grid(1,7));
    TB3 = P3_grid(1,5)+dx_dy_3*(P_XY_y-P3_grid(1,6));
    TD3 = P3_grid(1,9)+dx_dy_3*(P_XY_y-P3_grid(1,10));
    
    vox_1_select = [];
    vox_2_select = [];
    vox_3_select = [];
    
    
    tic
    for i = 1:256
        
        i
        
        TB1 = P1_grid((i-1)*256+1,5)+dy_dz_1*(P_vox_YZ(3,:)-P1_grid((i-1)*256+1,6));
        TD1 = P1_grid((i-1)*256+1,9)+dy_dz_1*(P_vox_YZ(3,:)-P1_grid((i-1)*256+1,10));
        
        TB2 = P2_grid((i-1)*256+1,5)+dx_dz_2*(P_vox_XZ(3,:)-P2_grid((i-1)*256+1,6));
        TD2 = P2_grid((i-1)*256+1,9)+dx_dz_2*(P_vox_XZ(3,:)-P2_grid((i-1)*256+1,10));
        
        TA3 = P3_grid((i-1)*256+1,4)+dy_dx_3*(P_vox_XY(1,:)-P3_grid((i-1)*256+1,3));
        TC3 = P3_grid((i-1)*256+1,8)+dy_dx_3*(P_vox_XY(1,:)-P3_grid((i-1)*256+1,7));

        col_select_1 = ind((P_vox_YZ(2,:)<=TB1)&(P_vox_YZ(2,:)>TD1));
        col_select_2 = ind((P_vox_XZ(1,:)<=TB2)&(P_vox_XZ(1,:)>TD2));
        row_select_3 = ind((P_vox_XY(2,:)<=TA3)&(P_vox_XY(2,:)>TC3));
        
        for j = 1:256
            
            if isempty(col_select_1) == 0
                
                TA1 = P1_grid((i-1)*256+j,4)+dz_dy_1*(P_vox_YZ(2,col_select_1)-P1_grid((i-1)*256+j,3));
                TC1 = P1_grid((i-1)*256+j,8)+dz_dy_1*(P_vox_YZ(2,col_select_1)-P1_grid((i-1)*256+j,7));
                
                vox_t_1 = ((P_vox_YZ(3,col_select_1)<=TA1)&(P_vox_YZ(3,col_select_1)>TC1)).*col_select_1;
                vox_t_1(vox_t_1==0) = [];
                
                if isempty(vox_t_1) == 0
                    
                    vox_1_select = [vox_1_select [vox_t_1; ((i-1)*256+j)*ones(1,length(vox_t_1))]];
                    
                end
                
            end
            
            if isempty(col_select_2) == 0
            
                TA2 = P2_grid((i-1)*256+j,4)+dz_dx_2*(P_vox_XZ(1,col_select_2)-P2_grid((i-1)*256+j,3));
                TC2 = P2_grid((i-1)*256+j,8)+dz_dx_2*(P_vox_XZ(1,col_select_2)-P2_grid((i-1)*256+j,7));
                
                vox_t_2 = ((P_vox_XZ(3,col_select_2)<=TA2)&(P_vox_XZ(3,col_select_2)>TC2)).*col_select_2;
                vox_t_2(vox_t_2==0) = [];
                
                if isempty(vox_t_2) == 0
                    
                    vox_2_select = [vox_2_select [vox_t_2; ((i-1)*256+j)*ones(1,length(vox_t_2))]];

                end
                
            end
            
            if isempty(row_select_3) == 0
                
                TB3 = P3_grid((i-1)*256+j,5)+dx_dy_3*(P_vox_XY(2,row_select_3)-P3_grid((i-1)*256+j,6));
                TD3 = P3_grid((i-1)*256+j,9)+dx_dy_3*(P_vox_XY(2,row_select_3)-P3_grid((i-1)*256+j,10));
                
                vox_t_3 = ((P_vox_XY(1,row_select_3)<=TB3)&(P_vox_XY(1,row_select_3)>TD3)).*row_select_3;
                vox_t_3(vox_t_3==0) = [];
                
                if isempty(vox_t_3) == 0
                    
                    vox_3_select = [vox_3_select [vox_t_3; ((i-1)*256+j)*ones(1,length(vox_t_3))]];
                    
                end
                
            end
            
        end
        
    end
    toc
    
    size(vox_1_select)
    size(vox_2_select)
    size(vox_3_select)
    
    size(unique(vox_1_select(1,:)))
    
    whos vox_1_select
    whos vox_2_select
    whos vox_3_select
    
    figure()
    hold on
    plot(P1_grid(1,:),P1_grid(2,:),'.','Color','k')
    plot(P_vox_YZ(2,vox_1_select(1,:)),P_vox_YZ(3,vox_1_select(1,:)),'.','Color','r')
    hold off
    axis equal
    
    figure()
    hold on
    plot(P2_grid(1,:),P2_grid(2,:),'.','Color','k')
    plot(P_vox_XZ(1,vox_2_select(1,:)),P_vox_XZ(3,vox_2_select(1,:)),'.','Color','r')
    hold off
    axis equal
    
    figure()
    hold on
    plot(P3_grid(1,:),P3_grid(2,:),'.','Color','k')
    plot(P_vox_XY(1,vox_3_select(1,:)),P_vox_XY(2,vox_3_select(1,:)),'.','Color','r')
    hold off
    axis equal
    
%     figure()
%     hold on
%     plot(P1_grid(:,1),P1_grid(:,2),'.','Color','k')
%     plot(P_vox_YZ(2,:),P_vox_YZ(3,:),'.','Color','r')
%     plot(P1_grid(:,3),P1_grid(:,4),'o','Color','m')
%     plot(P1_grid(:,5),P1_grid(:,6),'o','Color','y')
%     plot(P1_grid(:,7),P1_grid(:,8),'o','Color','c')
%     plot(P1_grid(:,9),P1_grid(:,10),'o','Color','k')
%     plot(P_YZ_y,TA1,'Color','b')
%     plot(P_YZ_y,TC1,'Color','b')
%     plot(TB1,P_YZ_z,'Color','g')
%     plot(TD1,P_YZ_z,'Color','g')
%     hold off
%     axis equal
    
%     figure()
%     hold on
%     plot(P2_grid(:,1),P2_grid(:,2),'.','Color','k')
%     plot(P_vox_XZ(1,:),P_vox_XZ(3,:),'.','Color','r')
%     plot(P2_grid(:,3),P2_grid(:,4),'o','Color','m')
%     plot(P2_grid(:,5),P2_grid(:,6),'o','Color','y')
%     plot(P2_grid(:,7),P2_grid(:,8),'o','Color','c')
%     plot(P2_grid(:,9),P2_grid(:,10),'o','Color','k')
%     plot(P_XZ_x,TA2,'Color','b')
%     plot(P_XZ_x,TC2,'Color','b')
%     plot(TB2,P_XZ_z,'Color','g')
%     plot(TD2,P_XZ_z,'Color','g')
%     hold off
%     axis equal

%     figure()
%     hold on
%     plot(P3_grid(:,1),P3_grid(:,2),'.','Color','k')
%     plot(P_vox_XY(1,:),P_vox_XY(2,:),'.','Color','r')
%     plot(P3_grid(:,3),P3_grid(:,4),'o','Color','m')
%     plot(P3_grid(:,5),P3_grid(:,6),'o','Color','y')
%     plot(P3_grid(:,7),P3_grid(:,8),'o','Color','c')
%     plot(P3_grid(:,9),P3_grid(:,10),'o','Color','k')
%     plot(P_XY_x,TA3,'Color','b')
%     plot(P_XY_x,TC3,'Color','b')
%     plot(TB3,P_XY_y,'Color','g')
%     plot(TD3,P_XY_y,'Color','g')
%     hold off
%     axis equal
    
    
    % Now convert the indices to sparse matrices:
    tic
    
    vox_1_mat = sparse(vox_1_select(1,:),vox_1_select(2,:),ones(1,size(vox_1_select,2)),257^3,256^2);
    vox_2_mat = sparse(vox_2_select(1,:),vox_2_select(2,:),ones(1,size(vox_2_select,2)),257^3,256^2);
    vox_3_mat = sparse(vox_3_select(1,:),vox_3_select(2,:),ones(1,size(vox_3_select,2)),257^3,256^2);
    
    sum_vox_mat = sum(vox_1_mat,2)+sum(vox_2_mat,2)+sum(vox_3_mat,2);
    
    disp('sum_vox_mat');
    size(sum_vox_mat)
    
    vox_1_mat((sum_vox_mat>0)&(sum_vox_mat<3),:) = 0;
    vox_2_mat((sum_vox_mat>0)&(sum_vox_mat<3),:) = 0;
    vox_3_mat((sum_vox_mat>0)&(sum_vox_mat<3),:) = 0;
    toc
    
    disp('maxima vox matrices');
    max(max(vox_1_mat))
    max(max(vox_2_mat))
    max(max(vox_3_mat))
    
%     
%     figure()
%     hold on
%     plot(P1_YZ(2,:),P1_YZ(3,:),'.','Color','k')
%     plot(P_vox_YZ(2,:),P_vox_YZ(3,:),'.','Color','r')
%     hold off
%     axis equal
%     
%     figure()
%     hold on
%     plot(P1_YZ(2,:),P1_YZ(3,:),'.','Color','r')
%     plot(P1_grid(:,1),P1_grid(:,5),'.','Color','m')
%     plot(P1_grid(:,2),P1_grid(:,6),'.','Color','y')
%     plot(P1_grid(:,3),P1_grid(:,7),'.','Color','c')
%     plot(P1_grid(:,4),P1_grid(:,8),'.','Color','k')
%     hold off
%     axis equal
%     
%     figure()
%     hold on
%     plot(P2_XZ(1,:),P2_XZ(3,:),'.','Color','r')
%     plot(P2_grid(:,1),P2_grid(:,5),'.','Color','m')
%     plot(P2_grid(:,2),P2_grid(:,6),'.','Color','y')
%     plot(P2_grid(:,3),P2_grid(:,7),'.','Color','c')
%     plot(P2_grid(:,4),P2_grid(:,8),'.','Color','k')
%     hold off
%     axis equal
%     
%     figure()
%     hold on
%     plot(P3_XY(1,:),P3_XY(2,:),'.','Color','r')
%     plot(P3_grid(:,1),P3_grid(:,5),'.','Color','m')
%     plot(P3_grid(:,2),P3_grid(:,6),'.','Color','y')
%     plot(P3_grid(:,3),P3_grid(:,7),'.','Color','c')
%     plot(P3_grid(:,4),P3_grid(:,8),'.','Color','k')
%     hold off
%     axis equal
    
    
%     figure()
%     hold on
%     plot(P1_YZ(2,:),P1_YZ(3,:),'.','Color','r')
%     plot([P1_YZ(2,1) P1_YZ(2,256)],[P1_YZ(3,1) P1_YZ(3,256)],'Color','k')
%     plot([P1_YZ(2,256) P1_YZ(2,256^2)],[P1_YZ(3,256) P1_YZ(3,256^2)],'Color','k')
%     hold off
%     axis equal
%     
%     figure()
%     hold on
%     plot(P2_XZ(1,:),P2_XZ(3,:),'.','Color','g')
%     plot([P2_XZ(1,1) P2_XZ(1,256)],[P2_XZ(3,1) P2_XZ(3,256)],'Color','k')
%     plot([P2_XZ(1,256) P2_XZ(1,256^2)],[P2_XZ(3,256) P2_XZ(3,256^2)],'Color','k')
%     hold off
%     axis equal
%     
%     figure()
%     hold on
%     plot(P3_XY(1,:),P3_XY(2,:),'.','Color','b')
%     plot([P3_XY(1,1) P3_XY(1,256)],[P3_XY(2,1) P3_XY(2,256)],'Color','k')
%     plot([P3_XY(1,256) P3_XY(1,256^2)],[P3_XY(2,256) P3_XY(2,256^2)],'Color','k')
%     hold off
%     axis equal
    
%     tic
%     vox_YZ_ind = dsearchn([P1_YZ(2,:)' P1_YZ(3,:)'],D1_YZ,[P_vox_YZ(2,:)' P_vox_YZ(3,:)'],Inf);
%     vox_XZ_ind = dsearchn([P2_XZ(1,:)' P2_XZ(3,:)'],D2_XZ,[P_vox_XZ(1,:)' P_vox_XZ(3,:)'],Inf);
%     vox_XY_ind = dsearchn([P3_XY(1,:)' P3_XY(2,:)'],D3_XY,[P_vox_XY(1,:)' P_vox_XY(2,:)'],Inf);
%     toc
%     
%     size(vox_YZ_ind)
%     size(vox_XZ_ind)
%     size(vox_XY_ind)
%     
%     % Step 4: Remove all infinite points from the analysis.
%     
%     vox_3D_ind = 1:(257^3);
%     
%     vox_grid = [vox_3D_ind; vox_YZ_ind; vox_XZ_ind; vox_XY_ind; P_vox];
%     
%     vox_grid(sum(vox_grid,2)==Inf) = [];
%     
%     size(vox_grid)
%     
%     % Step 5: Plot voxel grid
%     
%     figure()
%     hold on
%     plot3(vox_grid(:,5),vox_grid(:,6),vox_grid(:,7),'.','Color','k')
%     hold off
%     axis equal
%     xlabel('X')
%     ylabel('Y')
%     zlabel('Z')

end

function R_w = return_rot_mat(calib_param)
%
    theta = sqrt(calib_param(4)^2+calib_param(5)^2+calib_param(6)^2);
    omega = [ 0 -calib_param(6) calib_param(5); ...
              calib_param(6) 0 -calib_param(4); ...
              -calib_param(5) calib_param(4) 0 ];
    R_w = eye(3)+(sin(theta)/theta)*omega+((1-cos(theta))/theta^2)*(omega*omega);

end

function X_w = camera_to_world_projection(calib_param,X_uv)
%uses dlt to find the coordinates of the image in 3D
%takes uv image coordinates and calculates how the points in 3d space
%project to camera image
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
    
    C_inv = pinv(C); %inverse for a nonsquare matrix
    
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