%this uses 3D Orthogonal Distance Regression (ODR)
%at the moment I'm looking up how to use this fitting method. I found a
%sample code online and it seems to work well. Future research is required
%to verify how this method works

%from sushiel: check introduction to POD by anindy chatterjeedya for and
%introduction to svd
% Step 1: arrange points in a matrix format
points = [r c v];

% Step 2: find the mean of the points
avg = mean(points, 1);

% Step 3: subtract the mean from all points
subtracted = bsxfun(@minus, points, avg);

% Step 4 : perform SVD
[~, ~, V] = svd(subtracted);

% Step 5: find the direction vector 
%        (which is the right singular vector corresponding to the largest singular value)
direction = V(:, 1);

% Line is 'avg' and 'direction'
p0 = avg;
d = direction;

% Parametric equation: P = p0 + t*d
t=-100:1:170;
x_svd=avg(1)+direction(1)*t;
y_svd=avg(2)+direction(2)*t;
z_svd=avg(3)+direction(3)*t;
figure 
scatter3(x_svd,y_svd,z_svd)
hold on
scatter3(X(:), Y(:), Z(:), pointsize, im_33d(idx));
colormap(gray(256));