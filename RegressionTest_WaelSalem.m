%finds the best regression curve for a line in 2D
%issues: might be computationally expensive due to finding the inverse of a
%big matrix
x=1:1:20;
y=x+3*randn(length(x),1)'
plot(x,y)

A=[ones(length(x),1) x']
als=inv(transpose(A)*A)*transpose(A)*y'
%%
y1=als(1)+x*als(2)
hold on 
plot(x,y1)
%% 3D regression test of some data
im_3dtest=im_33d;
im_3dtest(im_3dtest<0.4)=0; 
im_3dtest(im_3dtest>=0.4)=nan;
idx = find(im_3dtest);
[X, Y, Z] = ind2sub(size(im_3dtest), idx);
pointsize = 20;
figure
scatter3(X(:), Y(:), Z(:), pointsize, im_3dtest(idx));
[x, y, z] = find(isnan(im_3dtest)); %point of the fly that I'm interested in are converted to nan and I find their indicies
hold on
scatter3(x,y,z)
[r,c,v] = ind2sub(size(im_3dtest),find(isnan(im_3dtest)));
scatter3(r,c,v) %this seems to give me the xyz coordinates, now time to fit a line
A_3D=[ones(length(r),1) r c ];
als_3D=inv(transpose(A_3D)*A_3D)*transpose(A_3D)*v;
%% works, i can fit a line to the fly, now time to look at alignments
%the only issue here may be speed of the function.
x_3d=linspace(min(r),max(r),100);
y_3d=linspace(min(c),max(c),100);
z_3d=als_3D(1)+x_3d*als_3D(2)+y_3d*als_3D(3);
figure
scatter3(X(:), Y(:), Z(:), pointsize, im_33d(idx));
colormap(gray(256));
hold on
scatter3(x_3d,y_3d,z_3d)
%%
z_3d_2=als_3D(1)+r*als_3D(2)+c*als_3D(3);

scatter3(r,c,z_3d_2)



