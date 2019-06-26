clc
clear all
root='S:\Public\Wael\Cal4\Test vid\TestVid1\';
calibration_dir='S:\Public\Wael\Cal4\DLT\';
%% reads an image
v_1 = VideoReader([root 'fly_calibrated_640x480_C001H001S0001.avi']);
v_2 = VideoReader([root 'fly_calibrated_640x480_C002H001S0001.avi']);
v_3 = VideoReader([root 'fly_calibrated_640x480_C003H001S0001.avi']);

img_1_t = readFrame(v_1);
img_2_t = readFrame(v_2);
img_3_t = readFrame(v_3);
%% reads DLTcoeff
load([calibration_dir 'DLT.txt']);%loads the DLT parameters (WS)
%% plots to check 

imshow(img_1_t)
[u,v] = getpts
%% do partial DLT
[m,b]=partialdlt(u,v,DLT(:,1),DLT(:,2))
t=1:1:400;
figure
imshow(img_2_t)
hold on
plot(t,t*m+b)
%% functions

function [m,b]=partialdlt(u,v,C1,C2)


% partialdlt takes as inputs a set of X,Y coordinates from one camera view
% and the DLT coefficients for the view and one additional view.  It
% returns the line coefficients m & b for Y=mX+b in the 2nd camera view.
% Under error-free DLT, the X,Y marker in the 2nd camera view must fall
% along the line given.
%
% Inputs:
%	u = X coordinate in camera 1
%	v = Y coordinate in camera 1
%	C1 = the 11 dlt coefficients for camera 1
%	C2 = the 11 dlt coefficients for camera 2
%
% Outputs:
%	m = slope of the line in camera 2
%	b = Y-intercept of the line in camera 2
%
% Initial version:
% Ty Hedrick, 11/25/03
%	6/14/04 - Ty Hedrick, updated to use the proper solution for x(i)

% pick 2 random Z (actual values are not important)
z(1)=500;
z(2)=-500;

% for each z predict x & y
x=z*NaN;
y=z*NaN;
for i=1:2
  Z=z(i);
  
  y(i)= -(u*C1(9)*C1(7)*Z + u*C1(9)*C1(8) - u*C1(11)*Z*C1(5) -u*C1(5) + ...
    C1(1)*v*C1(11)*Z + C1(1)*v - C1(1)*C1(7)*Z - C1(1)*C1(8) - ...
    C1(3)*Z*v*C1(9) + C1(3)*Z*C1(5) - C1(4)*v*C1(9) + C1(4)*C1(5)) / ...
    (u*C1(9)*C1(6) - u*C1(10)*C1(5) + C1(1)*v*C1(10) - C1(1)*C1(6) - ...
    C1(2)*v*C1(9) + C1(2)*C1(5));
  
  Y=y(i);
  
  x(i)= -(v*C1(10)*Y+v*C1(11)*Z+v-C1(6)*Y-C1(7)*Z-C1(8))/(v*C1(9)-C1(5));
end

% back the points into the cam2 X,Y domain
xy(1:2,1:2)=NaN;
for i=1:2
  xy(i,:)=dlt_inverse(C2(:),[x(i),y(i),z(i)]);
end

% get a line equation back, y=mx+b
m=(xy(2,2)-xy(1,2))/(xy(2,1)-xy(1,1));
b=xy(1,2)-m*xy(1,1);
end

function [uv] = dlt_inverse(c,xyz)

% function [uv] = dlt_inverse(c,xyz)
%
% This function reconstructs the pixel coordinates of a 3D coordinate as
% seen by the camera specificed by DLT coefficients c
%
% Inputs:
%  c - 11 DLT coefficients for the camera, [11,1] array
%  xyz - [x,y,z] coordinates over f frames,[f,3] array
%
% Outputs:
%  uv - pixel coordinates in each frame, [f,2] array
%
% Ty Hedrick

% write the matrix solution out longhand for Matlab vector operation over
% all points at once
uv(:,1)=(xyz(:,1).*c(1)+xyz(:,2).*c(2)+xyz(:,3).*c(3)+c(4))./ ...
  (xyz(:,1).*c(9)+xyz(:,2).*c(10)+xyz(:,3).*c(11)+1);
uv(:,2)=(xyz(:,1).*c(5)+xyz(:,2).*c(6)+xyz(:,3).*c(7)+c(8))./ ...
  (xyz(:,1).*c(9)+xyz(:,2).*c(10)+xyz(:,3).*c(11)+1);
end