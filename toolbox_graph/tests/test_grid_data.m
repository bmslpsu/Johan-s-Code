clear all;
close all;
clc;

addpath('C:/Users/jmelis/Desktop/new_flight_tracker/toolbox_graph/toolbox_graph');
addpath('C:/Users/jmelis/Desktop/new_flight_tracker/toolbox_graph/toolbox_graph/toolbox');
addpath('C:/Users/jmelis/Desktop/new_flight_tracker/toolbox_graph/toolbox_graph/off');

% test for data gridding
p = 7;
n = 128;
vertex = rand(2,p)*(n-1)+1;
face = delaunay(vertex(1,:), vertex(2,:));
v = vertex(1,:) .* vertex(2,:);     % X*Y function

% using the toolbox gridding (works for arbitray set of points)
M1 = griddata_arbitrary(face,vertex,v,n);

% using matlab gridding (works for grid sampling)
[Y,X] = meshgrid(1:n,1:n);
M2 = griddata(vertex(1,:), vertex(2,:), v, X,Y);