clear all;
close all;
clc;

% mov_folder = 'G:\Flyami_movies\Session_12_Dec_2017_10_35';
%mov_folder = 'G:\Flyami_movies\Session_11_Dec_2017_11_24';
mov_folder = 'G:\Flyami_movies\Session_11_Jan_2018_09_46';

for i = [1]
    
    mov_nr = i;
    
    tethered_flight_tracker(mov_folder,mov_nr);
    
end