% Replace 'your_directory_path' with the path to the directory you want to add
directory_path = 'C:\Siddharth\June-July_2023_Internship\MSS';

% Use the 'genpath' function to generate a list of 'directory_path' and all its subdirectories
all_subdirectories = genpath(directory_path);

% Add all the subdirectories to the MATLAB search path
addpath(all_subdirectories);

% To save this path for future MATLAB sessions, use 'savepath'
savepath;
