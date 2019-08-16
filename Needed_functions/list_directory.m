function [ dirlist ] = list_directory( dir_path, file_ext )
% LIST specified folder and specified file name extension.
% It can be UNIX or Windows path to specify. Without specifying
% an second parameter it uses all files in that folder.
%
% Michael Tesar,
% 2016, Ceske Budejovice
% Version 1.0
%
%    INPUT:
%       dir_path - a string containing path to folder containing files
%       file_ext - a string which filter files
%
%    OUTPUT:
%       dirlist - string or char array of files meet criteria
%
% Example: [all_mp3] = list_directory ('/Users/usr/music/','.mp3');
%          [files] = list_directory ('/Users/usr/Documents/');
%
% Access data: listDir(char(index)) return string of indexed file meets
% specified criteria.
% Hint: It also can be used not for file extension but also for substring
% filtering of files.
%% Directory listing
% ==================
% Declare global manual iterator
ii = 1;
% Check if both parameters are string variable
if ischar (dir_path) && ischar(file_ext)
     listDir  = dir(dir_path);          % List directory
     listName = {listDir.name}';        % Extract only names
     allFiles = char(listName(3:end));  % Delete navigation strings . and ..
     N = length(listName(3:end));       % Get the number of total files
     % Loop for all files
     for i = 1:N
         % Check if meet file extension criteria
         if strfind(allFiles(i,:), file_ext);
             % If do then add to another variable and ++ manual iterator
             dirlist(ii,:) = {char(allFiles(i,:))};
             ii = ii + 1;
         else
             continue
         end
     end
else
    % Both parameters needs to be strings
    error('Enter strings!')
end