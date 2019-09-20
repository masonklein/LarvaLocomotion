% Runs MinimalAnalysis7multi.m script a bunch of times
% Structure:
% 
% Ask user for how many FOLDERS (strains?) to look at
% 
% Get folder paths from user, and the save file name for each
%
% 
fprintf('---Minimal Analysis v7 batch program!---\n');


% User chooses the number of folders.  Matlab will analyze *all* the
% .bin files found in each folder.  So for doing strains where there were
% multiple experiments, put their .bin files in the same folder.  
clear global;
global NUMFOLDERS;
global DIRECTORYLIST;
global CURRENTDIRECTORY;
global SAVEFILENAMELIST;
global CURRENTSAVEFILENAME;
global CURRENTFOLDER;

NUMFOLDERS = input('How many folders do you want to analyze? ');

for i=1:NUMFOLDERS
    DIRECTORYLIST{i} = uigetdir;
    SAVEFILENAMELIST{i} = input('save file name for that folder: ','s');    
end

CURRENTFOLDER = 1;
while CURRENTFOLDER <= NUMFOLDERS
    CURRENTDIRECTORY = DIRECTORYLIST(CURRENTFOLDER);
    CURRENTSAVEFILENAME = SAVEFILENAMELIST(CURRENTFOLDER);
    
    MinimalAnalysis7multi; 
    
    clear;
    close all;
    global NUMFOLDERS;
    global DIRECTORYLIST;
    global CURRENTDIRECTORY;
    global SAVEFILENAMELIST;
    global CURRENTSAVEFILENAME;
    global CURRENTFOLDER;   
    CURRENTFOLDER = CURRENTFOLDER + 1;   
end    
