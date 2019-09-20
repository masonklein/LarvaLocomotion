%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%THOMAS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NOTE (07/09/2013): I don't think we want these commands here, it causes
%   problems when we have lots of data loaded because it clears all the 
%   variables.  The diary part might be useful but for now I'm removing
%   everything here.  (mjk)
%
% close all;
% clear;
% clc;
% diary diary;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%% BASIC LOCOMOTORY ANALYSIS VERSION 7 %%%%%%%%%%%%%%%%%%%
%
% Script loads data from a pair of .bin and .tim files; cleans tracks that
% are too short (in time or distance), too slow, or where the head/tail
% directions can't be determined very well; segments the tracks into runs
% and reorientations.  
%
% Then finds the mean larval speed, turning rate, head swing number, and 
% head swing angle.  Also prints the mean of these quantities for each
% track and generates a standard deviation and standard error accordingly.
%
% This latest version ... 
%
% [Updated 01/11/2012]
%
% [Updated again 03/10/2012 to add box-and-whisker numbers]
% 
% [Updated 01/27/2013 to include these changes:
%  * fixed standar d error for turn size, was using stdev/N instead of 
%    stdev/sqrt(N).
%  * fixed a mistake in the Excel file output where turn rate parameters
%    were mislabeled.
%  * added another measure for number of HS, the fraction of reorientations
%    with more than one head sweep.]
%
% [Updated 04/17/2013 to include these changes:
%  * added a feature that saves all the run durations (in seconds); 
%    the inverse of this number will be similar to the turn rate that
%    is already computed.]
%
% [Updated 06/30/2013 to blend with the adult fly analysis from the
%          Griffiths lab, to include these changes:
%  * now considering individual larvae instead of individual tracks, which
%    will be accomplished by only including tracks that are X minutes long
%    and that are active during some specific time Y 
%  * extract some "new" parameters, or at least labeled differently
%  * save information differently in a file, including the track point
%    locations
%  * get rid of a lot of the individual experiment saving stuff, just 
%    assume the strain and conditions are the same for each experiment
%  * remove the quartile/median/etc. material
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diary diary;
% mm / pixel conversion.  Use the values here or update depending on camera
% and plate position.  The easiest way to obtain this is to take a picture
% of a ruler sitting on the agar plate.  
lengthPerPixel = 1;       % (use if you want to use px/s)
%lengthPerPixel = 0.118;   % (with the 8X lens)
%lengthPerPixel = 0.0766; % (with the 12X lens)
lengthPerPixel = 0.0746;
%
% (NEW IN VERSION 7)
% ONE-TRACK-PER-ANIMAL RULES:
minTrackDuration = 240; % (in seconds)
maxTrackStartTime = 239; % (in seconds) 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% PART I: LOAD AND PREPARE DATA %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load experiment set; don't load if one already has been loaded
% NOTE: all of PART I will be skipped if it has already been run in the
%       current session of MATLAB
if (~exist('eset','var'))
    % Loads data from the .bin file, or perhaps multiple files, and stores
    % it in an ExperimentSet called "eset".  Asks the user for the location
    % of the .bin file(s).  Make sure the .bin files are all in the same
    % folder if you want to look at more than one.  
    eset = ExperimentSet.fromFiles();
    %
    % Defines an ESetCleaner called "ecl".  
    ecl = ESetCleaner;
    %
    %use this command to get report when analyzing data for the first time:
    %ecl.getReport(eset);
end

% Clean up "bad" tracks
existsAndDefault('cleanEset', 'true');
if (cleanEset)
    % Create an ESetCleaner called ecl.  (Redundant with the command above)
    ecl = ESetCleaner;
    
    % Properties of ecl.  These decide which tracks are thrown out.
    ecl.minHTValid = 0.95; % percentage of track where motion is in the same
                           % direction as a tail-to-head vector
    ecl.minDist = 85;      % distance (in pixels) that a track must span
    ecl.minSpeed = 0.50;   % average speed cutoff, in pixels/s (not mm/s)
    ecl.minPts = 800;      % minimum number of frames in the track 
                           % (e.g., 500 frames = 125 seconds at 4 Hz acq.)
    
    % This command performs the actual track cleaning:
    ecl.clean(eset);
    % 
    cleanEset = false;     % reset this to true if you want to re-clean    
end

% Fix head-tail orientation
existsAndDefault('fixht','true');
if (fixht)
    disp('fixing head-tail orientation');
    eset.executeTrackFunction('fixHTOrientation');
    fixht = false;
end

% Determine the speed thresholds used to decide when (1) a run ends and a
% new reorientation starts and (2) when a reorientation ends and a new run
% starts.  (These are not the same speeds).  These are found by looking at
% body orientation as a function of speed.  
existsAndDefault('autosetspeeds', true);
if (autosetspeeds)
    disp('setting segmentation speeds');
    eset.executeTrackFunction('setSegmentSpeeds');
    autosetspeeds = false;
end

% Segment the tracks (i.e., break them down into runs and reorientations)
existsAndDefault('segment', true);
if (segment)
    disp('segmenting tracks');
    eset.executeTrackFunction('segmentTrack');
    segment = false; 
end

% NEW PART: KEEP ONLY ONE TRACK PER ANIMAL
% (INTRODUCED IN VERSION 7)

for i=1:length(eset.expt)
    numTracks = length(eset.expt(i).track);
    %
    for j=1:numTracks
        startFrame = eset.expt(i).track(j).startFrame;
        numFrames = eset.expt(i).track(j).npts;
        startTime = eset.expt(i).track(j).dq.eti(1);
        duration = eset.expt(i).track(j).dq.eti(numFrames) - startTime;
        %
        if (startTime <= maxTrackStartTime && duration >= minTrackDuration)
            keepTrack(j) = true;
        else
            keepTrack(j) = false;
        end
        %
        if(isempty(eset.expt(i).track(j).run(1).runTime))
            keepTrack(j) = false;
        end
    end
    %
    % this command throws out the tracks we don't want:
    eset.expt(i).track = eset.expt(i).track(keepTrack);
    clear keepTrack
end
clear startFrame numFrames startTime duration;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% PART II: ANALYZE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% SIZE %%%%%
clear size;
k = 1;
for i=1:length(eset.expt)
   %
   for j=1:length(eset.expt(i).track)
      %
      size(k) = mean([eset.expt(i).track(j).pt.area]) * (lengthPerPixel^2);
      k=k+1;
      %
   end
   %
end

%%%%% SPEED %%%%%
clear speed speedmax ;
k=1;
for i=1:length(eset.expt)
   % 
   for j=1:length(eset.expt(i).track)
      speedVec = eset.expt(i).track(j).getDerivedQuantity('speed',false,'run');
      if(isempty(speedVec))
          speedVec = eset.expt(i).track(j).getDerivedQuantity('speed');
      end
      %
      speed(k) = mean(speedVec) * lengthPerPixel;
      speedmax(k) = max(speedVec) * lengthPerPixel;
      %
      k=k+1;
   end       
   % 
end
clear speedVec;

%%%%% ACCELERATION %%%%%
clear accel accelmax;
k=1;
for i=1:length(eset.expt)
   %
   for j=1:length(eset.expt(i).track)
      accelVec = eset.expt(i).track(j).getDerivedQuantity('acc',false,'run');
      if(isempty(accelVec))
          accelVec = eset.expt(i).track(j).getDerivedQuantity('acc');
      end
      accelVecX = accelVec(1,:);
      accelVecY = accelVec(2,:);
      accelVecMAG = sqrt( (accelVecX.*accelVecX)+(accelVecY.*accelVecY) );
      %
      accel(k) = mean(accelVecMAG) * lengthPerPixel;
      accelmax(k) = max(accelVecMAG) * lengthPerPixel;
      %
      k=k+1;
   end
   %
end
clear accelVec accelVecX accelVecY accelVecMAG;

%%%%% PATH LENGTH %%%%%
clear pathlength;
k=1;
for i=1:length(eset.expt)
   %
   for j=1:length(eset.expt(i).track)
      pathlengthVec = eset.expt(i).track(j).getDerivedQuantity('pathLength');
      %
      pathlength(k) = max(pathlengthVec) * lengthPerPixel;
      %
      k=k+1;
   end
   %
end
clear pathlengthVec;

%%%%% TRACK DURATION and RUN DURATION and RUN NUMBER %%%%%
clear trackduration runnumber runduration;
k=1;
for i=1:length(eset.expt)
   %
   for j=1:length(eset.expt(i).track)
      startFrame = eset.expt(i).track(j).startFrame;
      numFrames = eset.expt(i).track(j).npts;
      startTime = eset.expt(i).track(j).dq.eti(1);
      trackduration(k) = eset.expt(i).track(j).dq.eti(numFrames) - startTime;
      %
      runnumber(k) = length(eset.expt(i).track(j).run);
      %
      for m=1:length(eset.expt(i).track(j).run)
         runTimeVec(m) = eset.expt(i).track(j).run(m).runTime; 
      end
      runduration(k) = mean(runTimeVec);
      %
      k=k+1;
      clear runTimeVec;       
   end
   %
end
clear startFrame numFrames startTime;

%%%%% NUMBER OF HEAD SWEEPS %%%%%
clear HSnum HSnumfrac;
k=1;
for i=1:length(eset.expt)
   %
   for j=1:length(eset.expt(i).track)
      reoVec = eset.expt(i).track(j).reorientation;
      numHSvec = [reoVec.numHS];
      numHSvecgtrzero = numHSvec(numHSvec>0);
      numHSvecgtrone = numHSvec(numHSvec>1);
      %
      HSnum(k) = mean(numHSvec(numHSvec>0));
      HSnumfrac(k) = length(numHSvecgtrone)/length(numHSvecgtrzero);
      %
      k=k+1;
   end
end
clear reoVec numHSvec numHSvecgtrzero numHSvecgtrone;

%%%%% TURN SIZE %%%%%
clear HSsize;
k=1;
for i=1:length(eset.expt)
   %
   for j=1:length(eset.expt(i).track)
      reoVec = eset.expt(i).track(j).reorientation;
      %
      if(isempty(reoVec))
          HSsize(k) = NaN;
      else 
          numHSvec = [reoVec.numHS];
          reoVec2 = reoVec(numHSvec>0);
          %
          % lastHS = repmat(HeadSwing, size(reoVec2));
          for m=1:length(reoVec2)
            lastHS(m) = reoVec2(m).headSwing(end); 
          end
          lastHSmaxtheta = [lastHS.maxTheta];
          lastHSmaxtheta = abs(rad2deg(lastHSmaxtheta));
          %
          HSsize(k) = mean(lastHSmaxtheta);
      end
      %
      k=k+1;
   end
end
clear reoVec numHSvec reoVec2 lastHS lastHSmaxtheta;

%%%%% AVERAGES + UNCERTAINTIES %%%%%
numTracks = 0;
for i=1:length(eset.expt)
   numTracks = numTracks + length(eset.expt(i).track); 
end
%
trackduration_m = mean(trackduration);
trackduration_sd = std(trackduration);
trackduration_sem = trackduration_sd/sqrt(numTracks-1);
%
size_m = mean(size);
size_sd = std(size);
size_sem = size_sd/sqrt(numTracks-1);
%
speed_m = mean(speed);
speed_sd = std(speed);
speed_sem = speed_sd/sqrt(numTracks-1);
%
speedmax_m = mean(speedmax);
speedmax_sd = std(speedmax);
speedmax_sem = speedmax_sd/sqrt(numTracks-1);
%
accel_m = mean(accel);
accel_sd = std(accel);
accel_sem = accel_sd/sqrt(numTracks-1);
%
accelmax_m = mean(accelmax);
accelmax_sd = std(accelmax);
accelmax_sem = accelmax_sd/sqrt(numTracks-1);
%
pathlength_m = mean(pathlength);
pathlength_sd = std(pathlength);
pathlength_sem = pathlength_sd/sqrt(numTracks-1);
%
runnumber_m = mean(runnumber);
runnumber_sd = std(runnumber);
runnumber_sem = runnumber_sd/sqrt(numTracks-1);
%
runduration_m = mean(runduration);
runduration_sd = std(runduration);
runduration_sem = runduration_sd/sqrt(numTracks-1);
%
HSnum_m = nanmean(HSnum);
HSnum_sd = nanstd(HSnum);
HSnum_sem = HSnum_sd/sqrt(numTracks-1);
HSnumfrac_m = nanmean(HSnumfrac);
HSnumfrac_sd = nanstd(HSnumfrac);
HSnumfrac_sem = HSnumfrac_sd/sqrt(numTracks-1);
%
HSsize_m = nanmean(HSsize);
HSsize_sd = nanstd(HSsize);
HSsize_sem = HSsize_sd/sqrt(numTracks-1);
%
means = [trackduration_m size_m speed_m speedmax_m accel_m accelmax_m];
means = [means pathlength_m runnumber_m runduration_m HSnum_m HSnumfrac_m HSsize_m];
stds = [trackduration_sd size_sd speed_sd speedmax_sd accel_sd accelmax_sd];
stds = [stds pathlength_sd runnumber_sd runduration_sd HSnum_sd HSnumfrac_sd HSsize_sd];
sems = [trackduration_sem size_sem speed_sem speedmax_sem accel_sem accelmax_sem];
sems = [sems pathlength_sem runnumber_sem runduration_sem HSnum_sem HSnumfrac_sem HSsize_sem];
averages = [means;stds;sems];

% Convert the head sweep stuff to cells as a way to exclude NaNs and have
% them show up as blank cells in Excel (instead of 65535 because Matlab
% is stupid)
HSnum = num2cell(HSnum);
HSnumfrac = num2cell(HSnumfrac);
HSsize = num2cell(HSsize);
HSnum(cellfun(@isnan,HSnum))={[]};
HSnumfrac(cellfun(@isnan,HSnumfrac))={[]};
HSsize(cellfun(@isnan,HSsize))={[]};
%
trackduration = num2cell(trackduration);
size = num2cell(size);
speed = num2cell(speed);
speedmax = num2cell(speedmax);
accel = num2cell(accel);
accelmax = num2cell(accelmax);
pathlength = num2cell(pathlength);
runnumber = num2cell(runnumber);
runduration = num2cell(runduration);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% PART III: DISPLAY RESULTS %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

display('***************************************************************');

display(['Analyzed ' num2str(length(eset.expt)) ' experiments']);
for j=1:length(eset.expt)
    display(['  experiment #' num2str(j) ': ' eset.expt(j).fname]);
end

if(~isempty(eset.expt))
    
    display('***************************************************************');
    
    display('EXPERIMENTS COMBINED'  );
    display(['  ' num2str(numTracks) ' total larvae']);
    
    display('TRACK FILTERING');
    display(['  minHTValid = ' num2str(ecl.minHTValid,3)]);
    display(['  minDist = ' num2str(ecl.minDist,3) ' pixels']);
    display(['  minSpeed = ' num2str(ecl.minHTValid,3) ' pixels/s']);
    display(['  minPts = ' num2str(ecl.minPts,3) ' points']);
    display(['  [only tracks > ' num2str(minTrackDuration) 's, starting < ' num2str(maxTrackStartTime) ' into experiment]']);
    
    display('SIZE');
    display(['  ' num2str(size_m) ' mm^2']);
    display(['  (SD=' num2str(size_sd,3) ', SEM=' num2str(size_sem,3) ')']);
    
    display('SPEED');
    display(['  ' num2str(speed_m,3) ' mm/s']);
    display(['  (SD=' num2str(speed_sd,3) ', SEM=' num2str(speed_sem,3) ')']);
    
    display('ACCELERATION');
    display(['  ' num2str(accel_m,3) ' mm/s/s']);
    display(['  (SD=' num2str(accel_sd,3) ', SEM=' num2str(accel_sem,3) ')']);
    
    display('RUN DURATION');
    display(['  ' num2str(runduration_m,3) ' s']);
    display(['  (SD=' num2str(runduration_sd,3) ', SEM=' num2str(runduration_sem,3) ')']);
    
    display('HEAD SWEEP NUMBER');
    display(['  ' num2str(HSnum_m,3) ' sweeps/reorientation']);
    display(['  (SD=' num2str(HSnum_sd,3) ', SEM=' num2str(HSnum_sem,3) ')']);
    display(['  ' num2str(HSnumfrac_m,3) ' fraction of reorientations with > 1 sweep']);
    display(['  (SD=' num2str(HSnumfrac_sd,3) ', SE=' num2str(HSnumfrac_sem,3) ')']);
    
    display('HEAD SWEEP DEPTH');
    display(['  ' num2str(HSsize_m,3) ' deg (last sweep only)']);
    display(['  (SD=' num2str(HSsize_sd,3) ', SEM=' num2str(HSsize_sem,3) ')']);
       
end

display('***************************************************************');

diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART IV: SAVE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n\n - DATA SAVING...\n')
saveFileName = input('\nEnter file name (without `.` or `\\`) : ', 's');

folder=0;
while (folder ~= 1) && (folder ~= 2)
    fprintf('\nWhere do you want to save the results? \n')
    disp(['1 : current folder : ' pwd])
    disp('2 : choose a folder')
    folder=input('Type 1 or 2 : ');
end
 
if folder==2
    folder_name = uigetdir;
    saveFileName = [folder_name '\' saveFileName];
end

saveFileName1 = [saveFileName '_MAT'];
saveFileName2 = [saveFileName '_RESULTS'];
saveFileName3 = [saveFileName '_TEXT'];
saveFileName4 = [saveFileName '_TRACKS'];

% MATLAB FILE:
% Save these [twelve] variables to a file 
save(saveFileName1,'size','speed','speedmax','accel','accelmax','pathlength','trackduration','runnumber','runduration','HSnum','HSnumfrac','HSsize');
% The file can be renamed after it's saved.  To open the
% variables within the .m file, just drag the file into the workspace.
% They can be renamed there too, which would probably be a good idea for
% comparing these quantities across different experiments.  

% EXCEL FILE:
warning('off', 'MATLAB:xlswrite:AddSheet');
% Column Labels:
labels = {'larva','duration','size','<speed>','speedmax','<|accel|>','|accelmax|'};
labels = [labels {'pathlength','# runs','<run duration>'}];
labels = [labels {'<#HS/reo>','<frac #HS>1>','turnsize'}];
% Larva Labels:
larvalabels = '';
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        string = ['e' num2str(i) 't' num2str(j)];
        larvalabels = [larvalabels {string}];
    end
end
larvalabels = transpose(larvalabels);
larvalabelsrange = ['A2:A' num2str(numTracks+1)];
clear string;
% Data:
data = [trackduration;size;speed;speedmax;accel;accelmax;pathlength];
data = [data;runnumber;runduration;HSnum;HSnumfrac;HSsize];
data = transpose(data);
datarange = ['B2:M' num2str(numTracks+1)];
% Averages:
avglabels = {'mean';'sd';'sem'};
avglabelsrange = ['A' num2str(numTracks+3) ':A' num2str(numTracks+5)];
avgrange = ['B' num2str(numTracks+3) ':M' num2str(numTracks+5)];
% Write to File:
xlswrite(saveFileName2,labels);
xlswrite(saveFileName2,larvalabels,larvalabelsrange);
xlswrite(saveFileName2,data,datarange);
xlswrite(saveFileName2,avglabels,avglabelsrange);
xlswrite(saveFileName2,averages,avgrange);

% DIARY FILE:
movefile('diary',[saveFileName3 '.doc'],'f');

% TRACKS:
k = 1;
for i=1:length(eset.expt)
    for j=1:length(eset.expt(i).track)
        % Labels
        sheetName = ['Track' num2str(k) '(e' num2str(i) 't' num2str(j) ')'];
        labels = {'time','x','y','run?'};
        % Times:
        times = eset.expt(i).track(j).dq.eti;
        numPoints = length(times);
        % Positions:
        pos = eset.expt(i).track(j).getDerivedQuantity('sloc');
        xpos = pos(1,:);
        xpos = xpos*lengthPerPixel;
        ypos = pos(2,:);
        ypos = ypos*lengthPerPixel;
        % Run Y/N:
        runYN = eset.expt(i).track(j).isrun;
        % Combined Matrix:
        combined = [times;xpos;ypos;runYN];
        combined = transpose(combined);
        combinedrange = ['A2:D' num2str(numPoints+1)];
        % Write to Excel:
        xlswrite(saveFileName4,labels,sheetName);
        xlswrite(saveFileName4,combined,sheetName,combinedrange);
        %
        k=k+1;
    end
end
clear times pos xpos ypos runYN combined combinedrange;

fprintf('\ndone\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






