%% 1. Load 179 EphysUnits, and select 'good' 4s+ reach trials.
eu = EphysUnit.load('C:\SERVER\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate');
for iEu = 1:length(eu)
    trials = eu(iEu).makeTrials('press_spontaneous2');
    goodTrials = trials(trials.duration() > 4);
    eu(iEu).Trials.PressSpontaneous = goodTrials;
end

% Group EphysUnits by session (euByExp)
pa = Pawnalyzer2(eu, refEvent='reward');
euByExp = arrayfun(@(x) x.eu, pa.exp, UniformOutput=false);
clear pa iEu trials goodTrials

%% 2. Label paw tracjectories, one session at a time
% There's not enough RAM to load video clips from all 10 sessions, so we
% will label one session, save, and then load the next. Make sure to
% increment the value of iExp when you're done with a session

% 2.1 Read video clips from ONE session.
iExp = 1; % !!!Change this to 2, 3, 4, etc to work on other sessions.
clear pa, close all
pa = Pawnalyzer2(euByExp{iExp}, refEvent='press');
pa.getClips(nFramesBefore=15, nFramesAfter=0, keepData=false, trials='PressSpontaneous');

% 2.2 This opens the user interface. If you close it by mistake run this to 
% restart. Your labelled data should not be affected unless you exit MATLAB
pa.start();

% Some keybinds to remember
%   Left click for Contra (red circle), right click for Ipsi (yellow circle)
%   Click and drag to move the circle
%   Alt + Click and drag to precisely move the circle
%   Capslock + Click to turn on tracking mode
%   Mousewheel or arrow keys to change frame.
%   Shift + Mousewheel or Shift + arrow key to change trial.
%   Z or C to go to the first or last frame in the trial.
%   Ctrl+S to save, see step 2.3 below

% 2.3 Press "Ctrl+S" or type "pa.save()" in the command window to save your 
% progress. Unsaved progress will be lost when MATLAB is closed so SAVE OFTEN!
% The first time you save (for each session) the program will ask where you
% want to save the file, the default path should be fine, just make sure
% each session is saved to its own file so you don't overwrite your progress.

% 2.4 If MATLAB is restarted, rerun the script from the top until you get
% here, then type "pa.load()" in the command window, select your previously saved data. 

%% 3. You have finished labelling one sessions, now return to Step 2.1
% Make sure you change the value of 'iExp' on line 20. This determines
% which session is loaded.