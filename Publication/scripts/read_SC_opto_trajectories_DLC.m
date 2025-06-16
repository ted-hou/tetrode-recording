
euSC = EphysUnit.load('C:\SERVER\Units\TwoColor_SC\SingleUnit_NonDuplicate_NonDrift_SC');
%
expSC = CompleteExperiment3(euSC, cameras='r');
expSC.alignTimestamps(refEventNameArduino={'LEVER_PRESSED', 'OPTO1_ON', 'OPTO2_ON'}, refEventNameEphys={'PressOn', 'StimOn'}, trialDurationTolerance=0.6);
% exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON'}, refEventNameEphys={'StimOn'}, trialDurationTolerance=0.15);

%%
pSC.stimBluePowers = [500]*1e-6; 
pSC.stimRedPowers = [16000]*1e-6;
pSC.stimBlueDurations = [20]*1e-3;
pSC.stimRedDurations = [20]*1e-3;

clear groupsSC trajectoriesSC
groupsSC(length(expSC)) = struct(blue=[], red=[]);
trajectoriesSC(length(expSC)) = struct(HandCameraSide=[], HandOppositeSide=[], Jaw=[]);
for iExp = 1:length(expSC)
    groupsSC(iExp).blue = expSC(iExp).eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=pSC.stimBluePowers, duration=pSC.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsSC(iExp).red = expSC(iExp).eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=pSC.stimRedPowers, duration=pSC.stimRedDurations, location=[], wavelength=635));
    
    for bodypart = ["HandCameraSide", "HandOppositeSide", "Jaw"]
        for color = ["blue", "red"]
            [trajectoriesSC(iExp).(bodypart).(color).X, trajectoriesSC(iExp).(bodypart).(color).Y, trajectoriesSC(iExp).(bodypart).(color).L, trajectoriesSC(iExp).(bodypart).(color).t] = expSC(iExp).getTrajectoryByTrial('r', char(bodypart), trials=[groupsSC(iExp).(color).trials], window=[-0.5, 0.5], likelihoodThreshold=0.9, includeInvalid=true, alignTo='start', interp='previous');
        end
    end
end

clear iExp bodypart color

% Make session average
for bodypart = ["HandCameraSide", "HandOppositeSide", "Jaw"]
    for color = ["blue", "red"]
        for field = ["X", "Y", "L"]
            fielddata = arrayfun(@(traj) traj.(bodypart).(color).(field), trajectoriesSC(1:length(expSC)), UniformOutput=false);
            trajectoriesSC(length(expSC)+1).(bodypart).(color).(field) = cat(1, fielddata{:});
        end
        trajectoriesSC(length(expSC)+1).(bodypart).(color).t = trajectoriesSC(1).(bodypart).(color).t;
    end
end
clear bodypart color field fielddata

% %% Get video clips around time of stim (for manual validation of video/ephys alignment)
% clear clips
% for color = ["blue", "red"]
%     switch color
%         case "blue"
%             trials = [groupsBlue.trials];
%         case "red"
%             trials = [groupsRed.trials];
%     end
%     [clips.(color), t] = expSC.getVideoClip([trials.Start], 'r', numFramesBefore=15, numFramesAfter=15, bodyParts={'HandCameraSide', 'Jaw'}, file='C:\SERVER\desmond39\desmond39_20250522\desmond39_20250522_laser_1.mp4');
% end