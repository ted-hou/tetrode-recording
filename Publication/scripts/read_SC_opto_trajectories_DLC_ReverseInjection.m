
euSCRevInj = EphysUnit.load('C:\SERVER\Units\TwoColor_SC\ReverseInjection_Dummy');
%
expSCRevInj = CompleteExperiment3(euSCRevInj, cameras='r');
expSCRevInj.alignTimestamps(refEventNameArduino={'LEVER_PRESSED', 'OPTO1_ON', 'OPTO2_ON'}, refEventNameEphys={'PressOn', 'StimOn'}, trialDurationTolerance=0.8);
% exp.alignTimestamps(refEventNameArduino={'OPTO1_ON', 'OPTO2_ON'}, refEventNameEphys={'StimOn'}, trialDurationTolerance=0.15);

expSCRevInj = expSCRevInj(2:end);

%%
pSCRevInj.stimBluePowers = [16000]*1e-6; 
pSCRevInj.stimRedPowers = [16000]*1e-6;
pSCRevInj.stimBlueDurations = [250]*1e-3;
pSCRevInj.stimRedDurations = [250]*1e-3;

clear groupsSCRevInj trajectoriesSCRevInj
groupsSCRevInj(length(expSCRevInj)) = struct(blue=[], red=[]);
trajectoriesSCRevInj(length(expSCRevInj)) = struct(HandCameraSide=[], HandOppositeSide=[], Jaw=[]);
for iExp = 1:length(expSCRevInj)
    groupsSCRevInj(iExp).blue = expSCRevInj(iExp).eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=pSCRevInj.stimBluePowers, duration=pSCRevInj.stimBlueDurations, location=[], wavelength=[470, 473]));
    groupsSCRevInj(iExp).red = expSCRevInj(iExp).eu(1).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, selectBy=struct(power=pSCRevInj.stimRedPowers, duration=pSCRevInj.stimRedDurations, location=[], wavelength=635));
    
    for bodypart = ["HandIpsi", "HandCont", "Jaw"]
        for color = ["blue", "red"]
            [trajectoriesSCRevInj(iExp).(bodypart).(color).X, trajectoriesSCRevInj(iExp).(bodypart).(color).Y, trajectoriesSCRevInj(iExp).(bodypart).(color).L, trajectoriesSCRevInj(iExp).(bodypart).(color).t] = expSCRevInj(iExp).getTrajectoryByTrial('r', char(bodypart), trials=[groupsSCRevInj(iExp).(color).trials], window=[-0.5, 1], likelihoodThreshold=0.9, includeInvalid=true, alignTo='start', interp='previous');
        end
    end
end

clear iExp bodypart color

% Make session average
for bodypart = ["HandIpsi", "HandCont", "Jaw"]
    for color = ["blue", "red"]
        for field = ["X", "Y", "L"]
            fielddata = arrayfun(@(traj) traj.(bodypart).(color).(field), trajectoriesSCRevInj(1:length(expSCRevInj)), UniformOutput=false);
            trajectoriesSCRevInj(length(expSCRevInj)+1).(bodypart).(color).(field) = cat(1, fielddata{:});
        end
        trajectoriesSCRevInj(length(expSCRevInj)+1).(bodypart).(color).t = trajectoriesSCRevInj(1).(bodypart).(color).t;
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


