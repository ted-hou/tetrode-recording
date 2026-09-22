script_loadPreprocessedData_optrodeSNr_CoChR_VGATCRE_20260729;

%%
for iExp = 1:length(exp)
    eu0 = eu(expToEuIndices(iExp));
    % Find all non-stim reach/lick trials
    moveTimes = struct(press=eu0.EventTimes.FirstPress, lick=eu0.EventTimes.FirstLick);
    stimOn = eu0.EventTimes.LaserModBlueOn;
    stimOff = eu0.EventTimes.LaserModBlueOff;
    for trialType = ["press", "lick"]
        isCleanTrial = true(length(moveTimes.(trialType)), 1); % Not-stim adjacent
        for iTrial = 1:length(moveTimes.(trialType))
            tMove = moveTimes.(trialType)(iTrial);
            if any(isin(stimOn, tMove + [-5, 5]))
                isCleanTrial(iTrial) = false;
            end
        end
        moveTimes.(trialType) = moveTimes.(trialType)(isCleanTrial);

        % For each trial: get:
        % 1. movement onset time
        % 2. spike rate timecourse
        % 3. LFP timecourse (need to re-read raw data)
        for iTrial = 1:length(moveTimes.(trialType))

        end
    end
end