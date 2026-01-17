
%% Spike detection
folders = { ...
    'C:\SERVER\daisy27\daisy27_20250707', ...
    'C:\SERVER\daisy27\daisy27_20250715', ... % 16mW mirror 0, 100ms moved arm
    'C:\SERVER\daisy28\daisy28_20250714', ... % 16mW mirror 0, 100ms moved arm
    };

for iSession = 1:length(folders)
    try
        tr = TetrodeRecording;
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.ReadFiles(Duration=240, Channels=1:384, NumSigmas=4, NumSigmasReturn=1.5, NumSigmasReject=20, WaveformWindow=[-0.5, 1])
        tr.SaveNeuropixelIO()
        
        % Read detected spikes and NIDQ digital/analog channels
        % tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO();
        
        % Spike sort
        tr.LoadSpikes(1);
        tr.IterativeArtifactRemoval(1, MinSpikeRate=0.5, KIterative=4, KFinal=2, MaxIters=5, ...
            DimensionIterative=3, DimensionFinal=10, FeatureMethod='PCA', ClusterMethod='kmeans', ...
            WaveformWindow=[-0.5, 0.5]);
        tr.SaveSpikes(Path='Spikes_Dummy');
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end


% Convert to EphysUnits
folders = { ...
    'C:\SERVER\daisy27\daisy27_20250707', ...
    'C:\SERVER\daisy27\daisy27_20250715', ... % 16mW mirror 0, 100ms moved arm
    'C:\SERVER\daisy28\daisy28_20250714', ... % 16mW mirror 0, 100ms moved arm
    };

chunkSize = 1; % NumChannelsPerChunk
for iSession = 1:length(folders)
    try
        tr = TetrodeRecording();
        tr.SelectFiles(NeuropixelPath=folders{iSession});
        tr.LoadNeuropixelIO();
        tr.ParseNeuropixelIO();
        
        % IterativeArtifactRemoval (OnionPeeling): Load spikes
        for iChunk = 1:(384/chunkSize)
            channels = (iChunk-1)*chunkSize + 1 : iChunk*chunkSize;
            tr.LoadSpikes(channels, Path='Spikes_Dummy');

            if isempty(tr.Spikes) || isempty([tr.Spikes.Channel])
                tr.Spikes = [];
                continue
            end

            channels = [tr.Spikes.Channel];

            ar = AcuteRecording(tr, 'N/A');
            ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
            eu = EphysUnit(ar, readWaveforms=false, cullITI=false, savepath='C:\SERVER\Units\TwoColor_SC\ReverseInjection_Dummy', tr=tr);

            tr.Spikes = [];
            clear eu
        end
    catch ME
        warning('Could not process folder: %s', folders{iSession})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end
