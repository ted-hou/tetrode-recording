%%
load('\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\PawAnalysis\expName.mat')
load('\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\PawAnalysis\frames.mat')


while ~isempty(expNames)
    if ~isempty(frames{1})
        ps = PawnalyzerSession(expNames{1}, frames{1});
        ps.save();
    end
    expNames(1) = [];
    frames(1) = [];
end
clear ps

%%

clear, clc
load('C:\SERVER\PawAnalysis\expName.mat')

frames = cell(1, length(expNames));
for iExp = 1:length(expNames)
    fprintf(1, 'Extracting frames from session %i (%s)...', iExp, expNames{iExp})
    try
        tTic = tic();
        rs = RecordingSession(expNames{iExp});
        frames{iExp} = rs.extractFrames([rs.trials('press').Stop]);
        fprintf('Done (%.1fs)\n', toc(tTic));
    catch ME
        fprintf('\n')
        warning('Error when processing session %i (%s) - this one will be skipped.', iExp, expNames{iExp})
        warning('Error in program %s.\nTraceback (most recent at top):\n%s\nError Message:\n%s', mfilename, getcallstack(ME), ME.message)
    end
end
%%
frames = cat(4, frames{:});
implay(frames)