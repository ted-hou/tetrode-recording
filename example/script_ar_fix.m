%% Load
close all
% Load data if necessary
if ~exist('exp_D1', 'var')
    exp_D1 = readdir('C:\SERVER\Experiment_Galvo_D1Cre;DlxFlp;Ai80\AcuteRecording', 'D1');
end
if ~exist('exp_A2A', 'var')
    exp_A2A = readdir('C:\SERVER\Experiment_Galvo_A2ACre\AcuteRecording', 'A2A');
end

for i = 1:length(exp_D1.ar)
    exp_D1.ar(i).importProbeMap(exp_D1.sessionInfo(i).orientation, exp_D1.sessionInfo(i).ml, exp_D1.sessionInfo(i).dv, exp_D1.sessionInfo(i).ap);
    exp_D1.ar(i).save('ar_fixed_', 'C:\SERVER\Experiment_Galvo_D1Cre;DlxFlp;Ai80');
end

for i = 1:length(exp_A2A.ar)
    exp_A2A.ar(i).importProbeMap(exp_A2A.sessionInfo(i).orientation, exp_A2A.sessionInfo(i).ml, exp_A2A.sessionInfo(i).dv, exp_A2A.sessionInfo(i).ap);
    exp_A2A.ar(i).save('ar_fixed_', 'C:\SERVER\Experiment_Galvo_A2ACre');
end

%% Functions
function exp = readdir(fdir, label)
    exp.fdir = fdir;
    exp.ar = AcuteRecording.load(fdir);
    exp.sessionInfo = readSessionInfo(fdir);
    exp.crit = inferCriteria(exp.ar, [0.5, 0.4], 0.01);
    exp.label = label;
end

function sessionInfo = readSessionInfo(fdir)
    fileinfo = dir(sprintf('%s\\sessionInfo_*.mat', fdir));
    for i = 1:length(fileinfo)
        S(i) = load(sprintf('%s\\%s', fileinfo(i).folder, fileinfo(i).name), 'sessionInfo');
    end
    sessionInfo = [S.sessionInfo];
end

function crit = inferCriteria(ar, varargin)
    p = inputParser();
    p.addRequired('AcuteRecording', @(x) isa(x, 'AcuteRecording'));
    p.addOptional('Light', [0.5, 0.4], @isnumeric);
    p.addOptional('Duration', 0.01, @(x) isnumeric(x) && length(x) == 1);
    p.parse(ar, varargin{:});
    ar = p.Results.AcuteRecording;

    crit(length(ar)) = struct('light', [], 'duration', []);
    for i = 1:length(ar)
        critFound = false;
        for light = p.Results.Light(:)'
            if any([ar(i).conditions.light] == light & [ar(i).conditions.duration] == p.Results.Duration)
                crit(i).light = light;
                critFound = true;
                break
            end
        end
        crit(i).duration = p.Results.Duration;
        assert(critFound, 'Criteria (Light=%s, Duration=%s) not matched for experiment %i: %s.', num2str(p.Results.Light), num2str(p.Results.Duration), i, ar(i).expName)
    end
end