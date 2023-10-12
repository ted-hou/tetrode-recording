%% Load whitelisted BehaviorSession objects.
whitelist = ...
[ ...
    {'daisy2'   }; ...
    {'daisy3'   }; ...
    {'daisy4'   }; ...
    {'daisy5'   }; ...
    {'daisy7'   }; ...
    {'daisy8'   }; ...
    {'daisy9'   }; ...
    {'daisy10'  }; ...
    {'daisy12'  }; ...
    {'daisy13'  }; ...
    {'daisy14'  }; ...
    {'daisy15'  }; ...
    {'daisy16'  }; ...
    {'desmond10'}; ...
    {'desmond11'}; ...
    {'desmond12'}; ...
    {'desmond13'}; ...
    {'desmond15'}; ...
    {'desmond16'}; ...
    {'desmond17'}; ...
    {'desmond18'}; ...
    {'desmond20'}; ...
    {'desmond21'}; ...
    {'desmond22'}; ...
    {'desmond23'}; ...
    {'desmond24'}; ...
    {'desmond25'}; ...
    {'desmond26'}; ...
    {'desmond27'}; ...
];

S = load('C:\SERVER_PRIVATE\data\bs.mat');
bs = S.bs;
clear S
bs = cellfun(@(x) x(x.isvalid()), bs, 'UniformOutput', false);
animalNames = cellfun(@(bs) bs.animalName, bs, UniformOutput=false);
isWhitelisted = ismember(animalNames, whitelist);
bs = bs(isWhitelisted);
nAnimals = length(bs);

% Step 3. Process data
% Extract data
interval = cellfun(@(bs) arrayfun(@(bs) mean(bs.getParams('INTERVAL_MIN')), bs), bs, 'UniformOutput', false);
maxInterval = cellfun(@(x) max(x), interval);
lastInterval = cellfun(@(x) x(end), interval);
nPressTrials = cellfun(@(bs) bs.countTrials('press'), bs, UniformOutput=false);
nPressTrialsAggr = cat(1, nPressTrials{:});
nPressTrialsAggr(nPressTrialsAggr == 0) = [];
fprintf(1, '%g animals, %g to %g press trials, median %g.\n', length(bs), prctile(nPressTrialsAggr, 5), prctile(nPressTrialsAggr, 95), median(nPressTrialsAggr))

intervalAggr = cat(1, interval{:});

pressTimes = cellfun(@(x) arrayfun(@(y) y.getEventTimesRelative('LEVER_PRESSED', RequireTrialType='lever'), x, 'UniformOutput', false, 'ErrorHandler', @(varargin) []), bs, 'UniformOutput', false);
pressTimesAggr = cat(1, pressTimes{:});
nPressTrialsActual = nonzeros(cellfun(@length, pressTimesAggr));

nPressTrialsEarly = cellfun(@(bs) bs.countTrials('press', 'early'), bs, UniformOutput=false);
nPressTrialsCorrect = cellfun(@(bs) bs.countTrials('press', 'correct'), bs, UniformOutput=false);
nPressTrialsNoMove = cellfun(@(bs) bs.countTrials('press', 'nomove'), bs, UniformOutput=false);

nPressTrialsEarly = nonzeros(cat(1, nPressTrialsEarly{:}));
nPressTrialsCorrect = nonzeros(cat(1, nPressTrialsCorrect{:}));
nPressTrialsNoMove = nonzeros(cat(1, nPressTrialsNoMove{:}));

fprintf(1, '%g to %g early press trials, median %g.\n', prctile(nPressTrialsEarly, 5), prctile(nPressTrialsEarly, 95), median(nPressTrialsEarly))
fprintf(1, '%g to %g correct press trials, median %g.\n', prctile(nPressTrialsCorrect, 5), prctile(nPressTrialsCorrect, 95), median(nPressTrialsCorrect))
fprintf(1, '%g to %g nomove press trials, median %g.\n', prctile(nPressTrialsNoMove, 5), prctile(nPressTrialsNoMove, 95), median(nPressTrialsNoMove))


% Extract first press/lick times
minTrials = 30;
nPressTrials = cellfun(@(x) x.countTrials('press'), bs, 'UniformOutput', false);
nLickTrials = cellfun(@(x) x.countTrials('lick'), bs, 'UniformOutput', false);
selPress = cellfun(@(x) x>minTrials, nPressTrials, UniformOutput=false);
selLick = cellfun(@(x) x>minTrials, nLickTrials, UniformOutput=false);
bsPress = cellfun(@(x, y) x(y), bs, selPress, UniformOutput=false);
bsLick = cellfun(@(x, y) x(y), bs, selLick, UniformOutput=false);
pressTimes = cellfun(@(x) arrayfun(@(y) y.getEventTimesRelative('LEVER_PRESSED', RequireTrialType='lever'), x, 'UniformOutput', false, 'ErrorHandler', @(varargin) []), bsPress, 'UniformOutput', false);
lickTimes = cellfun(@(x) arrayfun(@(y) y.getEventTimesRelative('LICK', RequireTrialType='lick'), x, 'UniformOutput', false, 'ErrorHandler', @(varargin) []), bsLick, 'UniformOutput', false);


% Collect data from certain days
daysPress = 1:50;
daysLick = 1:30;
ndaysPress = length(daysPress);
ndaysLick = length(daysLick);
ndays = max(ndaysPress, ndaysLick);
pt = cell(nAnimals, ndays);
lt = cell(nAnimals, ndays);
for ia = 1:nAnimals
    for id = 1:ndays
        try
            pt{ia, id} = pressTimes{ia}{daysPress(id)};
        catch
            pt{ia, id} = [];
        end
    end
end
for ia = 1:nAnimals
    for id = 1:ndays
        try
            lt{ia, id} = lickTimes{ia}{daysLick(id)};
        catch
            lt{ia, id} = [];
        end
    end
end
[mPress, nSessionsPress] = max(cellfun(@isempty, pt), [], 2);
[mLick, nSessionsLick] = max(cellfun(@isempty, lt), [], 2);
nSessionsPress = nSessionsPress - 1;
nSessionsLick = nSessionsLick - 1;
nSessionsPress(mPress == 0 & nSessionsPress == 0) = size(pt, 2);
nSessionsLick(mLick == 0 & nSessionsLick == 0) = size(lt, 2);

nAnimalsPress = nnz(nSessionsPress ~= 0);
nAnimalsLick = nnz(nSessionsLick ~= 0);

ptcat = cell(1, ndaysPress);
ltcat = cell(1, ndaysLick);
for id = 1:ndaysPress
    ptcat{id} = cat(1, pt{:, id});
end
for id = 1:ndaysLick
    ltcat{id} = cat(1, lt{:, id});
end
nTrialsPress = cellfun(@length, ptcat);
nTrialsLick = cellfun(@length, ltcat);