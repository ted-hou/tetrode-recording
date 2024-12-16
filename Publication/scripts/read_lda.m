%% Read result summaries
load('C:\SERVER\Units\lda_pressVsLick_20241212.mat')
load('C:\SERVER\Units\lda_reach2tgt_20241213.mat')

%% (SLOW) Read full results containing model predictions from each permutation
% % 99CI is saved in the summaries, reload this if we need to calculate new CIs
% load('C:\SERVER\Units\lda_pressVsLick_fullBootData_20241212.mat')
% load('C:\SERVER\Units\lda_reach2tgt_fullBootData_20241213.mat')


%% Quick summary (99% CI, mean) of bootstrap for lick vs reach
% clear dfBootStats;
% 
% Y = cell(length(sr), 1);
% for iExp = 1:length(sr)
%     nPress = size(resp(iExp).press, 1);
%     nLick = size(resp(iExp).lick, 1);
%     Y{iExp} = vertcat(repmat("press", [nPress, 1]), repmat("lick", [nLick, 1]));
% end
% Y = cat(1, Y{:});
% 
% isPress = Y == "press";
% isLick = Y == "lick";
% 
% dfBootStats.press.X = transpose(squeeze(mean(dfBoot(isPress, :, :), 1, 'omitnan')));
% dfBootStats.press.mu = mean(dfBootStats.press.X, 1, 'omitnan');
% dfBootStats.press.ci = quantile(dfBootStats.press.X, [0.01, 0.99], 1);
% 
% dfBootStats.lick.X = transpose(squeeze(mean(dfBoot(isLick, :, :), 1, 'omitnan')));
% dfBootStats.lick.mu = mean(dfBootStats.lick.X, 1, 'omitnan');
% dfBootStats.lick.ci = quantile(dfBootStats.lick.X, [0.01, 0.99], 1);
% 
% dfBootStats.all.X = transpose(squeeze(mean(dfBoot, 1, 'omitnan')));
% dfBootStats.all.mu = mean(dfBootStats.all.X, 1, 'omitnan');
% dfBootStats.all.ci = quantile(dfBootStats.all.X, [0.01, 0.99], 1);
% 
% % save('C:\SERVER\Units\lda_pressVsLick_20241212.mat', 'pLDA', 'likelihood', 'sr', 'resp', 't', 'goodExpNames', 'nUnits', 'dfBootStats')


%% Quick summary (99% CI, mean) of bootstrap for reach 2tgt
% clear dfBootStats2tgt;
% 
% Y = cell(length(sr2tgt), 1);
% for iExp = 1:length(sr2tgt)
%     nContraOut = size(resp2tgt(iExp).contraOut, 1);
%     nContraIn = size(resp2tgt(iExp).contraIn, 1);
%     Y{iExp} = vertcat(repmat("lateral", [nContraOut, 1]), repmat("medial", [nContraIn, 1]));
% end
% Y = cat(1, Y{:});
% 
% isContraOut = Y == "lateral";
% isContraIn = Y == "medial";
% 
% dfBootStats2tgt.contraOut.X = transpose(squeeze(mean(dfBoot2tgt(isContraOut, :, :), 1, 'omitnan')));
% dfBootStats2tgt.contraOut.mu = mean(dfBootStats2tgt.contraOut.X, 1, 'omitnan');
% dfBootStats2tgt.contraOut.ci = quantile(dfBootStats2tgt.contraOut.X, [0.01, 0.99], 1);
% 
% dfBootStats2tgt.contraIn.X = transpose(squeeze(mean(dfBoot2tgt(isContraIn, :, :), 1, 'omitnan')));
% dfBootStats2tgt.contraIn.mu = mean(dfBootStats2tgt.contraIn.X, 1, 'omitnan');
% dfBootStats2tgt.contraIn.ci = quantile(dfBootStats2tgt.contraIn.X, [0.01, 0.99], 1);
% 
% dfBootStats2tgt.all.X = transpose(squeeze(mean(dfBoot2tgt, 1, 'omitnan')));
% dfBootStats2tgt.all.mu = mean(dfBootStats2tgt.all.X, 1, 'omitnan');
% dfBootStats2tgt.all.ci = quantile(dfBootStats2tgt.all.X, [0.01, 0.99], 1);
% 
% % save('C:\SERVER\Units\lda_reach2tgt_20241213.mat', 'pLDA', 'likelihood2tgt', 'sr2tgt', 'resp2tgt', 't', 'expNames2tgt', 'euExpIndex', 'goodExpIndices2tgt', 'nUnits2tgt', 'dfBootStats2tgt')
% 
% clear Y iExp nContraOut nContraIn Y isContraOut isContraIn
