%% 2tgt
p.bootAlpha = 0.01;
p.nboot = 100000;
p.responseWindowDirReach = [-0.1, 0.2];
nUnits = length(euReachDir2Tgt);
boot.press2tgt(1) = struct(target='contra-out', h=NaN(nUnits, 1), muDiffCI=NaN(nUnits, 2), muDiffObs=NaN(nUnits, 1));
boot.press2tgt(2) = struct(target='contra-in', h=NaN(nUnits, 1), muDiffCI=NaN(nUnits, 2), muDiffObs=NaN(nUnits, 1));

p.minNumTrials = 4;
assert(p.minNumTrials == 4)
ETA = trajCombined2tgt.eta;
N = horzcat(ETA.N);
selUnits = all(N >= p.minNumTrials, 2);

c.hasPress2tgt = selUnits;

[boot.press2tgt(1).h(c.hasPress2tgt), boot.press2tgt(1).muDiffCI(c.hasPress2tgt, :), boot.press2tgt(1).muDiffObs(c.hasPress2tgt)] = bootstrapMoveResponse( ...
    euReachDir2Tgt(c.hasPress2tgt), 'press_spontaneous_lateral', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowDirReach, trials=trajCombined2tgt.trials{1}(c.hasPress2tgt), correction=trajCombined2tgt.correction{1}(c.hasPress2tgt));
[boot.press2tgt(2).h(c.hasPress2tgt), boot.press2tgt(2).muDiffCI(c.hasPress2tgt, :), boot.press2tgt(2).muDiffObs(c.hasPress2tgt)] = bootstrapMoveResponse( ...
    euReachDir2Tgt(c.hasPress2tgt), 'press_spontaneous_medial', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
    responseWindow=p.responseWindowDirReach, trials=trajCombined2tgt.trials{2}(c.hasPress2tgt), correction=trajCombined2tgt.correction{2}(c.hasPress2tgt));

fprintf(1, '\nAll done\n')


%% Report bootstraped movement response direction
assert(nnz(isnan(boot.press2tgt(1).h(c.hasPress2tgt))) == 0)
assert(nnz(isnan(boot.press2tgt(2).h(c.hasPress2tgt))) == 0)

c.isPressUp2tgt{1} = boot.press2tgt(1).h == 1 & c.hasPress2tgt;
c.isPressUp2tgt{2} = boot.press2tgt(2).h == 1 & c.hasPress2tgt;

c.isPressDown2tgt{1} = boot.press2tgt(1).h == -1 & c.hasPress2tgt;
c.isPressDown2tgt{2} = boot.press2tgt(2).h == -1 & c.hasPress2tgt;

c.isPressResponsive2tgt{1} = c.isPressUp2tgt{1} | c.isPressDown2tgt{1};
c.isPressResponsive2tgt{2} = c.isPressUp2tgt{2} | c.isPressDown2tgt{2};

fprintf(1, ['%g total SNr units (baseline spike rate > %g):\n' ...
    '\t%g with %d+ press trials (2tgt, either direction) ;\n' ...
    '\t%g up, %g down for contra-out;\n' ...
    '\t%g up, %g down for contra-in;\n' ...
    '\t%g up/up, %g down/down;\n' ...
    '\t%g up/down, %g down/up;\n' ...
    ], ...
    length(euReachDir2Tgt), p.minSpikeRate, ...
    nnz(c.hasPress2tgt), p.minNumTrials, ...
    nnz(c.isPressUp2tgt{1}), nnz(c.isPressDown2tgt{1}), ...
    nnz(c.isPressUp2tgt{2}), nnz(c.isPressDown2tgt{2}), ...
    nnz(c.isPressUp2tgt{1} & c.isPressUp2tgt{2}), nnz(c.isPressDown2tgt{1} & c.isPressDown2tgt{2}), ...
    nnz(c.isPressUp2tgt{1} & c.isPressDown2tgt{2}), nnz(c.isPressDown2tgt{1} & c.isPressUp2tgt{2}))

%% 4tgt
p.bootAlpha = 0.01;
p.nboot = 100000;
p.responseWindowDirReach = [-0.1, 0.2];
nUnits = length(euReachDir4Tgt);
boot.press4tgt(1) = struct(target='contra-out', h=NaN(nUnits, 1), muDiffCI=NaN(nUnits, 2), muDiffObs=NaN(nUnits, 1));
boot.press4tgt(2) = struct(target='contra-front', h=NaN(nUnits, 1), muDiffCI=NaN(nUnits, 2), muDiffObs=NaN(nUnits, 1));
boot.press4tgt(3) = struct(target='contra-in', h=NaN(nUnits, 1), muDiffCI=NaN(nUnits, 2), muDiffObs=NaN(nUnits, 1));
boot.press4tgt(4) = struct(target='ipsi-front', h=NaN(nUnits, 1), muDiffCI=NaN(nUnits, 2), muDiffObs=NaN(nUnits, 1));

p.minNumTrials4tgt = 4;
assert(p.minNumTrials4tgt == 4)

N = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
IPAW = [1, 1, 1, 3];
assert(p.minNumTrials4tgt == 4)
selUnits = N{1, IPAW(1)} >= p.minNumTrials4tgt & N{2, IPAW(2)} >= p.minNumTrials4tgt & N{3, IPAW(3)} >= p.minNumTrials4tgt & N{4, IPAW(4)} >= p.minNumTrials4tgt;

c.hasPress4tgt = selUnits;

for iTarget = 1:4
    [boot.press4tgt(iTarget).h(c.hasPress4tgt), boot.press4tgt(iTarget).muDiffCI(c.hasPress4tgt, :), boot.press4tgt(iTarget).muDiffObs(c.hasPress4tgt)] = bootstrapMoveResponse( ...
        euReachDir4Tgt(c.hasPress4tgt), 'press_spontaneous_lateral', nboot=p.nboot, alpha=p.bootAlpha, withReplacement=false, oneSided=false, ...
        responseWindow=p.responseWindowDirReach, trials=trajCombined.trials{iTarget, IPAW(iTarget)}(c.hasPress4tgt), correction=trajCombined.correction{iTarget, IPAW(iTarget)}(c.hasPress4tgt));
end

fprintf(1, '\nAll done\n')


%% Report bootstraped movement response direction
assert(nnz(isnan(boot.press4tgt(1).h(c.hasPress4tgt))) == 0)
assert(nnz(isnan(boot.press4tgt(2).h(c.hasPress4tgt))) == 0)
assert(nnz(isnan(boot.press4tgt(3).h(c.hasPress4tgt))) == 0)
% assert(nnz(isnan(boot.press4tgt(4).h(c.hasPress4tgt))) == 0)

for iTarget = 1:4
    c.isPressUp4tgt{iTarget} = boot.press4tgt(iTarget).h == 1 & c.hasPress4tgt;
    c.isPressDown4tgt{iTarget} = boot.press4tgt(iTarget).h == -1 & c.hasPress4tgt;
    c.isPressResponsive4tgt{iTarget} = c.isPressUp4tgt{iTarget} | c.isPressDown4tgt{iTarget};
end
fprintf('%g total SNr units with %i+ trials (baseline spike rate > %g):\n', nnz(c.hasPress4tgt), p.minNumTrials4tgt, p.minSpikeRate)
for iTarget = 1:4
    fprintf('target=%s, nUp=%i, nDown=%i;\n', boot.press4tgt(iTarget).target, nnz(c.isPressUp4tgt{iTarget}), nnz(c.isPressDown4tgt{iTarget}))
end

fprintf('contra-out vs contra-in, up/up=%i, down/down=%i;\n', nnz(c.isPressUp4tgt{1} & c.isPressUp4tgt{3}), nnz(c.isPressDown4tgt{1} & c.isPressDown4tgt{3}))
fprintf('contra-front vs ipsi-front, up/up=%i, down/down=%i;\n', nnz(c.isPressUp4tgt{2} & c.isPressUp4tgt{4}), nnz(c.isPressDown4tgt{2} & c.isPressDown4tgt{4}))
