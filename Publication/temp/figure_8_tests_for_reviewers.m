p.stimThreshold = 2;

% Calculate onset timing
etaRed = SNr_SCRetro.eta.stimRed(["100", "500", "2000+"]);
etaBlue = SNr_SCRetro.eta.stimBlue(["100", "500", "2000"]);
clear tOnset
tOnset.red = NaN(length(eu), 3);
tOnset.blue = NaN(length(eu), 3);
for iPower = 1:length(etaRed)
    t = etaRed(iPower).t;
    X = etaRed(iPower).X(:, t>=0);
    t = t(t>=0);
    for iEu = 1:size(X, 1)
        iStart = strfind([false, X(iEu, :)>2], [0, 1, 1]);
        if isempty(iStart)
            tOnset.red(iEu, iPower) = NaN;
        else
            iStart = iStart(1) - 1;
            if iStart > 0
                tOnset.red(iEu, iPower) = t(iStart);
            else
                tOnset.red(iEu, iPower) = 0;
            end
        end
    end
end

for iPower = 1:length(etaBlue)
    t = etaBlue(iPower).t;
    X = etaBlue(iPower).X(:, t>=0);
    t = t(t>=0);
    for iEu = 1:size(X, 1)
        iStart = strfind([false, X(iEu, :)>p.stimThreshold], [0, 1, 1]);
        if isempty(iStart)
            tOnset.blue(iEu, iPower) = NaN;
        else
            iStart = iStart(1) - 1;
            if iStart > 0
                tOnset.blue(iEu, iPower) = t(iStart);
            else
                tOnset.blue(iEu, iPower) = 0;
            end
        end
    end
end

% Calculate stim response
metaRed = cat(2, SNr_SCRetro.meta.stimRed.values{:});
metaBlue = cat(2, SNr_SCRetro.meta.stimBlue.values{:});
metaRed = array2table(metaRed, VariableNames=SNr_SCRetro.meta.stimRed.keys);
metaBlue = array2table(metaBlue, VariableNames=SNr_SCRetro.meta.stimBlue.keys);

cStim.isRed = metaRed{:, "2000+"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue100 = metaBlue{:, "100"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue = metaBlue{:, "500"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue500 = metaBlue{:, "500"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue2000 = metaBlue{:, "2000"}>p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlueNotRed = cStim.isBlue500 & ~cStim.isRed & SNr_SCRetro.c.isSNr;
cStim.isRedNotBlue = cStim.isRed & ~cStim.isBlue2000 & SNr_SCRetro.c.isSNr;
fprintf("%i (of %i) neurons with DV < -3.8\n", nnz(SNr_SCRetro.c.isSNr), length(SNr_SCRetro.c.isSNr));
fprintf("isRed = %i, isBlue = %i, isBlueNotRed = %i, isRedNotBlue = %i\n", nnz(cStim.isRed), nnz(cStim.isBlue500), nnz(cStim.isBlueNotRed), nnz(cStim.isRedNotBlue));

cStim.isRedNeg = metaRed{:, "2000+"}<-p.stimThreshold & SNr_SCRetro.c.isSNr;
cStim.isBlue500Neg = metaBlue{:, "500"}<-p.stimThreshold & SNr_SCRetro.c.isSNr;
fprintf("isRedNeg = %i, isBlueNeg = %i, isEitherNeg = %i (%.1f%%)\n", nnz(cStim.isRedNeg), nnz(cStim.isBlue500Neg), nnz(cStim.isRedNeg | cStim.isBlue500Neg), nnz(cStim.isRedNeg | cStim.isBlue500Neg)/nnz(SNr_SCRetro.c.isSNr)*100);


%%


negEuIndices = [find(cStim.isBlue500Neg)', find(cStim.isRedNeg)'];

% close all
p.sz = 1.5;
fig = figure(Units='inches', Position=[1, 0.5, 4.7, 10]);
tl = tiledlayout(fig, length(negEuIndices), 2, TileSpacing='tight');
ax = gobjects(length(negEuIndices), 2);
for iUnit = 1:length(negEuIndices)
    ax(iUnit, 1) = nexttile(tl);
    ax(iUnit, 2) = nexttile(tl);
    iEu = negEuIndices(iUnit);
    EphysUnit.plotRaster(ax(iUnit, 1), SNr_SCRetro.rd.stim(iEu), xlim=[-0.02, 0.05], filterByTwoColorConditions=struct(wavelength=[470, 473], power=["*1e6>25"]), mergeStimTrains=true, sz=p.sz);
    EphysUnit.plotRaster(ax(iUnit, 2), SNr_SCRetro.rd.stim(iEu), xlim=[-0.02, 0.05], filterByTwoColorConditions=struct(wavelength=[590, 593, 635], power=["*1e6>25"]), mergeStimTrains=true, sz=p.sz);

    delete(ax(iUnit, 1).Legend);
    delete(ax(iUnit, 2).Legend);
    title(ax(iUnit, 1), sprintf('SNr unit %i', iUnit))
    title(ax(iUnit, 2), '')
    if iUnit <= nnz(cStim.isBlue500Neg)
        title(ax(iUnit, 1), "blue neg", Interpreter='none')
        title(ax(iUnit, 2), "")
    else
        title(ax(iUnit, 2), "red neg", Interpreter='none')
        title(ax(iUnit, 1), "")
    end
end
xlim(ax, [-0.02, 0.05])
xticks(ax, [0, 20, 50]*1e-3)
xticklabels(ax, ["0", "20", "50"])
xlabel(ax, '')
ylabel(ax, '')
xticklabels(ax(1, :), ["", "", ""])
fontsize(ax, p.fontSize, 'points')
xlabel(tl, 'Time from opto onset (ms)', FontSize=p.fontSize)
ylabel(tl, 'Trial', FontSize=p.fontSize)

% ax(1, 2).Legend.Location = 'eastoutside';
% ax(2, 2).Legend.Location = 'eastoutside';
% ax(1, 2).Legend.Position = [0.704491725768321,0.604253468207187,0.234027777777778,0.223958333333333];

% hLetter = text(ax(1, 1), 0, 0, 'd', FontSize=16, FontName='Arial', FontWeight='bold', Units='inches');
% ax(1).Units = 'inches';
% hLetter.HorizontalAlignment = 'right';
% hLetter.VerticalAlignment = 'top';
% hLetter.Position = [-0.4, ax(1).Position(4) + 0.25, 0];
% 
% copygraphics(fig, BackgroundColor='none', ContentType='vector')
% clear i j fig tl ax egUnitNames iUnit iEu hLetter


%% 
cStim.isFastRed = tOnset.red(:, 2) < median(tOnset.red(:, 2), 'omitnan');
cStim.isFastBlue = tOnset.blue(:, 2) < median(tOnset.blue(:, 2), 'omitnan');

sel = cStim.isFastBlue & cStim.isReachSCProjecting;
fprintf("%i/%i (%.1f%%)\n", nnz(sel & SNr_SCRetro.c.isPressDown), nnz(sel), 100*nnz(sel & SNr_SCRetro.c.isPressDown)/nnz(sel))

sel = cStim.isFastRed & cStim.isLickSCProjecting;
fprintf("%i/%i (%.1f%%)\n", nnz(sel & SNr_SCRetro.c.isLickDown), nnz(sel), 100*nnz(sel & SNr_SCRetro.c.isLickDown)/nnz(sel))

