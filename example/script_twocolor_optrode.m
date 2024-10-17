%% Load TetrodeRecording objects and convert to EphysUnits
sessionInfo(1) = struct('name', 'daisy19_20240410',     'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(2) = struct('name', 'daisy19_20240411',     'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(3) = struct('name', 'daisy20_20240416',     'ml', -2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(4) = struct('name', 'desmond36_20240419',   'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(5) = struct('name', 'desmond36_20240425',   'ml', +2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');
sessionInfo(6) = struct('name', 'desmond37_20240415',   'ml', -2300, 'ap', 600, 'dv', -3000, 'duraOffset', -950, 'facing', 'left', 'model', '128DN');

for iExp = 2:length(sessionInfo)
    try
        files = dir(sprintf('Y:\\tr_sorted_%s*.mat', sessionInfo(iExp).name));
        files = arrayfun(@(x) [x.folder, x.name], files, UniformOutput=false);
        tr = TetrodeRecording.BatchLoadSimple(files);
        ar = AcuteRecording(tr, 'N/A');
        ar.binMoveResponse(tr, 'none', Window=[-1, 0], Store=true);
        eu = EphysUnit(ar, readWaveforms=true, cullITI=false, savepath='Y:\Units\TwoColor_Optrode', tr=tr, sessionInfo=sessionInfo(iExp));
        clear tr ar eu
    catch
        warning('Failed to process session %i', iExp)
    end
end
clear iExp
%% Load EU
eu = EphysUnit.load('C:\SERVER\Units\TwoColor_Optrode\NonDuplicate_SingleUnit_Good');

%% Make rasters
rd = eu.getRasterData('stimtwocolor', window=[-0.1, 0.4], durErr=1e-2, shutterDelay=0.01);
%% Make rasters and save to disk
fig = figure(Units='inches', Position=[0, 0, 6, 8]);
ax = axes(fig);
if ~exist('Y:\Figures\TwoColor_Optrode', 'dir')
    mkdir('Y:\Figures\TwoColor_Optrode')
end
for iEu = 1:length(eu)
    cla(ax)
    EphysUnit.plotRaster(ax, rd(iEu), xlim=[-0.1, 0.3], sz=2);
    print(fig, sprintf('Y:\\Figures\\TwoColor_Optrode\\%s.png', eu(iEu).getName()), '-dpng');
end

%% Make PE-ISI and save to disk
clear p
p.isiWindow = [-0.5, 0.5];
p.isiRes = 1e-3;
p.xlim = [-0.1, 0.1];
p.path = 'Y:\Figures\TwoColor_Optrode_PEISI\stage';
if ~exist(p.path, 'dir')
    mkdir(p.path)
end
for iEu = 1:length(eu)
% for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[0, 0, 6, 8]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'});
    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = isi;
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes);
        normSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
    end
    imagesc(ax, 1e3*t, [], normSR)
    xlim(ax, p.xlim*1e3)
    clim(ax, [0, 50])
    ax.YAxisLocation = 'right';
    colormap(ax, 'hot')
    h = colorbar(ax, 'westoutside');
    h.Label.String = '\Deltasp/s';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    title(ax, eu(iEu).getName, Interpreter="none")
    print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    close(ax.Parent)
end
%% 'wavelength', 'power', 'location'
for iEu = 1:length(eu)
    ax = axes(figure(Units='inches', Position=[0, 0, 6, 8]));
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'location'});
    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = isi;
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes);
        normSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0), 'omitnan');
    end
    imagesc(ax, 1e3*t, [], normSR)
    xlim(ax, p.xlim*1e3)
    clim(ax, [0, 50])
    ax.YAxisLocation = 'right';
    colormap(ax, 'hot')
    h = colorbar(ax, 'westoutside');
    h.Label.String = '\Deltasp/s';
    yticks(ax, 1:length(groups));
    yticklabels(ax, {groups.label})
    xlabel(ax, 'Time from opto onset (ms)')
    title(ax, eu(iEu).getName, Interpreter="none")
    % print(ax.Parent, sprintf('%s\\%s.png', p.path, eu(iEu).getName), '-dpng')
    % close(ax.Parent)
end


%% {'wavelength', 'power'} plot for all units merged
p.isiWindow = [-0.5, 0.5];
p.isiRes = 1e-3;
p.xlim = [-0.1, 0.1];

GROUPS = cell(length(eu), 1);
ISI = cell(length(eu), 1);
NORMSR = [];
for iEu = 1:length(eu)
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'}, selectBy=struct(wavelength=[], location=[], duration=[], power=[25, 50, 100]*1e-6));
    
    % Make a custom group containing all powers >500uW, group by wavelength
    % groupsHighPower473 = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'}, selectBy=struct(wavelength=473, location=[], duration=[], power=[500, 1000, 2000]*1e-6));
    % groupsHighPower593 = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power'}, selectBy=struct(wavelength=593, location=[], duration=[], power=[500, 1000, 2000]*1e-6));
    % 
    % groups(end + 1) = struct( ...
    %     trials=[groupsHighPower473.trials], ...
    %     pulseIndices=[groupsHighPower473.pulseIndices], ...
    %     label=' 473nm 500+uW', ...
    %     wavelength=473, ...
    %     power=500e-6);
    % groups(end + 1) = struct( ...
    %     trials=[groupsHighPower593.trials], ...
    %     pulseIndices=[groupsHighPower593.pulseIndices], ...
    %     label=' 593nm 500+uW', ...
    %     wavelength=593, ...
    %     power=500e-6);
    % [A, I] = sort([groups.wavelength] + [groups.power]);
    % groups = groups(I);

    isi = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    normSR = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)));
    for iGrp = 1:length(groups)
        [isi(iGrp, :), t] = eu(iEu).getMeanPEISI('stimtwocolor', groups(iGrp).trials, window=p.isiWindow, resolution=p.isiRes, shutterDelay=0.01);
        normSR(iGrp, :) = 1./isi(iGrp, :) - mean(1./isi(iGrp, t<0&t>-0.1), 'omitnan');
    end
    GROUPS{iEu} = groups;
    ISI{iEu} = isi;
    if isempty(NORMSR)
        NORMSR = NaN(length(groups), length(p.isiWindow(1):p.isiRes:p.isiWindow(2)), length(eu));
    end
    NORMSR(:, :, iEu) = normSR;
end


%%
p.fontSize=9;

animalNames = eu.getAnimalName();
c.isA2AChrimsonR = ismember(animalNames, {'daisy19', 'daisy20'});
c.isD1ChR2 = ismember(animalNames, {'desmond36', 'desmond37'});
minDeltaSR = 5;

egUnitNames = { ...
    'daisy19_20240411_Channel78_Unit1', ...
    % 'desmond36_20240419_Channel47_Unit2', ...
    'desmond37_20240415_Channel119_Unit1', ...
    };

SEL = {c.isA2AChrimsonR, c.isD1ChR2};
LABEL = {'A2A-ChrimsonR', 'D1-ChR2'};
YTINC = [50, 100];
fig = figure(Units='inches', Position=[0, 0, 6, 3.5]);

layout.top.h = 3;
layout.bottom.h = 2;

tl = tiledlayout(fig, layout.top.h + layout.bottom.h, 1);
tlTop = tiledlayout(tl, 1, 2);
tlTop.Layout.Tile = 1; tlTop.Layout.TileSpan = [layout.top.h, 1];

tlBottom = tiledlayout(tl, 1, 2);
tlBottom.Layout.Tile = 1 + layout.top.h; tlBottom.Layout.TileSpan = [layout.bottom.h, 1];

for i = 1:2
    iEu = find(strcmpi(eu.getName(), egUnitNames{i}));
    assert(~isempty(iEu) && length(iEu) == 1)
    ax = nexttile(tlTop);
    rd = eu(iEu).getRasterData('stimtwocolor', window=[-0.5, 0.5], durErr=1e-2, shutterDelay=0.01, sort=false);
    groups = eu(iEu).groupTwoColorStimTrials({'wavelength', 'power', 'duration'}, ...
        selectBy=struct(wavelength=[], location=[], duration=[0.2], power=[25, 50, 100]*1e-6));
    hold(ax, 'on')
    everyNth = 1;
    timescale = 1e3;
    sz = 2.5;
    
    iPulseStart = 1;
    sortedI = rd.I;
    hPatch = gobjects(length(groups), 1);
    for iGrp = 1:length(groups)
        inGroup = ismember(rd.I, groups(iGrp).pulseIndices);
        uniqueI = unique(rd.I(inGroup));
        iPulseEnd = iPulseStart + length(uniqueI) - 1;
        newI = iPulseStart:iPulseEnd;
        newI = newI(randperm(length(newI)));
        sortedI(inGroup) = changem(rd.I(inGroup), newI, uniqueI);
    
        switch groups(iGrp).wavelength
            case 473
                color = [0.2, 0.2, 0.8];
            case 593
                color = [0.8, 0.2, 0.2];
        end
        powers = sort(unique([groups.power]), 'ascend');
        iPower = find(powers == groups(iGrp).power);
        alpha = 0.1 + 0.5*((iPower - 1)./(length(powers) - 1));
        partialDesc = sprintf('%inm, %gµW', groups(iGrp).wavelength, groups(iGrp).power*1e6);
    
        scatter(ax, rd.t(inGroup) * timescale, sortedI(inGroup), sz, 'k', 'filled', DisplayName='spikes')
        hPatch(iGrp) = patch(ax, [0, groups(iGrp).duration*1e3, groups(iGrp).duration*1e3, 0], [iPulseStart, iPulseStart, iPulseEnd, iPulseEnd], color, ...
                                FaceAlpha=alpha, EdgeColor=color, EdgeAlpha=alpha, DisplayName=partialDesc);
        iPulseStart = iPulseEnd + 1;
    end
    if i == 2
        lgd = legend(ax, hPatch);
        lgd.Layout.Tile = 'east';
    end
    xlim(ax, [-50, 100])
    ylabel('Trial')
    ylim(ax, [0, iPulseEnd+1])

    yt = 0:YTINC(i):iPulseEnd;
    yt(1) = 1;
    if round(yt(end)./100) == round(iPulseEnd./100)
        yt(end) = iPulseEnd;
    else
        yt(end + 1) = iPulseEnd;
    end
    yticks(ax, yt)
    
    ax.YAxis.Direction = 'reverse';
    title(ax, sprintf('%s (example unit)', LABEL{i}), Interpreter="none")

    % imagesc(ax, 1e3*t, [], NORMSR(:, :, iEu))
    % xlim(ax, [-50, 100])
    % clim(ax, [0, 50])
    % ax.YAxisLocation = 'right';
    % colormap(ax, 'hot')
    % % h = colorbar(ax, 'westoutside');
    % % h.Label.String = '\Deltasp/s';
    % yticks(ax, 1:length(groups));
    % yticklabels(ax, {groups.label})
    % % xlabel(ax, 'Time from opto onset (ms)')
    % title(ax, sprintf('%s (example unit)', LABEL{i}), Interpreter="none")
    % fontsize(ax, p.fontSize, 'points')
end

for i = 1:2
    ax = nexttile(tlBottom);
    selOpsin = SEL{i}(:);
    selResponsive = squeeze(any(abs(max(NORMSR(:, t>0.01&t<0.05, :), [], 2)) > minDeltaSR, 1));
    selResponsive = selResponsive(:);
    imagesc(ax, 1e3*t, [], mean(NORMSR(:, :, selResponsive & selOpsin), 3, 'omitnan'))
    xlim(ax, [-50, 100])
    clim(ax, [0, 50])
    ax.YAxisLocation = 'right';
    colormap(ax, 'hot')
    if i == 2
        h = colorbar(ax, 'westoutside');
        h.Label.String = '\Deltasp/s';
        h.Layout.Tile = 'east';
    end
    yticks(ax, 1:length(GROUPS{1}));
    label = arrayfun(@(grp) sprintf('%inm, %gµW', grp.wavelength, grp.power*1e6), GROUPS{1}, UniformOutput=false);
    yticklabels(ax, {GROUPS{1}.label})
    % xlabel(ax, 'Time from opto onset (ms)')
    title(ax, sprintf('%s (%i/%i responsive)', LABEL{i}, nnz(selResponsive & selOpsin), nnz(selOpsin)), Interpreter="none")
    fontsize(ax, p.fontSize, 'points')
end
xlabel(tlBottom, 'Time from opto onset (ms)', FontSize=p.fontSize)