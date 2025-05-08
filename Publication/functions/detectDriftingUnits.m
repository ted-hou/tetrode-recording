function isDrifting = detectDriftingUnits(eu, varargin)
p = inputParser();
p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
p.addParameter('smoothWindow', 300, @isnumeric) % Smoothing window width in seconds (see smoothdata())
p.addParameter('tolerance', 0.05, @(x) isnumeric(x) && x <= 1 && x >= 0)
p.addParameter('spikeRateThreshold', 15, @isnumeric) % If smoothed spike rate dips below 
p.addParameter('includeITI', false, @islogical)

p.parse(eu, varargin{:})
eu = p.Results.eu;
smoothWindow = p.Results.smoothWindow;
tolerance = p.Results.tolerance;
spikeRateThreshold = p.Results.spikeRateThreshold;
includeITI = p.Results.includeITI;

fig = figure();
ax(1) = subplot(2, 1, 1);
ax(2) = subplot(2, 1, 2);
hold(ax, 'on')
isDrifting = false(size(eu));
for iEu = 1:length(eu)
    if isempty(eu(iEu).SpikeCountTimestamps)
        sampleRate = 1/0.1;
        [sc, t] = eu(iEu).getSpikeCounts(0.1);
        if ~includeITI
            [t, sc, ~] = eu(iEu).cullITIData(t, sc, 'all', 'extendedWindow', [0, 0]);
        end
    else
        sampleRate = 1./eu(iEu).SpikeCountStats.resolution;
        t = eu(iEu).SpikeCountTimestamps;
        sc = eu(iEu).SpikeCounts;
    end
    x = smoothdata(sc.*sampleRate, 'movmean', smoothWindow.*sampleRate);
    if nnz(x < spikeRateThreshold)/length(x) > tolerance
        
        plot(ax(2), t, x)
        isDrifting(iEu) = true;
    else
        plot(ax(1), t, x, Color=[0 0 0 0.2])
    end
end
hold(ax, 'off')
yticks(ax, [0, 15, 40])
xlabel(ax, 'Time in session (s)')
ylabel(ax, 'Spike rate (sp/s)')
title(ax(1), 'Non drifting')
title(ax(2), 'Drifting')

fprintf('%i drifting, %i good, out of %i total.\n', nnz(isDrifting), nnz(~isDrifting), length(eu))