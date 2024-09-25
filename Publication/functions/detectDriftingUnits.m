function isDrifting = detectDriftingUnits(eu, varargin)
p = inputParser();
p.addRequired('eu', @(x) isa(x, 'EphysUnit'))
p.addParameter('smoothWindow', 300, @isnumeric) % Smoothing window width in seconds (see smoothdata())
p.addParameter('tolerance', 0.05, @(x) isnumeric(x) && x <= 1 && x >= 0)
p.addParameter('spikeRateThreshold', 15, @isnumeric) % If smoothed spike rate dips below 

p.parse(eu, varargin{:})
eu = p.Results.eu;
smoothWindow = p.Results.smoothWindow;
tolerance = p.Results.tolerance;
spikeRateThreshold = p.Results.spikeRateThreshold;


ax = axes(figure);
hold(ax, 'on')
isDrifting = false(size(eu));
for iEu = 1:length(eu)
    sampleRate = 1./eu(iEu).SpikeCountStats.resolution;
    t = eu(iEu).SpikeCountTimestamps;
    x = smoothdata(eu(iEu).SpikeCounts.*sampleRate, 'movmean', smoothWindow.*sampleRate);
    if nnz(x < spikeRateThreshold)/length(x) > tolerance
        plot(t, x)
        isDrifting(iEu) = true;
    else
        plot(t, x, Color=[0 0 0 0.2])
    end
end
hold(ax, 'off')
yticks(ax, [0, 15, 40])
xlabel('Time in session (s)')
ylabel('Spike rate (sp/s)')

fprintf('%i drifting, %i good, out of %i total.\n', nnz(isDrifting), nnz(~isDrifting), length(eu))