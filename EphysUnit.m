
classdef EphysUnit < handle
	properties
		ExpName = ''
		Channel
		Unit
		MeanSpikeRate
		MeanSpikeRateTimestamps
		SpikeRate
		SpikeRateTimestamps
		Kernel
		KernelTimestamps
		NumTrials
		SpikeTimes
		PressTimes
		CueTimes
	end

	methods
		%% EhpysUnit: Constructor from PETH struct
		function obj = EphysUnit(varargin)
			if nargin == 0
				return
			end

			p = inputParser();
			addOptional(p, 'PETH', struct(), @isstruct);
			addOptional(p, 'kernel', 'gaussian', @ischar);
			addParameter(p, 'lambda1', 5, @isnumeric);
			addParameter(p, 'lambda2', 10, @isnumeric);
			addParameter(p, 'sigma', 0.1, @isnumeric);
			addParameter(p, 'resolution', 0.001, @isnumeric);
			addParameter(p, 'kernelWindow', 1, @isnumeric);
			addParameter(p, 'trialWindow', [-4, 0], @isnumeric);
			parse(p, varargin{:});
			PETH = p.Results.PETH;

			% Initialize obj array
			obj(length(PETH)) = EphysUnit();

			% Populate objects
			for i = 1:length(PETH)
				obj(i).ExpName = PETH(i).ExpName;
				obj(i).Channel = PETH(i).Channel;
				obj(i).Unit = PETH(i).Cluster;
				tr = TetrodeRecording.BatchLoad({obj(i).ExpName});
				obj(i).loadSingleTrialData(tr, PETH(i).Channel, PETH(i).Cluster);
				[obj(i).SpikeRate, obj(i).SpikeRateTimestamps, obj(i).Kernel, obj(i).KernelTimestamps] = obj(i).spikeConv(p.Results.kernel, 'sigma', p.Results.sigma, 'lambda1', p.Results.lambda1, 'lambda2', p.Results.lambda2, 'resolution', p.Results.resolution, 'kernelWindow', p.Results.kernelWindow, 'trialWindow', p.Results.trialWindow);
			end
		end

		% Convolve with kernel: exponential, y(t) = e^(-lambda1*t) - e^(-lambda2*t). b > a, or gaussian: N(0, sigma)
		function varargout = spikeConv(obj, varargin)
			p = inputParser();
			addOptional(p, 'kernel', 'gaussian', @ischar);
			addParameter(p, 'lambda1', 5, @isnumeric);
			addParameter(p, 'lambda2', 10, @isnumeric);
			addParameter(p, 'sigma', 0.1, @isnumeric);
			addParameter(p, 'resolution', 0.001, @isnumeric);
			addParameter(p, 'kernelWindow', 1, @isnumeric);
			addParameter(p, 'trialWindow', [-4, 0], @isnumeric);
			parse(p, varargin{:});
			lambda1 = p.Results.lambda1;
			lambda2 = p.Results.lambda2;
			sigma = p.Results.sigma;
			kernelType = p.Results.kernel;
			resolution = p.Results.resolution;
			kernelWindow = p.Results.kernelWindow;
			trialWindow = p.Results.trialWindow;

			switch kernelType
				case 'exponential'
					kernelWindow = [0, kernelWindow];
					tKernel = kernelWindow(1):resolution:kernelWindow(2);
					kernel = exp(-lambda1*tKernel) - exp(-lambda2*tKernel);
				case 'gaussian'
					kernelWindow = kernelWindow * [-0.5, 0.5];
					tKernel = kernelWindow(1):resolution:kernelWindow(2);
					kernel = normpdf(tKernel, 0, sigma);
			end
			kernel = kernel / sum(kernel) / resolution;

			% Preallocate spikeDensity x trial
			trialWindow = trialWindow + diff(kernelWindow) * [-0.5, 0.5];
			t = trialWindow(1) + diff(kernelWindow)/2 : resolution : trialWindow(2) - diff(kernelWindow)/2;
			numSamples = length(t) - 1;

			for i = 1:length(obj)
				spikeDensity = zeros(numSamples, obj(i).NumTrials);

				% Do the convolution
				for iTrial = 1:obj(i).NumTrials
					tGlobal = obj(i).PressTimes(iTrial) + trialWindow(1) : resolution : obj(i).PressTimes(iTrial) + trialWindow(2);
					spikeCounts = histcounts(obj(i).SpikeTimes, tGlobal);
					spikeDensity(:, iTrial) = conv(spikeCounts, kernel, 'valid');
				end

				obj(i).SpikeRate = spikeDensity;
				obj(i).SpikeRateTimestamps = t(1:end-1);
				obj(i).Kernel = kernel;
				obj(i).KernelTimestamps = tKernel;
			end
			if length(obj) == 1
				varargout = {spikeDensity, t(1:end-1), kernel, tKernel};
			else
				varargout = {};
			end
		end

		function plot(obj, varargin)
			p = inputParser();
			addParameter(p, 'maxShownResults', Inf, @isnumeric);
			addParameter(p, 'stagger', 0.3, @isnumeric);
			addParameter(p, 'orderBy', 'moveLatency', @ischar); % 'moveLatency', 'random', 'original'
			parse(p, varargin{:});
			maxShownResults = p.Results.maxShownResults;
			stagger = p.Results.stagger;
			orderBy = p.Results.orderBy;

			for i = 1:length(obj)
				if maxShownResults > 0
					numResults = min(maxShownResults, obj(i).NumTrials);

					f = figure('Units', 'normalized', 'OuterPosition', [0, 0, 0.33, 1]);
					ax1 = subplot(1, 2, 1);
					ax2 = subplot(1, 2, 2);

					% Use the same seed every time. Plot 50 random curves
					staggerHeight = stagger * median(obj(i).MeanSpikeRate);
					staggerOffset = repmat(linspace(0, staggerHeight*(numResults - 1), numResults), length(obj(i).SpikeRateTimestamps), 1);

					switch lower(orderBy)
						case 'original'
							order = 1:numResults;
						case 'random'
							rng(12345); % Use the same seed
							order = randi(obj(i).NumTrials, numResults, 1);
						case 'movelatency'
							latency = obj(i).PressTimes - obj(i).CueTimes;
							[latency, order] = sort(latency, 'descend');
							order = order(1:numResults);
							latency = latency(1:numResults);
						otherwise
							error('Unsupported "orderBy" argument.');
					end

					plot(ax1, obj(i).SpikeRateTimestamps, obj(i).SpikeRate(:, order) + staggerOffset, 'k');
					title(ax1, sprintf('Single trial spike rates\n(%i trials shown)', numResults), 'Interpreter', 'none')
					xlabel(ax1, 'Time (s) to lever press')
					ylabel(ax1, 'a.u.')
					yticks(ax1, [])
					yticklabels(ax1, [])

					% Plot trial start times
					if strcmpi(orderBy, 'moveLatency')
						hold(ax1, 'on');
						plot(ax1, -latency, staggerOffset(1, :), 'r')
						hold(ax1, 'off')
					end
					xlim(ax1, obj(i).SpikeRateTimestamps([1, length(obj(i).SpikeRateTimestamps)]))

					% Plot trial average
					plot(ax2, obj(i).MeanSpikeRateTimestamps, obj(i).MeanSpikeRate, 'k');
					title(ax2, sprintf('Mean spike rate\n(%i trials)', obj(i).NumTrials), 'Interpreter', 'none')
					xlim(ax2, ax1.XLim);
					xlabel(ax2, 'Time (s) to lever press')
					ylabel(ax2, 'sp/s')

					suptitle(obj(i).getUnitDisplayName());
				end
			end
		end
		
		function plotKernel(obj)
			for i = 1:length(obj)
				f = figure();
				ax = axes(f);
				plot(ax, obj(i).KernelTimestamps, obj(i).Kernel);
				title(ax, sprintf('Kernel\n%s', obj(i).getUnitDisplayName()), 'Interpreter', 'none')
				xlabel(ax, 'Time (s)')
				ylabel(ax, 'Value')
			end
		end

		function t = toTable(obj)
			t = struct2table(toStruct(obj));
		end

		function s = toStruct(obj)
			w = warning('off','MATLAB:structOnObject');
			s = arrayfun(@struct, obj);
			warning(w);
		end
	end

	% Private methods
	methods (Access = {})
		%% getSingleTrialData: Get single trial data
		function loadSingleTrialData(obj, tr, channel, cluster)
			[obj.MeanSpikeRate, obj.MeanSpikeRateTimestamps, obj.NumTrials, obj.SpikeRate, obj.SpikeTimes, obj.PressTimes, obj.CueTimes] = tr.PETHistCounts(...
								channel, 'Cluster', cluster,...
								'Event', 'PressOn', 'Exclude', 'LickOn',...
								'TrialLength', 6, 'ExtendedWindow', 1, 'SpikeRateWindow', 100);
		end

		function displayName = getUnitDisplayName(obj)
			displayName = sprintf('%s Chn%i Unit%i', obj.ExpName, obj.Channel, obj.Unit);
		end
	end
end