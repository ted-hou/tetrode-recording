p = inputParser;
			addRequired(p, 'Channels', @isnumeric);
			addOptional(p, 'Reference', 'CueOn', @ischar);
			addOptional(p, 'Event', 'PressOn', @ischar);
			addOptional(p, 'Exclude', 'LickOn', @ischar);
			addParameter(p, 'Clusters', [], @isnumeric);
			addParameter(p, 'MinTrialLength', 0, @isnumeric);
			addParameter(p, 'Bins', 5, @isnumeric);
			addParameter(p, 'BinMethod', 'percentile', @ischar);
			addParameter(p, 'SpikeRateWindow', 100, @isnumeric);
			addParameter(p, 'ExtendedWindow', [0, 0], @isnumeric);
			addParameter(p, 'SelectedSampleIndex', [], @isnumeric);
			addParameter(p, 'XLim', [], @isnumeric);
			addParameter(p, 'LineStyle', '-', @ischar);
			addParameter(p, 'Ax', []);
			parse(p, channels, varargin{:});
			channels 			= p.Results.Channels;
			reference 			= p.Results.Reference;
			event 				= p.Results.Event;
			exclude 			= p.Results.Exclude;
			clusters 			= p.Results.Clusters;
			minTrialLength 		= p.Results.MinTrialLength;
			nBins 				= p.Results.Bins;
			binMethod 			= p.Results.BinMethod;
			spikeRateWindow 	= p.Results.SpikeRateWindow;
			extendedWindow 		= p.Results.ExtendedWindow;
			selectedSampleIndex	= p.Results.SelectedSampleIndex;
			xRange 				= p.Results.XLim;			
			lineStyle			= p.Results.LineStyle;			
			ax 					= p.Results.Ax;

			if ~isempty(reference)
				referenceDisplayName = reference;
				reference = obj.DigitalEvents.(reference);
			else
				reference = [];
			end
			if ~isempty(event)
				eventDisplayName = event;
				event = obj.DigitalEvents.(event);
			else
				event = [];
			end
			if ~isempty(exclude)
				exclude = obj.DigitalEvents.(exclude);
			else
				exclude = [];
			end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			nBins = numel(obj.DigitalEvents.CueOn); % ah addded, delete this!!!
			%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


			reference = sort(reference);
			event = sort(event);
			[reference, event] = TetrodeRecording.FindFirstInTrial(reference, event, exclude);

			% Bin trials according to trial length (t(event) - t(reference))
			trialLength = event - reference;
			switch binMethod
				case 'percentile'
					edges = prctile(trialLength(trialLength > minTrialLength), 0:(100/nBins):100);
				case 'equal'
					edges = linspace(max(minTrialLength, min(trialLength)), max(trialLength), nBins + 1);
				otherwise
					error('Unrecognized bin method.')
			end
			[NTrials, ~, bins] = histcounts(trialLength, edges);

			

			for iChannel = channels
				[sampleIndex, spikes, trials] = obj.GetSpikesByTrial(iChannel, 'Reference', reference, 'Event', event, 'Clusters', clusters,...
					'Window', extendedWindow, 'WindowReference', 'StartAndEnd');

				if ~isempty(selectedSampleIndex)
					selected = ismember(sampleIndex, selectedSampleIndex);
					sampleIndex = sampleIndex(selected);
					spikes 		= spikes(selected);
					trials 		= trials(selected);
				end

				if isempty(ax)
					hAxes = axes(figure());
				else
					hAxes = ax;
					axes(hAxes);
				end
				hold on
				for iBin = 1:nBins
					inBin = ismember(trials, find(bins == iBin));
					cutOff = -edges(iBin);
					spikesRelative = spikes(inBin) - event(trials(inBin)); % Spike times relative to event
					spikesRelative = spikesRelative(spikesRelative > (cutOff + extendedWindow(1)) & spikesRelative < extendedWindow(2));
					% Windows for estimating spike rate
					thisEdges = (cutOff + extendedWindow(1)):(spikeRateWindow/1000):extendedWindow(2);
					if length(thisEdges) < 3
						continue
					end
					% if thisEdges(end) < 0
					% 	thisEdges(end + 1) = 0;
					% end
					thisSpikeRate = histcounts(spikesRelative, thisEdges);
					thisSpikeRate = (1000*thisSpikeRate/spikeRateWindow)/NTrials(iBin);
					thisCenters = (thisEdges(1:end - 1) + thisEdges(2:end))/2;

% 					[thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(iBin, 'Styles', {lineStyle});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AH changed%%%%%%%%%%%%%%%%%%%
                    [thisColor, thisStyle] = TetrodeRecording.GetColorAndStyle(1, 'Styles', {lineStyle});
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
					plot(hAxes, thisCenters, thisSpikeRate, [thisColor, thisStyle],...
						'DisplayName', ['[', num2str(-cutOff, 2), ' s, ', num2str(edges(iBin + 1), 2), ' s] (', num2str(NTrials(iBin)), ' trials)'],...
						'LineWidth', 1);
				end
				xlabel(hAxes, ['Time relative to ', eventDisplayName, ' (s)']);
				ylabel(hAxes, 'Mean firing rate (Spikes/s)')
				legend(hAxes, 'Location', 'Best');
				title(hAxes, 'Peri-event time histogram');
				if ~isempty(xRange)
					xlim(hAxes, xRange);
				end				
				hold off
			end