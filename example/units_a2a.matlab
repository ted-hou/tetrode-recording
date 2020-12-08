% Channel id's are indices for ptr.SelectedChannels. Not actual channel id.

batchPlotList = {...
	% Depth 3.8
	'desmond20', 20201022, 2, 1;... % Big unit, flat
	'desmond20', 20201022, 2, 2;... % Big unit, flat
	'desmond20', 20201022, 6, 1;... % Big unit, very late off, fast excitation
	'desmond20', 20201023, 1, 1;... % Small unit, off into late on, fast excitation
	'desmond20', 20201023, 2, 1;... % Small unit, flat, fast excitation
	'desmond20', 20201023, 4, 1;... % Big unit, late off, fast excitation
	'desmond20', 20201023, 4, 2;... % Big unit, 15Hz/flat, no effect
	'desmond20', 20201023, 7, 1;... % Big unit, off, fast excitation
	'desmond20', 20201026, 4, 1;... % Big unit, late off, fast excitation
	'desmond20', 20201026, 4, 2;... % Big unit, 15Hz/flat, no effect
	'desmond20', 20201026, 8, 1;... % Big unit, late off, fast excitation

	% Depth 4.12
	'desmond20', 20201109, 1, 1;... % Medium unit, on, fast excitation maybe?
	'desmond20', 20201109, 1, 2;... % Big unit, on, fast excitation into inhibition
	'desmond20', 20201109, 3, 1;... % Big unit, off into very late on, fast excitation
	'desmond20', 20201109, 4, 1;... % Big unit, early on, no effect maybe?
	'desmond20', 20201109, 4, 2;... % Small unit, early on, fast excitation
	'desmond20', 20201109, 5, 1;... % Big unit, early on, fast excitation
	'desmond20', 20201109, 5, 2;... % Small unit, very early on, fast excitation
	'desmond20', 20201109, 6, 1;... % Medium unit, early on, fast excitation into inhibition, source duplicate maybe?
	'desmond20', 20201109, 6, 2;... % Medium unit, early on, fast excitation
	'desmond20', 20201109, 7, 1;... % Big unit, on, fast excitation maybe?
	'desmond20', 20201109, 7, 2;... % Big unit, on, no effect maybe?
	'desmond20', 20201109, 8, 1;... % Small unit, on, fast excitation
	'desmond20', 20201109, 11, 1;... % Small unit, on, fast excitation into inhibition, duplicate maybe?
	'desmond20', 20201109, 12, 1;... % Medium unit, on, fast excitation into inhibition, duplicate maybe?
	'desmond20', 20201109, 12, 1;... % Small unit, on, fast excitation
	'desmond20', 20201109, 13, 1;... % Medium unit, on, fast excitation
	'desmond20', 20201109, 15, 1;... % Medium unit, on, fast excitation
	'desmond20', 20201109, 15, 2;... % Small unit, on, fast excitation into inhibition, duplicate maybe?

	% Depth 4.28
	'desmond20', 20201116, 3, 1;... % Big unit, on, no effect maybe?
	'desmond20', 20201116, 5, 1;... % Medium unit, flat into very late off, excitation
	'desmond20', 20201116, 8, 1;... % Small unit, on, excitation
	'desmond20', 20201116, 9, 1;... % Big unit, on into late off, excitation
	'desmond20', 20201116, 10, 1;... % Small unit, on, excitation
	'desmond20', 20201116, 10, 1;... % Small unit, on, excitation
	'desmond20', 20201118, 2, 1;... % Big unit, off, excitation
	'desmond20', 20201118, 4, 1;... % Big unit, on, slow inhibition maybe?, source duplicate maybe?
	'desmond20', 20201118, 4, 2;... % Big unit, late off, excitation, source duplicate maybe?
	'desmond20', 20201118, 5, 1;... % Big unit, on, slow inhibition maybe?, duplicate maybe?
	'desmond20', 20201118, 5, 2;... % Big unit, late off, excitation, duplicate maybe?
	'desmond20', 20201118, 8, 1;... % Big unit, very late on, slow inhibition maybe?

	% Depth 4.44
	'desmond20', 20201120, 2, 1;... % Medium unit, on, excitation
	'desmond20', 20201201, 2, 1;... % Big unit, on, excitation
	'desmond20', 20201201, 4, 1;... % Big unit, on, excitation
	'desmond20', 20201201, 5, 1;... % Big unit, on, excitation
	'desmond20', 20201201, 5, 2;... % Big unit, on, excitation
	'desmond20', 20201201, 7, 1;... % Big unit, off, no effect
	'desmond20', 20201201, 8, 1;... % Big unit, on, excitation

	% Depth 4.52
	'desmond20', 20201202, 2, 1;... % Big unit, on, excitation
	'desmond20', 20201202, 3, 1;... % Big unit, on, excitation
	'desmond20', 20201202, 4, 1;... % Big unit, on, excitation maybe?
	'desmond20', 20201202, 4, 2;... % Big unit, on into off, excitation
	'desmond20', 20201202, 5, 2;... % Medium unit, on, excitation maybe?
	'desmond20', 20201202, 7, 1;... % Medium unit, on late, no effect
	'desmond20', 20201202, 8, 1;... % Big unit, on late, excitation
	'desmond20', 20201202, 8, 1;... % Big unit, on, no effect
	'desmond20', 20201202, 8, 2;... % Small unit, on, excitation
	'desmond20', 20201202, 9, 1;... % Small unit, on, excitation
	'desmond20', 20201202, 9, 2;... % Medium unit, on, excitation
	'desmond20', 20201202, 10, 1;... % Medium unit, on, excitation
	'desmond20', 20201202, 11, 2;... % Small unit, on into off, excitation

	};

% 
sessions = {...
	'desmond20', 20201022;...
	'desmond20', 20201023;...
	'desmond20', 20201026;...
	'desmond20', 20201027;... % All short trials

};

needsManualSorting = {...
	'desmond20', 20201116, 2;...
	'desmond20', 20201116, 8;...
	'desmond20', 20201116, 10;...
	'desmond20', 20201118, 7;...
	'desmond20', 20201201, 4;..

};