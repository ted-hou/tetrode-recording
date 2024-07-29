%% Load units
eu = EphysUnit.load('\\research.files.med.harvard.edu\neurobio\NEUROBIOLOGY SHARED\Assad Lab\Lingfeng\Data\Units\acute_3cam_reach_direction_2tgts\SingleUnits_NonDuplicate');

%% Plot exponential pdf with lambda=1/20 (mean=1/lambda)
lambda = 1/20;
t = 0:0.1:30;
p = lambda*exp(-lambda*t);

fig = figure();
ax = axes(fig);
hold(ax, 'on')
plot(ax, t, p.*length(ipi))
xlabel(ax, 'Timeout (s)')
ylabel(ax, 'Probability')


%% Tally Inter-press intervals from all 10 sessions, make histogram
expNames = {eu.ExpName};
[~, ia] = unique(expNames);

ipi = [];
for iEu = ia'
    disp(iEu)
    pressTimes = eu(iEu).EventTimes.Press;
    ipi = [ipi, diff(pressTimes)];
end

fig = figure();
ax = axes(fig);
histogram(ax, ipi, 0:30)
xlabel(ax, 'Inter-press interval (s)')
ylabel(ax, 'Num of trials')
