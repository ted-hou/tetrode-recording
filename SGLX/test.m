clear

% headWindow = [3000, 3020];
headWindow = [3600, 3660];

imec.binName = 'desmond38_20250407_g0_t0.imec0.ap.bin';
imec.path = '\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\desmond38\desmond38_20250407\desmond38_20250407\desmond38_20250407_g0';
imec.meta = SGLX_readMeta.ReadMeta(imec.binName, imec.path);
imec.fs = str2double(imec.meta.imSampRate);
imec.analogChannels = 1:384;

nidq.binName = 'desmond38_20250407_g0_t0.nidq.bin';
nidq.path = '\\research.files.med.harvard.edu\neurobio\Assad Lab\Lingfeng\Data\desmond38\desmond38_20250407\desmond38_20250407\desmond38_20250407_g0';
nidq.meta = SGLX_readMeta.ReadMeta(nidq.binName, nidq.path);
nidq.fs = str2double(nidq.meta.niSampRate);
nidq.analogChannels = 1:2;

imec.head.data = SGLX_readMeta.ReadBin(floor(imec.fs*headWindow(1)), floor(imec.fs*diff(headWindow)), imec.meta, imec.binName, imec.path);
imec.head.data = SGLX_readMeta.GainCorrectIM(imec.head.data, imec.analogChannels, imec.meta);
imec.head.analog = imec.head.data(imec.analogChannels, :);
imec.head.digital = SGLX_readMeta.ExtractDigital(imec.head.data, imec.meta, 1, 6);
% Perform common median referecing and bandpass filtering
imec.head.analog = bandpass(imec.head.analog', [250, 7500], imec.fs)';
imec.head.analog = imec.head.analog - median(imec.head.analog, 1);

nidq.head.data = SGLX_readMeta.ReadBin(floor(nidq.fs*headWindow(1)), floor(nidq.fs*diff(headWindow)), nidq.meta, nidq.binName, nidq.path);
nidq.head.data = SGLX_readMeta.GainCorrectNI(nidq.head.data, nidq.analogChannels, nidq.meta);
nidq.head.analog = nidq.head.data(nidq.analogChannels, :);
nidq.head.digital = SGLX_readMeta.ExtractDigital(nidq.head.data, nidq.meta, 1, 0:7);


imec.head.t = headWindow(1) + (0:1/imec.fs:floor(imec.fs*diff(headWindow)-1)/imec.fs);
nidq.head.t = headWindow(1) + (0:1/nidq.fs:floor(nidq.fs*diff(headWindow)-1)/nidq.fs);

%%
fig = figure();

% xl = [3009.2, 3009.4];
% xl = [3621, 3622];
xl = [3653.5, 3654.5]; 


ax = gobjects(2, 1);
ax(1) = subplot(2, 1, 1);
hold(ax(1), 'on')
plot(ax(1), imec.head.t, imec.head.analog(1:384, :)*10000 + (1:384)');
plot(ax(1), imec.head.t, imec.head.digital);
xlim(ax(1), xl)
ylim(ax(1), [0, 31])

ax(2) = subplot(2, 1, 2);
hold(ax(2), 'on')
plot(ax(2), nidq.head.t, nidq.head.analog + 8+[1; 2]);
plot(ax(2), nidq.head.t, nidq.head.digital + uint8((0:7))');
xlim(ax(2), xl)