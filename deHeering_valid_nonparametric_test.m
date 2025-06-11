% Validates the non-parametric frequency tagging analysis method by 
% comparing it against the parametric approach from de Heering and Rossion 
% 2015, using EEG data to test for significant responses at 1.2 Hz across 
% 32 channels and generating topographic plots and scatter plots to 
% visualize the comparison.
% 
% Joel Frohlich
% Last updated 10.06.2025


files = dir('./deHeering_Rossion_EEGData/Experiment_1/**/*.mat');
nch = 32;
fs = 250;
Nsur = 200;
order = nan(32,15);
Nshift = 1000;
numSub = length(files);
testfreq = 1.2;

fres = 0.05; % frequency resolution
binperhz = 1/fres;
N = fs*binperhz;

% Note: channel order has bee n confirmed in the .lw5 files

% Design the Butterworth filter
Fc1 = 0.5; % Lower cutoff frequency
Fc2 = 20;  % Upper cutoff frequency
Wn = [2*Fc1/fs, 2*Fc2/fs]; % Normalized cutoff frequencies for bandpass
[b,a] = butter(5, Wn, 'bandpass');

foi = linspace(0,round(fs)/2,N/2);
[~,where] = min(abs(foi - testfreq));
flank = round(10*(binperhz/20))+1;

for ich = 1:nch
    alldat = nan(numSub,5000);
    for ifile = 1:length(files)
        eeg = load(sprintf('%s/%s',files(ifile).folder,files(ifile).name));
        eegdat = squeeze(eeg.data);
        alldat(ifile,:) = zscore(filtfilt(b,a,eegdat(ich,:)));
    end

    tmp = abs(fft(mean(alldat),N));
    amp_spectrum = tmp(1:N/2);
    left = amp_spectrum(where-flank:where-2);
    right = amp_spectrum(where+2:where+flank);
    try
        SNR = amp_spectrum(where)/mean([left right]);
        sigma = std([left right]);
        mu = mean([left right]);
    catch
        SNR = amp_spectrum(where)/mean([left' right']);
        sigma = std([left' right']);
        mu = mean([left' right']);
    end
    Z = (amp_spectrum(where) - mu)/sigma;


    %% frequency transform of surrogates

    shiftedData = nan(size(alldat,1),size(alldat,2),Nshift);
    parfor isub = 1:numSub
        rng(isub, 'Twister'); % Set the random seed for reproducibility
        shiftedData(isub,:,:) = surrogate(alldat(isub,:), Nshift, 'IAAFT1', 0, fs)';
        fprintf('Channel %i out of %i, %1.2f%% complete\n',ich, nch, isub/numSub*100)
    end

    parfor ishift = 1:Nshift
        x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
        left = x_spectrum(where-flank:where-2);
        right = x_spectrum(where+2:where+flank);
        SNRshift(ishift) = x_spectrum(where)/mean([left right]);
        fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
    end

    % non-parametric
    Pval(ich) = sum(SNRshift > SNR)/Nshift;
    fprintf('Frequency tagging effect, non-parametric test, P = %1.2f\n',Pval(ich))

end

save deHeering_validation 

%% Compare results with de Heering et al. 2015

channel_labels = {'Fp1', 'AF3', 'F7', 'F3', 'FC1', 'FC5', 'T7', 'C3', ...
                  'CP1', 'CP5', 'P7', 'P3', 'Pz', 'PO3', 'O1', 'Oz', ...
                  'O2', 'PO4', 'P4', 'P8', 'CP6', 'CP2', 'C4', 'T8', ...
                  'FC6', 'FC2', 'F4', 'F8', 'AF4', 'Fp2', 'Fz', 'Cz'};

T0 = readtable('./deHeering_Rossion_EEGData/deHeering_results.csv');
order = sortref(T0.Electrode,channel_labels');
T = T0(order,:);
T.Pval = 1 - normcdf(T.Zscore);
flr = min(T.Pval(T.Pval>0));
T.Pval(T.Pval==0) = flr;
Pval(Pval==0) = 1/(Nshift+1);

P_heer = -log10(T.Pval); % Frohlich et al. approach
P_froh = -log10(Pval); % de Heering and Rossion approach

myfigure2

x = [-1 -1 10 10];
y = [1.3 10 10 1.3];

% Plot the rectangle
fill(x, y, 'y', 'facealpha',0.3,'edgecolor','none'); % 'r' specifies the color red

x = [1.3 1.3 10 10];
y = [-1 10 10 -1];

% Plot the rectangle
fill(x, y, 'c', 'facealpha',0.3,'edgecolor','none'); % 'r' specifies the color red

scatter(P_froh,P_heer,100,'linewidth',2)
mylsline

legend({'P < 0.05 (Parametric','P < 0.05 (Non-parametric)','Channels','Least squares fit'},...
    'AutoUpdate','off','Box','Off','fontsize',17,'location','northeastoutside')

ylabel('-log10(P) Parametric (de Heering and Rossion)')
xlabel('-log10(P) Non-parametric (Frohlich et al.)')
[r,p] = corr(-log10(Pval)',-log10(T.Pval));
title(sprintf('r = %1.2f, p = %1.1e',r,p),'fontsize',18)

text(P_froh,P_heer,T.Electrode,'fontsize',12,...
    'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
axis([-0.2 3.25 -0.2 5])
makefigpretty

print('-dpng','Nonparametric_validation_scatter_plot.png')
print('-dsvg','Nonparametric_validation_scatter_plot.svg')

% show channels that are significant according to both approaches

idx = find(P_froh' >= -log10(0.05) & P_heer >= -log10(0.05));

fprintf('The following channels are sig. with both approaches: ')
for ich = 1:length(idx)
    fprintf('%s ',string(channel_labels{idx(ich)}))
end
fprintf('\n')

idx = find(P_froh' >= -log10(0.05) & ~(P_heer >= -log10(0.05)));

fprintf('The following channels are sig. ONLY with the nonparametric method: ')
for ich = 1:length(idx)
    fprintf('%s ',string(channel_labels{idx(ich)}))
end
fprintf('\n')

idx = find(~(P_froh' >= -log10(0.05)) & P_heer >= -log10(0.05));

fprintf('The following channels are sig. ONLY with the parametric method: ')
for ich = 1:length(idx)
    fprintf('%s ',string(channel_labels{idx(ich)}))
end
fprintf('\n')


%% Topoplots

% Ensure FieldTrip is properly initialized
if ~exist('ft_defaults', 'file')
    error('FieldTrip is not properly initialized. Please run ft_defaults first.');
end

% Load the layout file explicitly first
try
    cfg_layout = [];
    cfg_layout.layout = 'biosemi32.lay';
    layout = ft_prepare_layout(cfg_layout);
catch layout_error
    error('Failed to load layout file. Error: %s', layout_error.message);
end

% Create the data structure with careful initialization
dummy = [];
dummy.time = 0;                   % Single time point
dummy.label = layout.label(1:32); % Use labels from layout to ensure match
dummy.avg = zeros(32, 1);         % Initialize with zeros first
dummy.dimord = 'chan_time';       % Required dimension order

% Assign the p-values
if ~exist('P_froh', 'var')
    error('P_froh variable not found in workspace');
end

if numel(P_froh) ~= 32
    error('P_froh must contain exactly 32 values (one per channel). Current size: %d', numel(P_froh));
end

% Assign p-values, ensuring correct orientation
dummy.avg = P_froh';

% Create minimal configuration to reduce potential conflicts
cfg = [];
cfg.layout = layout;              % Use the already loaded layout
cfg.parameter = 'avg';
cfg.style = 'fill';
cfg.marker = 'on';
cfg.markersize = 8;
cfg.comment = 'no';
cfg.gridscale = 300; 
rwb = customcolormap_preset('red-white-blue');
rwb = rwb(size(rwb,1)/2+1:end,:);
cfg.colormap = rwb;  
cfg.markersymbol  = 'o';
cfg.markersize = 18;
cfg.marker = 'labels';
cfg.markercolor = 'y';
cfg.numcontour = 6; 

CMAX = 2;

% Create figure
myfigure
try
    ft_topoplotER_JF(cfg, dummy);
    mycolorbar
    clim([0 CMAX])
    title('Non-parametric P-values','fontsize',18);
catch
    % something bounces here after it's already plotted, just continue with
    % the plot labels
    mycolorbar
    clim([0 CMAX])
    title('Non-parametric P-values','fontsize',18);
end
c.Label.String = '-log10(P)';
c.Label.FontSize = 24;
print('-dpng','nonparametric_validation_topo.png')
print('-dsvg','nonparametric_validation_topo.svg')

dummy.avg = P_heer';
myfigure
try
    ft_topoplotER_JF(cfg, dummy);
    mycolorbar
    clim([0 CMAX])
    title('Parametric P-values','fontsize',18);
catch
    % something bounces here after it's already plotted, just continue with
    % the plot labels
    mycolorbar
    clim([0 CMAX])
    title('Parametric P-values','fontsize',18);
end
c.Label.String = '-log10(P)';
c.Label.FontSize = 24;
print('-dpng','parametric_validation_topo.png')
print('-dsvg','pparametric_validation_topo.svg')





