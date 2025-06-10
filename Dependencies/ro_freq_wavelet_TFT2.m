function [pow,foi] = ro_freq_wavelet_TFT2(dat,cfg)

% Modification by JF to use linearly spaced wavelets (14.11.2024)

% This function is modified from Joerg Hipp's Wavelet transform code so that
% output is not averaged across time (gives a true time-frequency transform)
% Morlet wavelet transform of continous electrophysiological data
%
% Input
% - dat ... [channels x samples], set invalid data sections to NaN assuming the input signal was in uV
% - cfg ... struct, see script itself for default parameters, cfg.fsample
%           has to be set.
% Output
% 
% Pow = output power spectrum. The total power in uV^2 for channel ich over a 
%   frequency range [f1 f2] can be obtained with idf = foi> f1 & foi<f2 and
%       if cfg.norm is 'linear',   trapz(foi(idf),pow(ich,idf))
%       if cfg.norm is 'log2',   trapz(log2(foi(idf)),pow(ich,idf))

% 9.7.2016, Joerg Hipp
% 30.5.2017, Joerg Hipp, return foi_target instead of foi
% 02.02.2019, Joel Frohich, modified so that no averaging across time
% [later data], Joel Frohlich, modified so that frequency is linearly
% spaced [NOT log spaced]

% senity checks

dbstop if error

if ~isfield(cfg,'fsample'), error('specific cfg.fsample!'), end
n_sens=size(dat,1); n_sample=size(dat,2);
if n_sens>n_sample, fprintf('Warning, number of channles (%i) larger than number of samples (%i), check orientation of data matrix\n',n_sens,n_sample), end

if ~isfield(cfg,'verbose'),            cfg.verbose=1;                        end % if true, plot progress in command line
if ~isfield(cfg,'oct_bw'),             cfg.oct_bw = 0.5;                     end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'oct_delta'),          cfg.oct_delta = cfg.oct_bw/4;         end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'foi_start'),          cfg.foi_start = 2;                    end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'foi_end'),            cfg.foi_end = 32;                     end % wavelet frequency resolution -> q=f/sigma_f=
if ~isfield(cfg,'window_shift'),       cfg.window_shift=0.25;                end % fraction the window is shifted of total window size
if ~isfield(cfg,'kernel_width'),       cfg.kernel_width=5;                   end
if ~isfield(cfg,'allow_fraction_nan'), cfg.allow_fraction_nan=1;             end
if ~isfield(cfg,'norm');               cfg.norm='linear';                    end % 'linear', 'log2'
if ~isfield(cfg,'freq_shift_factor');  cfg.freq_shift_factor=1;              end

%% spectral parameter
foi        = cfg.foi_start:cfg.oct_delta:cfg.foi_end; % [Hz] frequencies to analyze
foi_target = foi;
foi        = foi*cfg.freq_shift_factor;
foi_tmp    = (cfg.foi_start)-(cfg.oct_delta/2):cfg.oct_delta:((cfg.foi_end)+cfg.oct_delta/2);
foi_tmp    = foi_tmp*cfg.freq_shift_factor;
foi_delta  = [foi_tmp(1:end-1);foi_tmp(2:end)]';          % [Hz] range for integration
foi_min    = 2*foi/(2^cfg.oct_bw+1);                      % arithmetic mean
foi_max    = 2*foi/(2^-cfg.oct_bw+1);
sigma_freq = (foi_max-foi_min)/(2*sqrt(2*log(2)));        % std in freq domain
sigma_time = 1./(2*pi*sigma_freq);                        % std in time domain

if cfg.verbose, tic, fprintf('Morelet Wavelet transform ['), end
for ifoi=1:length(foi)
    if cfg.verbose, fprintf('.'), end
    % convolution kernel
    %n_win = round(cfg.kernel_width*max(sigma_time)*cfg.fsample+1);
    n_win = cfg.fsample;
    n_shift = round(n_win*cfg.window_shift);
    t = ((1:n_win)-n_win/2-0.5)/cfg.fsample;
    z = t./sigma_time(ifoi);
    TAPER = exp(-(1/2)*z.^2);
    TAPER = TAPER/sqrt(sum(abs(TAPER).^2));
    iEXP = exp(1i*2*pi*foi(ifoi)*t);
    KERNEL = (TAPER.*iEXP).';
    
    % collect info on nan sections
    idx_up = find(diff([0,isnan(sum(dat))])==1);
    idx_down = find(diff([0,isnan(sum(dat))])==-1);
    nan_width = zeros(1,size(dat,2));
    for icnt=1:length(idx_up)-1
        nan_width(idx_up(icnt):idx_down(icnt)) = idx_down(icnt)-idx_up(icnt);
    end
    if length(idx_up)>length(idx_down)
        nan_width(idx_up(end):end) = length(nan_width)+1-idx_up(end);
    end
    
    % memory allocation
    DAT      = nan(n_sens,length(1:n_shift:n_sample-n_win+1));
    frac_nan = nan(1,size(DAT,2));
    % convolution
    cnt = 0;
    for isection = 1:n_shift:size(dat,2)-n_win+1
        section = double(dat(:,isection:isection+n_win-1));
        %nan_width_section = nan_width(isection:isection+n_win-1);
        n_nan = sum(isnan(section(1,:)));
        cnt=cnt+1;
        frac_nan(cnt) = n_nan/size(section,2);
        if n_nan==0
            DAT(:,cnt) = section*KERNEL*sqrt(2/cfg.fsample);
        %elseif n_nan<size(section,2)*cfg.allow_fraction_nan & ...
        %        max(nan_width_section)<size(section,2)*cfg.allow_fraction_nan
        else
            idx_valid = find(~isnan(section(1,:)));
            KERNEL_tmp = KERNEL(idx_valid)/sqrt(sum(abs(TAPER(idx_valid)).^2));
            DAT(:,cnt) = section(:,idx_valid)*KERNEL_tmp*sqrt(2/cfg.fsample);
        %else
        %    DAT(:,cnt) = nan(size(section,1),1);
        end % valid section
    end
    
    % derive measures for frequency-transformed data
    idx_valid = find(~isnan(DAT(1,:)));
    DAT = DAT(:,idx_valid);
    
    % get 2D time-frequency representation for all channels
    try
        pow(:,:,ifoi)        = abs(DAT).^2;                % [uV^2*Hz^-1] assuming the input signal was in uV
    catch
        pow(:,1:size(DAT,2),ifoi)        = abs(DAT).^2;
    end

end % loop foi
if cfg.verbose, fprintf('] %.1f sec\n',toc), end

% Normalization
if strcmp(cfg.norm,'log2'), error('this version of the function does not support log normalization'), end


