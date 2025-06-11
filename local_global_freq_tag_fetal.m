% Performs frequency tagging analysis on fetal MEG data for local-global 
% experiments, testing for neural responses at the stimulation frequency 
% (5/3 Hz) and a control frequency (3 Hz), while examining correlations 
% with demographic and physiological variables through bootstrapped 
% resampling.
%
% Joel Frohlich
% Last updated 11.06.2025

clearvars
load Data_traces.mat 
load Alltable.mat

close all
dbstop if error
try
    parpool
catch
    fprintf('Parpool already open\n')
end

fsample = 610.3516;
alpha_thresh = 0.05;
Nshift = 1000;
ntotal = 81; % total recordings across longitudinal visits

testfreq = 5/3*(round(fsample)/fsample); % do it this way because fft only takes integer sampling rates
control = 3*(round(fsample)/fsample); % do it this way because fft only takes integer sampling rates;  control frequency in Hz (from the other experiment)
T2 = readtable('Overview_Participants.xlsx');
T3 = readtable('dataset_differencesT3T4.csv');

%%% Common data %%%
condition = 'A';
norm = true;
common = true;

[s4A, d4A, HRVsA, HRVdA, AGEa, SEXa, REa, T0] = ...
    process_table(Datatable, T2, T3, condition, norm, common);

condition = 'B';

[s4B, d4B, HRVsB, HRVdB, AGEb, SEXb, REb] = ...
    process_table(Datatable, T2, T3, condition, norm, common);

% sanity check
assert(all(AGEa==AGEb))
assert(all(SEXa==SEXb))
assert(all(REa==REb))

AGE = AGEa;
SEX = SEXa;
RE = REa;

%%% Different data %%%

condition = 'A'; % just local standards
common = false;
% ssss and sssS
[s4Ax, d4Ax, HRVsAx, HRVdAx, AGEax, SEXax, REax,Tax] = ...
    process_table(Datatable, T2, T3, condition, norm, common);

condition = 'B'; % just local deviants
common = false;
% sssd and sssD
[s4Bx, d4Bx, HRVsBx, HRVdBx, AGEbx, SEXbx, REbx,Tbx] = ...
    process_table(Datatable, T2, T3, condition, norm, common);

T = [T0; Tax; Tbx];


AGE = [AGE; AGEax; AGEbx];
SEX = [SEX; SEXax; SEXbx];
RE = [RE; REax; REbx];
HRV = nanmean([ HRVsA HRVdA;  HRVsAx HRVdAx;  HRVsBx HRVdBx],2);
assert(all(isnan(setdiff(HRVsA,HRVdB))),'Problem with HRV, double check')
assert(all(isnan(setdiff(HRVdA,HRVsB))),'Problem with HRV, double check')

% frequency tagging

% % concatenate different conditions
m800 = round(0.8*fsample);
% NOTE - Don't z-score because the data are already normalized
alldat = [s4A(:,1:end-m800),d4A(:,1:end-m800),s4B(:,1:end-m800),d4B(:,1:end-m800); ... % [ssss sssD sssd sssS] data in both blocks
    s4Ax(:,1:end-m800),d4Bx(:,1:end-m800),s4Ax(:,1:end-m800),d4Bx(:,1:end-m800); ... % [ssss ssssD ssss sssD] A block
    s4Bx(:,1:end-m800),d4Ax(:,1:end-m800),s4Bx(:,1:end-m800),d4Ax(:,1:end-m800)]; % [sssd sssS sssd sssS] B block

assert(size(alldat,1) == ntotal,'wrong number of total recordings')

avgdat = mean(alldat,1);
numSub = size(alldat,1);

% stat test of amplitude at f Hz
Fs = round(fsample);
fres = 0.05; % frequency resolution
binperhz = 1/fres; 
N = Fs*binperhz;
foi = linspace(0,Fs/2,N/2);
tmp = abs(fft(avgdat,N));
amp_spectrum = tmp(1:N/2);

ratio = round(fsample)/fsample;

myfigure2
plot(foi,amp_spectrum,'linewidth',2)
xticks([1:8].*ratio)
xticklabels(1:8)
xlim([1/20 8*ratio])
xlabel('Frequency (Hz)')
ylabel('Amplitude (A.U.)')
title('Fetal MEG FFT amplitude spectrum, Experiment 2','fontsize',22)
plot(ones(1,100).*(testfreq),linspace(0,1e5,100),'k--')
legend({'Amplitude spectrum','Stimulation frequency'},'fontsize',16,...
    'location','northeast','box','off')
makefighandsome
xlim([0.5 8])
ylim([0 30000])
print('-dpng','fMEG_amplitude_spectrum_localglobal_fft.png')
print('-dsvg','fMEG_amplitude_spectrum_localglobal_fft.svg')

flank = round(10*(binperhz/20))+1;
[~,where] = min(abs(foi - testfreq));
left = amp_spectrum(where-flank:where-2);
right = amp_spectrum(where+2:where+flank);
try
    sigma = std([left right]);
    mu = mean([left right]);
catch
    sigma = std([left' right']);
    mu = mean([left' right']);
end
Z = (amp_spectrum(where) - mu)/sigma;
p_one_tailed = 1 - normcdf(Z);
Frq = unique(round(diff(foi),2));


% parametric
fprintf('Frequency tagging effect (%1.3f Hz freq. res.), Z = %1.3f, P = %1.3f\n',Frq,Z,p_one_tailed)

%% Control test

[~,where] = min(abs(foi - control));
left = amp_spectrum(where-flank:where-2);
right = amp_spectrum(where+2:where+flank);
try
    sigma = std([left right]);
    mu = mean([left right]);
catch
    sigma = std([left' right']);
    mu = mean([left' right']);
end
xZ = (amp_spectrum(where) - mu)/sigma;
xp_one_tailed = 1 - normcdf(xZ);
xFrq = unique(round(diff(foi),2));

% parametric
fprintf('CONTROL: Frequency tagging effect (%1.3f Hz freq. res.), Z = %1.3f, P = %1.3f\n',xFrq,xZ,xp_one_tailed)



%% plot reponse
myfigure2
time = linspace(0,9.6,length(avgdat));
plot(time,avgdat,'linewidth',2)
beep = 0.2:0.6:9.2;
for ibeep = 1:length(beep)
    plot(repmat(beep(ibeep),1,100),linspace(25,150,100),'k:')
    if ibeep == 1 
        legend({'MEG signal','Stimulus onset'},'autoupdate','off',...
            'fontsize',18,'location','northeast')
        legend box off
    end
end

xlabel('time (sec)')
ylabel('Field strength (AU)')
ylim([25 150])
title('Grand-averaged fetal MEG signal, Experiment 2','fontsize',22)
makefighandsome

print('-dsvg','Neural_response_tones_LG_fetus.svg')

%% Correlational analysis

% extract fetal sex data
Tvar = Alltable(~Alltable.BornYet,:);
whichones =[];
for isub = 1:length(T.ID)
    whichones(isub) = find(T.ID(isub) == Tvar.ID,1,'first');
end

assert(all(SEX == Tvar.XChrom(whichones)),'wrong sex data');
MAGE = Tvar.MomsAge(whichones);
WGHT = Tvar.BirthWeight(whichones);
BMI = Tvar.BMI(whichones);

nsub = size(alldat,1);
nboot = 1000;
nrs = 500;

ages = nan(nboot, nrs);
sexS = nan(nboot, nrs);
mageS = nan(nboot, nrs);
wghtS = nan(nboot, nrs);
bmiS = nan(nboot, nrs);
snr = nan(nboot, nrs);
zz = nan(nboot, nrs);

ctype = 'pearson';

parfor iboot = 1:nboot
    rng(iboot,"twister")
    for irs = 1:nrs
        NDX = randi(nsub,1,nsub);
        tmp = abs(fft(mean(alldat(NDX,:)),N));
        ages(iboot,irs) = mean(AGE(NDX));
        hrvS(iboot,irs) = nanmean(HRV(NDX));
        sexS(iboot,irs) = mean(SEX(NDX))-1;
        mageS(iboot,irs) = mean(MAGE(NDX));
        wghtS(iboot,irs) = mean(WGHT(NDX));
        bmiS(iboot,irs) = mean(BMI(NDX));

        amp_spectrum = tmp(1:N/2);
        [~,where] = min(abs(foi - testfreq));
        left = amp_spectrum(where-11:where-2);
        right = amp_spectrum(where+2:where+11);
        snr(iboot,irs) = amp_spectrum(where)/mean([left right]);
        sigma = std([left right]);
        mu = mean([left right]);
        zz(iboot,irs) = (amp_spectrum(where) - mu)/sigma;
    end
    fprintf('Bootstrapping is %1.2f%% complete\n',iboot/nboot*100)
end

% calculate each correlation
Rhrv = diag(corr(hrvS',zz','rows','complete','type',ctype));
Rage = diag(corr(ages',zz','rows','complete','type',ctype));
Rsex = diag(corr(sexS',zz','rows','complete','type',ctype));
Rmage = diag(corr(mageS',zz','rows','complete','type',ctype));
Rwght = diag(corr(wghtS',zz','rows','complete','type',ctype));
Rbmi = diag(corr(bmiS',zz','rows','complete','type',ctype));

%% Compute parametric P-values

[Page, Rm_age, Rci_age] = getP(Rage,nrs);
[Phrv, Rm_hrv, Rci_hrv] = getP(Rhrv,nrs);
[Psex, Rm_sex, Rci_sex] = getP(Rsex,nrs);
[Pmage,Rm_mage, Rci_mage] = getP(Rmage,nrs);
[Pwght,Rm_wght, Rci_wght] = getP(Rwght,nrs);
[Pbmi, Rm_bmi, Rci_bmi] = getP(Rbmi,nrs);

fprintf('Experimental correlation with age, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    Rm_age,Rci_age(1),Rci_age(2),Page)
fprintf('Experimental correlation with HRV, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    Rm_hrv,Rci_hrv(1),Rci_hrv(2),Phrv)
fprintf('Experimental correlation with sex, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    Rm_sex,Rci_sex(1),Rci_sex(2),Psex)
fprintf('Experimental correlation with mother''s age, r = %1.2f (CI: r = %1.2f - %1.2f),P = %1.3f (two-tailed)\n',...
    Rm_mage,Rci_mage(1),Rci_mage(2),Pmage)
fprintf('Experimental correlation with birth weight, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    Rm_wght,Rci_wght(1),Rci_wght(2),Pwght)
fprintf('Experimental correlation with BMI, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    Rm_bmi,Rci_bmi(1),Rci_bmi(2),Pbmi)


% Finally, repeat the correlational analysis using a control frequency

[~,cf] = min(abs(foi - control));

agex = nan(nboot, nrs);
sexx = nan(nboot, nrs);
magex = nan(nboot, nrs);
wghtx = nan(nboot, nrs);
bmix = nan(nboot, nrs);
snrx = nan(nboot, nrs);
zzx = nan(nboot, nrs);

ctype = 'pearson';

parfor iboot = 1:nboot
    rng(iboot,"twister")
    for irs = 1:nrs
        NDX = randi(nsub,1,nsub);
        tmp = abs(fft(mean(alldat(NDX,:)),N));
        agex(iboot,irs) = mean(AGE(NDX));
        hrvx(iboot,irs) = nanmean(HRV(NDX));
        sexx(iboot,irs) = mean(SEX(NDX))-1;
        magex(iboot,irs) = mean(MAGE(NDX));
        wghtx(iboot,irs) = mean(WGHT(NDX));
        bmix(iboot,irs) = mean(BMI(NDX));

        amp_spectrum = tmp(1:N/2);
        [~,cf] = min(abs(foi - control));
        left = amp_spectrum(cf-11:cf-2);
        right = amp_spectrum(cf+2:cf+11);
        snrx(iboot,irs) = amp_spectrum(cf)/mean([left right]);
        sigma = std([left right]);
        mu = mean([left right]);
        zzx(iboot,irs) = (amp_spectrum(cf) - mu)/sigma;
    end
    fprintf('Bootstrapping is %1.2f%% complete\n',iboot/nboot*100)
end

% calculate each correlation
xRhrv = diag(corr(hrvS',zzx','rows','complete','type',ctype));
xRage = diag(corr(ages',zzx','rows','complete','type',ctype));
xRsex = diag(corr(sexS',zzx','rows','complete','type',ctype));
xRmage = diag(corr(mageS',zzx','rows','complete','type',ctype));
xRwght = diag(corr(wghtS',zzx','rows','complete','type',ctype));
xRbmi = diag(corr(bmiS',zzx','rows','complete','type',ctype));

%% Compute parametric P-values

[xPage, xRm_age, xRci_age] = getP(xRage,nrs);
[xPhrv, xRm_hrv, xRci_hrv] = getP(xRhrv,nrs);
[xPsex, xRm_sex, xRci_sex] = getP(xRsex,nrs);
[xPmage,xRm_mage, xRci_mage] = getP(xRmage,nrs);
[xPwght,xRm_wght, xRci_wght] = getP(xRwght,nrs);
[xPbmi, xRm_bmi, xRci_bmi] = getP(xRbmi,nrs);

fprintf('CONTROL correlation with age, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    xRm_age,xRci_age(1),xRci_age(2),xPage)
fprintf('CONTROL correlation with HRV, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    xRm_hrv,xRci_hrv(1),xRci_hrv(2),xPhrv)
fprintf('CONTROL correlation with sex, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    xRm_sex,xRci_sex(1),xRci_sex(2),xPsex)
fprintf('CONTROL correlation with mother''s age, r = %1.2f, (CI: r = %1.2f - %1.2f),  P = %1.3f (two-tailed)\n',...
    xRm_mage,xRci_mage(1),xRci_mage(2),xPmage)
fprintf('CONTROL correlation with birth weight, r = %1.2f, (CI: r = %1.2f - %1.2f) ,P = %1.3f (two-tailed)\n',...
    xRm_wght,xRci_wght(1),xRci_wght(2),xPwght)
fprintf('CONTROL correlation with BMI, r = %1.2f (CI: r = %1.2f - %1.2f), P = %1.3f (two-tailed)\n',...
    xRm_bmi,xRci_bmi(1),xRci_bmi(2),xPbmi)


%% frequency transform of surrogates

[~,where] = min(abs(foi - testfreq));
shiftedData = nan(size(alldat,1),size(alldat,2),Nshift);
SNRshift = nan(1,Nshift);
for isub = 1:numSub
    rng(isub, 'Twister'); % Set the random seed for reproducibility
    shiftedData(isub,:,:) = surrogate(alldat(isub,:), Nshift, 'IAAFT1', 0, fsample)';
    fprintf('%1.2f%% complete\n',isub/numSub*100)
end

parfor ishift = 1:Nshift
    x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
    left = x_spectrum(where-flank:where-2);
    right = x_spectrum(where+2:where+flank);
    sigma = std([left' right']);
    mu = mean([left' right']);
    SNRshift(ishift) = (x_spectrum(where) - mu)/sigma;
    fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
end


% non-parametric
Pval = sum(SNRshift >= Z)/Nshift;
fprintf('Frequency tagging effect, non-parametric test, P = %1.2f\n',Pval)

%% control test
[~,where] = min(abs(foi - control));
xSNRshift = nan(1,Nshift);
parfor ishift = 1:Nshift
    x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
    left = x_spectrum(where-flank:where-2);
    right = x_spectrum(where+2:where+flank);
    sigma = std([left' right']);
    mu = mean([left' right']);
    xSNRshift(ishift) = (x_spectrum(where) - mu)/sigma;
    fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
end

% non-parametric
xPval = sum(xSNRshift >= xZ)/Nshift;
fprintf('CONTROL: Frequency tagging effect, non-parametric test, P = %1.2f\n',xPval)

%% Output table

P = [p_one_tailed Pval Page Pwght Pbmi Pmage Psex Phrv];
stat = [Z nan Rm_age Rm_wght Rm_bmi Rm_mage Rm_sex Rm_hrv];
cilo = [nan nan Rci_age(1) Rci_wght(1) Rci_bmi(1) Rci_mage(1) Rci_sex(1) Rci_hrv(1)];
cihi = [nan nan Rci_age(2) Rci_wght(2) Rci_bmi(2) Rci_mage(2) Rci_sex(2) Rci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tout = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tout,'LG_fetal_stats_out.csv')

% Control test stats

P = [xp_one_tailed xPval xPage xPwght xPbmi xPmage xPsex xPhrv];
stat = [xZ nan xRm_age xRm_wght xRm_bmi xRm_mage xRm_sex xRm_hrv];
cilo = [nan nan xRci_age(1) xRci_wght(1) xRci_bmi(1) xRci_mage(1) xRci_sex(1) xRci_hrv(1)];
cihi = [nan nan xRci_age(2) xRci_wght(2) xRci_bmi(2) xRci_mage(2) xRci_sex(2) xRci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tctr = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tctr,'LG_fetal_CONTROL_stats_out.csv')


%% save everything
close all
save(sprintf('local_global_freqtag_fetal_datamax_%s',date),'-v7.3')



%% plot with smoothing 

% Wavelet parameters
cfg = [];
cfg.fsample = Fs;
cfg.foi_start = 0.5;
cfg.foi_end = 8;
cfg.window_shift=0.5;
cfg.oct_bw = 0.2;
cfg.oct_delta = cfg.oct_bw/4;
cfg.verbose = false;


% % FOOOF
% settings=[];
% settings.peak_width_limits = [0.5 2];
% settings.max_n_peaks = 2;
% settings.aperiodic_mode = 'fixed';
% fa = 0.5;
% fb = 10;
% f_range = [fa fb];

% Note: we MUST zscore (or at least demean) the signal, otherwise the
% wavelet transform breaks when the input signal isn't oscillating about
% 0 ... the raw signal is in units percent change and never crosses 0.
[pow_all,foi2] = ro_freq_wavelet_TFT2(zscore(avgdat),cfg);
%%fooof_results = fooof(foi', squeeze(mean(pow_all,2)), f_range, settings,true);

myfigure2
plot(foi2,log10(squeeze(mean(pow_all,2))),'linewidth',2)
ylim([-2 0])
xlabel('Frequency (Hz)')
ylabel('Amplitude (A.U.)')
title('Fetal MEG power spectrum, Experiment 2','fontsize',22)
plot(ones(1,100).*testfreq,linspace(min(log10(pow)),max(log10(pow)),100),'k--')
makefighandsome
xticks([1:8].*ratio)
xticklabels(1:8)
xlim([0.5*ratio 8*ratio])
xlabel('Frequency (Hz)')
print('-dpng','fMEG_power_spectrum_localglobal_smooth')
print('-dsvg','fMEG_power_spectrum_localglobal_smooth')


%% helper functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[P,Rmed,Rci] = getP(R,n)

Rmed = median(R); 

% get the parametric p-value
% Calculate the t-statistic
t = (Rmed * sqrt(n - 2)) / sqrt(1 - Rmed^2);

% Degrees of freedom
df = n - 2;

% Calculate the p-value for a two-tailed test
P = 2 * (1 - tcdf(abs(t), df));

Rci = prctile(R,[2.5 97.5]);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[standard, deviant, HRVs, HRVd, AGE, SEX, RE,T] = ...
    process_table(Datatable, T2, T3, condition, norm, common)


assert(islogical(common))
assert(islogical(norm))

switch condition
    case 'A'
        idx1 = all(Datatable.Blocklabel == 'A (ssss)',2); % ssss trials
        idx2 = all(Datatable.Blocklabel == 'B (sssd)',2); % sssS trials
    case 'B'
        idx1 = all(Datatable.Blocklabel == 'B (sssd)',2); % sssd trials
        idx2 = all(Datatable.Blocklabel == 'A (ssss)',2); % sssD trials
    otherwise 
        error('Unrecognized condition')
end

T = Datatable(idx1,:);
Tdev = Datatable(idx2,:);

tmp1 = T.ID*100 + T.GA;
T.tmpID = tmp1;
tmp2 = Tdev.ID*100 + Tdev.GA;
Tdev.tmpID = tmp2;

switch common % are we working on the common data?
    case true
        [~,Ia,Ib] = intersect(tmp1,tmp2); %find common data

        T = T(Ia,:);
        Tdev = Tdev(Ib,:);
        assert(all(T.tmpID==Tdev.tmpID),'Data don''t match between blocks')

    case false
        [~,Xb] = setdiff(tmp2,tmp1);
        [~,Xa] = setdiff(tmp1,tmp2);
        T = T(Xa,:);
        Tdev = Tdev(Xb,:);
        assert(isempty(intersect(T.tmpID,Tdev.tmpID)),'These are not the data that are missing in one condition')
end


standard = [];
deviant = [];
HRVa = nan(size(T,1),1);
HRVb = nan(size(T,1),1);
AGE = nan(size(T,1),1);
RE = nan(size(T,1),1);
SEX = nan(size(T,1),1);


for irow = 1:size(T,1)
    switch norm
        case true
            idx = T.global_standards_normalized(irow)';
        case false
            idx = T.global_standards(irow)';
    end
    standard = [standard; idx{:}];
    AGE(irow) = T.GA(irow);
    RE(irow) = T.ID(irow);

    % Add sex as a covariate
    for jrow = 1:size(T2,1)
        if T2.ID(jrow) == T.ID(irow)
            switch T2.FetalSex{jrow}
                case 'm'
                    SEX(irow) = 1;
                    break
                case 'f'
                    SEX(irow) = 2;
                    break
                otherwise
                    continue
            end
        end
    end

    % Add SDNN
    for jrow = 1:size(T3,1)
        if T3.ID(jrow) == T.ID(irow) && T3.GA(jrow) == T.GA(irow) ...
                && contains(T3.Block{jrow},'A')
            HRVa(irow) = T3.SDNN(jrow);
        elseif T3.ID(jrow) == T.ID(irow) && T3.GA(jrow) == T.GA(irow) ...
                && contains(T3.Block{jrow},'B')
            HRVb(irow) = T3.SDNN(jrow);
        end
    end
    if mod(irow,10) == 0, fprintf('STEP 1 %1.1f%% complete\n',irow/size(T,1)*100), end
end

clear idx

for irow = 1:size(Tdev,1)
    switch norm
        case true
            idx = Tdev.global_deviants_normalized(irow)';
        case false
            idx = Tdev.global_deviants(irow)';
    end
    deviant = [deviant; idx{:}];
end

switch condition
    case 'A'
        HRVs = HRVa;
        HRVd = HRVb;
    case 'B'
        HRVs = HRVb;
        HRVd = HRVa;
end

end

