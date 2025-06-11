% Performs frequency tagging analysis on fetal MEG data from statistical
% learning experiments, testing for neural responses at 3 Hz versus a 
% control frequency (5/3 Hz), while examining correlations with 
% demographic and physiological variables through bootstrapped resampling 
% and combining data from R (random) and S (structured) experimental 
% conditions.
%
% Joel Frohlich
% 11.06.2025

clearvars
r = load('stat_learning_output_RMS_27-May-2025.mat');
s = load('stat_learning_output_RMS_embedded_pattern_27-May-2025.mat');

hrv = load('stat_learning_hrv_26-Mar-2025.mat');

r.hrv = nan(1,length(r.SID));
s.hrv = nan(1,length(s.SID));

% get HRV as SDNN
% R condition (random tones)
for i = 1:length(r.SID)
    id = sprintf('SL%02d-R', r.SID(i));
    for j = 1:length(hrv.meg_recording)
        if ~isempty(hrv.meg_recording{j}) && contains(hrv.meg_recording{j},id)
            r.hrv(i) = hrv.hrv(j);
        end
    end
end
% S condition (structured tone patterns)
for i = 1:length(s.SID)
    id = sprintf('SL%02d-S', s.SID(i));
    for j = 1:length(hrv.meg_recording)
        if ~isempty(hrv.meg_recording{j}) && contains(hrv.meg_recording{j},id)
            s.hrv(i) = hrv.hrv(j);
        end
    end
end


%%
tmp = load('stat_learning_demo_vars');
fsample = 610.3516;
Td = tmp.Td; clear tmp

for irow = 1:size(Td,1)
    for jrow = 1:size(r.SID,1)
        if Td.ID(irow) == r.SID(jrow)
            r.sex(jrow) = Td.sex(irow);
            r.bmi(jrow) = Td.BMI(irow);
            r.wght(jrow) = Td.birth_weight(irow);
            r.bage(jrow) = Td.birth_GA(irow);
            r.mage(jrow) = Td.maternal_age(irow);
        end
    end

    for jrow = 1:size(s.SID,1)
        if Td.ID(irow) == s.SID(jrow)
            s.sex(jrow) = Td.sex(irow);
            s.bmi(jrow) = Td.BMI(irow);
            s.wght(jrow) = Td.birth_weight(irow);
            s.bage(jrow) = Td.birth_GA(irow);
            s.mage(jrow) = Td.maternal_age(irow);
        end
    end
end

% identify subjects in each condition 
whoR = r.SID;
whoS = s.SID;

[~,where1,where2] = intersect(whoS,whoR);
[~,there1] = setdiff(whoS,whoR);
[~,there2] = setdiff(whoR,whoS); 

assert(all(s.ages(where1) == r.ages(where2)),'data aren''t matched')
assert(all(s.sex(where1) == r.sex(where2)),'data aren''t matched')

AGE = [s.ages(where1)' s.ages(there1)' r.ages(there2)'];
SEX = [s.sex(where1) s.sex(there1) r.sex(there2)];
BMI = [s.bmi(where1) s.bmi(there1) r.bmi(there2)];
BAGE = [s.bage(where1) s.bage(there1) r.bage(there2)];
MAGE = [s.mage(where1) s.mage(there1) r.mage(there2)];
WGHT = [s.wght(where1) s.wght(there1) r.wght(there2)];
HRV = [s.hrv(where1) s.hrv(there1) r.hrv(there2)];
SID = [s.SID(where1)' s.SID(there1)' r.SID(there2)'];

% if a subject doesn't have usable data from one condition, then make the
% most of data from the usable condition
alldat = [s.avgdat1(where1,:) r.avgdat1(where2,:); ...
    s.avgdat1(there1,:) s.avgdat1(there1,:); ...
    r.avgdat1(there2,:) r.avgdat1(there2,:)];

% check that we have the correct number of recordings, one per subject
assert(size(alldat,1) == length(unique(SID)),'Wrong number of MEG recordings')

% Normalize signal by zscore
alldat = zscore(alldat')';

% number of iterations for non-parametric testing
Nshift = 1000;

testfreq = 3*(round(fsample)/fsample); % do it this way because fft only takes integer sampling rates; 
control = 5/3*(round(fsample)/fsample); % do it this way because fft only takes integer sampling rates

% sampling rate
Fs0 = 610.3516;
Fs = 610; % down-sample target
numSub = size(alldat,1);

%% stat test of amplitude at f Hz
N = round(Fs)*20;
foi = linspace(0,round(Fs)/2,N/2);
tmp = abs(fft(mean(alldat),N));
amp_spectrum = tmp(1:N/2);
ratio = round(Fs)/Fs;

myfigure2
plot(foi,amp_spectrum,'linewidth',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude (A.U.)')
title('MEG FFT amplitude spectrum, Experiment 1','fontsize',22)
plot(ones(1,100).*(testfreq),linspace(0,3000,100),'k--')
makefighandsome
xlim([0.5 8*ratio])
ylim([0 150])
print('-dpng','fMEG_amplitude_spectrum__fft.png')
print('-dsvg','fMEG_amplitude_spectrum__fft.svg')


[~,where] = min(abs(foi - testfreq));
left = amp_spectrum(where-11:where-2);
right = amp_spectrum(where+2:where+11);
sigma = std([left right]);
mu = mean([left right]);
Z = (amp_spectrum(where) - mu)/sigma;
p_one_tailed = 1 - normcdf(Z);
Frq = unique(round(diff(foi),2));

% parametric
fprintf('Frequency tagging effect (%1.3f Hz freq. res.), Z = %1.3f, P = %1.3f\n',Frq,Z,p_one_tailed)

%% Control test

[~,where] = min(abs(foi - control));
left = amp_spectrum(where-11:where-2);
right = amp_spectrum(where+2:where+11);
sigma = std([left right]);
mu = mean([left right]);
xZ = (amp_spectrum(where) - mu)/sigma;
xp_one_tailed = 1 - normcdf(xZ);
xFrq = unique(round(diff(foi),2));

% parametric
fprintf('CONTROL: Frequency tagging effect (%1.3f Hz freq. res.), Z = %1.3f, P = %1.3f\n',xFrq,xZ,xp_one_tailed)

%% plot reponse
myfigure2
time = linspace(0,8,length(alldat)); 
plot(time,mean(alldat),'linewidth',2)
beep = linspace(0,8,24); 
for ibeep = 1:length(beep)
    plot(repmat(beep(ibeep),1,100),linspace(-0.8, 0.8,100),'k:')
    if ibeep == 1 
        legend({'MEG signal','Stimulus onset'},'autoupdate','off',...
            'fontsize',18,'location','northeast')
        legend box off
    end
end


xlabel('time (sec)')
ylabel('Field strength (AU)')
ylim([-0.601 0.601])
makefighandsome

print('-dsvg','Neural_response_tones_SL_fetus_RMS.svg')

%% Bootstrapped correlations

nsub = size(alldat,1);
nboot = 1000;
nrs = 500;

ages = nan(nboot, nrs);
sexS = nan(nboot, nrs);
mageS = nan(nboot, nrs);
wghtS = nan(nboot, nrs);
bmiS = nan(nboot, nrs);
hrvS = nan(nboot, nrs);
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

%% calculate each correlation
Rhrv = diag(corr(hrvS',zz','rows','complete','type',ctype));
Rage = diag(corr(ages',zz','rows','complete','type',ctype));
Rsex = diag(corr(sexS',zz','rows','complete','type',ctype));
Rmage = diag(corr(mageS',zz','rows','complete','type',ctype));
Rwght = diag(corr(wghtS',zz','rows','complete','type',ctype));
Rbmi = diag(corr(bmiS',zz','rows','complete','type',ctype));

% Compute parametric P-values

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

%% Finally, repeat the correlational analysis using a control frequency

[~,cf] = min(abs(foi - control));

agex = nan(nboot, nrs);
sexx = nan(nboot, nrs);
magex = nan(nboot, nrs);
wghtx = nan(nboot, nrs);
bmix = nan(nboot, nrs);
hrvx = nan(nboot, nrs);
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

%% calculate each correlation

xRhrv = diag(corr(hrvx',zzx','rows','complete','type',ctype));
xRage = diag(corr(agex',zzx','rows','complete','type',ctype));
xRsex = diag(corr(sexx',zzx','rows','complete','type',ctype));
xRmage = diag(corr(magex',zzx','rows','complete','type',ctype));
xRwght = diag(corr(wghtx',zzx','rows','complete','type',ctype));
xRbmi = diag(corr(bmix',zzx','rows','complete','type',ctype));

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

shiftedData = nan(size(alldat,1),size(alldat,2),Nshift);
parfor isub = 1:numSub
    rng(isub, 'Twister'); % Set the random seed for reproducibility
    shiftedData(isub,:,:) = surrogate(alldat(isub,:), Nshift, 'IAAFT1', 0, Fs)';
    fprintf('%1.2f%% complete\n',isub/numSub*100)
end


[~,where] = min(abs(foi - testfreq));
SNRshift = nan(1,Nshift);
parfor ishift = 1:Nshift
    x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
    left = x_spectrum(where-11:where-2);
    right = x_spectrum(where+2:where+11);
    sigma = std([left right]);
    mu = mean([left right]);
    SNRshift(ishift) = (x_spectrum(where) - mu)/sigma;
    fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
end


% non-parametric
Pval = sum(SNRshift >= Z)/Nshift;
fprintf('Frequency tagging effect, non-parametric test, P = %1.2f\n',Pval)


%% Control test
[~,where] = min(abs(foi - control));
xSNRshift = nan(1,Nshift);
parfor ishift = 1:Nshift
    x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
    left = x_spectrum(where-11:where-2);
    right = x_spectrum(where+2:where+11);
    sigma = std([left right]);
    mu = mean([left right]);
    xSNRshift(ishift) = (x_spectrum(where) - mu)/sigma;
    fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
end

% non-parametric
xPval = sum(xSNRshift >= xZ)/Nshift;
fprintf('CONTROL: Frequency tagging effect, non-parametric test, P = %1.2f\n',xPval)

%% Output table

if ~exist('Phrv','var')
    load('stat_learn_hrv','Phrv','Rm_hrv','Rci_hrv','xPhrv','xRm_hrv','xRci_hrv')
end

P = [p_one_tailed Pval Page Pwght Pbmi Pmage Psex Phrv];
stat = [Z nan Rm_age Rm_wght Rm_bmi Rm_mage Rm_sex Rm_hrv];
cilo = [nan nan Rci_age(1) Rci_wght(1) Rci_bmi(1) Rci_mage(1) Rci_sex(1) Rci_hrv(1)];
cihi = [nan nan Rci_age(2) Rci_wght(2) Rci_bmi(2) Rci_mage(2) Rci_sex(2) Rci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tout = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tout,'SL_fetal_stats_out.csv')

% Control test stats

P = [xp_one_tailed xPval xPage xPwght xPbmi xPmage xPsex xPhrv];
stat = [xZ nan xRm_age xRm_wght xRm_bmi xRm_mage xRm_sex xRm_hrv];
cilo = [nan nan xRci_age(1) xRci_wght(1) xRci_bmi(1) xRci_mage(1) xRci_sex(1) xRci_hrv(1)];
cihi = [nan nan xRci_age(2) xRci_wght(2) xRci_bmi(2) xRci_mage(2) xRci_sex(2) xRci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tctr = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tctr,'SL_fetal_CONTROL_stats_out.csv')

%% save everything
close all
save(sprintf('stat_learning_freqtag_fetal_%s',date))

%% plot data stack

myfigure2

% plot with smoothing

% Wavelet parameters
cfg = [];
cfg.fsample = Fs;
cfg.foi_start = 0.5;
cfg.foi_end = 8;
cfg.window_shift=0.025;
cfg.oct_bw = 0.2;
cfg.oct_delta = cfg.oct_bw/4;
cfg.verbose = false;

[pow_all,foi] = ro_freq_wavelet_TFT2(mean(alldat),cfg);

plot(foi,log10(squeeze(mean(pow_all))),'linewidth',2)
xlabel('Frequency (Hz)')
ylabel('Power (A.U.)')
title('Fetal MEG power spectrum, Experiment 1','fontsize',22)
plot(ones(1,100).*3,linspace(-3.25,-1.5,100),'k--')
makefighandsome
xlim([0.5 8])
ylim([-3.25 -2.4])
print('-dpng','fMEG_power_spectrum_stat_learning_smooth')
print('-dsvg','fMEG_power_spectrum_stat_learning_smooth')


%% helper function


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