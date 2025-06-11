% Frequency tagging analysis on neonatal MEG data from local-global 
% experiments, testing for neural responses at the stimulation frequency 
% (5/3 Hz) versus a control frequency (3 Hz), while examining correlations
% with demographic and physiological variables through bootstrapped 
% resampling.
% 
% Joel Frohlich
% 10.06.2025

clearvars  
load time_traces_data_manuscript.mat       
load Alltable.mat

T = Alltable(logical(Alltable.BornYet),:);
ntotal = 20; % total number of neonatal recordings


subs = 1:max(T.ID);
cnt = 0;
cnta = 0;
cntb = 0;
NDX = nan(4,1);
NDXA = nan(2,1);
NDXB = nan(2,1);
subjectIDs = [];
fsample = 610.3516;
m800 = round(0.8*fsample);

testfreq = 5/3*(round(fsample)/fsample); % do it this way because fft only takes integer sampling rates
control = 3*(round(fsample)/fsample); % do it this way because fft only takes integer sampling rates;  % control frequency in Hz (from the other experiment)
Nshift = 1000;


for isub = subs
    % filter on subjects with good data in all four conditions
    nx = sum(T.ID == isub);
    
    if nx == 4
        cnt = cnt + 1;
        NDX(:,cnt) = find(contains(data_tab.ID, sprintf('0%i-',subs(isub))) | ...
            contains(data_tab.ID, sprintf('o%i-',subs(isub))));
        subjectIDs = [subjectIDs subs(isub)];
    elseif nx == 2
        idx = find(contains(data_tab.ID, sprintf('0%i-',subs(isub))) | ...
            contains(data_tab.ID, sprintf('o%i-',subs(isub))));
        tmp = data_tab.condition(idx);
        
        if strcmp(tmp{1},'local_dev') % block B (trained on sssd)
            cntb = cntb + 1;
            subjectIDs = [subjectIDs subs(isub)];
            NDXB(:,cntb) = idx;
        elseif strcmp(tmp{1},'sandard') % [sic], block A (trained on ssss)
            cnta = cnta + 1;
            subjectIDs = [subjectIDs subs(isub)]; 
            NDXA(:,cnta) = idx;
        else
            error('Unexpected condition')
        end
    elseif nx ~= 0
        error('Unexpected data')
    end
end

%%
ssss = [];
sssd = [];
sssS = [];
sssD = [];

n4 = size(NDX,2); % who has data in all four conditions?

for isub = 1:size(NDX,2)
    ssss = [ssss; data_tab.data(NDX(1,isub),:)];
    sssD = [sssD; data_tab.data(NDX(2,isub),:)];
    sssd = [sssd; data_tab.data(NDX(3,isub),:)];
    sssS = [sssS; data_tab.data(NDX(4,isub),:)];
end

for isub = 1:size(NDXA,2)
    ssss = [ssss; data_tab.data(NDXA(1,isub),:)];
    sssD = [sssD; data_tab.data(NDXA(2,isub),:)];
end

for isub = 1:size(NDXB,2)
    sssd = [sssd; data_tab.data(NDXB(1,isub),:)];
    sssS = [sssS; data_tab.data(NDXB(2,isub),:)];
end

assert(size(sssd,1) == size(sssS,1),'unbalanced data')
assert(size(sssD,1) == size(ssss,1),'unbalanced data')

whichones =[];
for isub = 1:length(subjectIDs)
    whichones(isub) = find(subjectIDs(isub) == T.ID,1,'first');
end

AGE = T.GA(whichones);
HRV = 10.^T.SDNN(whichones);

alldat = [ssss(1:n4,1:end-m800) sssD(1:n4,1:end-m800) sssd(1:n4,1:end-m800) sssS(1:n4,1:end-m800); ...
    ssss(n4+1:end,1:end-m800) sssD(n4+1:end,1:end-m800) ssss(n4+1:end,1:end-m800) sssD(n4+1:end,1:end-m800); ...
    sssd(n4+1:end,1:end-m800) sssS(n4+1:end,1:end-m800) sssd(n4+1:end,1:end-m800) sssS(n4+1:end,1:end-m800)];

alldat = zscore(alldat')';
assert(size(alldat,1) == ntotal,'Wrong number of recordings')
avgdat = mean(alldat);

figure
plot(avgdat)
print('-dpng','neo_amplitude_spectrum_localglobal_fft.png')
print('-dsvg','neo_amplitude_spectrum_localglobal_fft.svg')


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
ylabel('Amplitude (Z-score)')
title('Neonatal MEG FFT amplitude spectrum, Experiment 2','fontsize',22)
plot(ones(1,100).*(testfreq),linspace(0,1e5,100),'k--')
legend({'Amplitude spectrum','Stimulation frequency'},'fontsize',16,...
    'location','northeast','box','off')
makefighandsome
xlim([0.5 8])
ylim([0 1000])
print('-dpng','neo_fMEG_amplitude_spectrum_localglobal_fft.png')
print('-dsvg','neo_fMEG_amplitude_spectrum_localglobal_fft.svg')

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

% Control test

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


%% Visualize time-averged signal, this will help us confirm that it's not an instantaneous artifact

myfigure2
time = linspace(0,9.6,length(avgdat));
plot(time,avgdat,'linewidth',2)
beep = 0.2:0.6:9.2;
for ibeep = 1:length(beep)
    plot(repmat(beep(ibeep),1,100),linspace(min(avgdat)*1.5,max(avgdat)*1.5,100),'k:')
    if ibeep == 1 
        legend({'MEG signal','Stimulus onset'},'autoupdate','off',...
            'fontsize',18,'location','northeast')
        legend box off
    end
end

xlabel('time (sec)')
ylabel('Field strength (AU)')
title('Grand-averaged neonatal MEG signal, Experiment 2','fontsize',22)
makefighandsome

print('-dsvg','Neural_response_tones_LG_newborn.svg')


%% correlate SNR with age and HRV

assert(mean(round(HRV)) > 1, 'HRV is still log-transformed!')

% extract fetal sex data
where1 =[];
for isub = 1:length(subjectIDs)
    where1(isub) = find(subjectIDs(isub) == T.ID,1,'first');
end
where2 =[];
for isub = 1:length(subjectIDs)
    where2(isub) = find(subjectIDs(isub) == Alltable.ID,1,'first');
end


%%% Update the table with BMIs that were missing
IDn = unique(Alltable.ID(logical(Alltable.BornYet)));
IDf = unique(Alltable.ID(~logical(Alltable.BornYet)));
IDx = setdiff(IDn,IDf);
% missing BMI data -- the first column gives the subject ID, the second
% column gives the maternal BMI (prepregnancy) 
bmidat = [34 65/1.68^2; 26 56/1.66^2; 27 61/1.67^2; 6 72.5/1.66^2];

for i = 1:length(IDx)
    idx = find(bmidat(:,1) == IDx(i));
    for j = 1:size(Alltable,1)
        if Alltable.ID(j) == IDx(i)
            Alltable.BMI(j) = bmidat(bmidat(:,1)==IDx(i),2)
        end
    end
end

SEX = T.XChrom(where1);
MAGE = T.MomsAge(where1);
WGHT = T.BirthWeight(where1);
BMI = Alltable.BMI(where2);


assert(~any(BMI==0),'missing bmi')

nsub = size(alldat,1);
nboot = 1000;
nrs = 500;

ages = nan(nboot, nrs);
sexS = nan(nboot, nrs);
bageS = nan(nboot, nrs);
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

%% calculate each correlation

Rhrv = diag(corr(hrvS',zz','rows','complete','type',ctype));
Rage = diag(corr(ages',zz','rows','complete','type',ctype));
Rsex = diag(corr(sexS',zz','rows','complete','type',ctype));
Rmage = diag(corr(mageS',zz','rows','complete','type',ctype));
Rwght = diag(corr(wghtS',zz','rows','complete','type',ctype));
Rbmi = diag(corr(bmiS',zz','rows','complete','type',ctype));

% Compute non-parametric P-values

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


%% Finally, repeat the correlational analysis using a control frequency (3 Hz)

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

% Compute non-parametric P-values

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


%% DEBUG %% remove this later
% Output table
Pval = nan;
P = [p_one_tailed Pval Page Pwght Pbmi Pmage Psex Phrv];
stat = [Z nan Rm_age Rm_wght Rm_bmi Rm_mage Rm_sex Rm_hrv];
cilo = [nan nan Rci_age(1)  Rci_wght(1) Rci_bmi(1) Rci_mage(1) Rci_sex(1) Rci_hrv(1)];
cihi = [nan nan Rci_age(2)  Rci_wght(2) Rci_bmi(2) Rci_mage(2) Rci_sex(2) Rci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tout = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tout,'LG_newborn_stats_out_DEBUG.csv')

% Control test stats
xPval= nan;
P = [xp_one_tailed xPval xPage xPwght xPbmi xPmage xPsex xPhrv];
stat = [xZ nan xRm_age xRm_wght xRm_bmi xRm_mage xRm_sex xRm_hrv];
cilo = [nan nan xRci_age(1) xRci_wght(1) xRci_bmi(1) xRci_mage(1) xRci_sex(1) xRci_hrv(1)];
cihi = [nan nan xRci_age(2) xRci_wght(2) xRci_bmi(2) xRci_mage(2) xRci_sex(2) xRci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tctr = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tctr,'LG_newborn_CONTROL_stats_out_DEBUG.csv')
return
%%%
%% frequency transform of surrogates

shiftedData = nan(size(alldat,1),size(alldat,2),Nshift);
parfor isub = 1:nsub
    rng(isub, 'Twister'); % Set the random seed for reproducibility
    shiftedData(isub,:,:) = surrogate(alldat(isub,:), Nshift, 'IAAFT1', 0, fsample)';
    fprintf('%1.2f%% complete\n',isub/nsub*100)
end


[~,where] = min(abs(foi - testfreq));
SNRshift = nan(1,Nshift);
parfor ishift = 1:Nshift
    x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
    left = x_spectrum(where-flank:where-2);
    right = x_spectrum(where+2:where+flank);
    sigma = std([left right]);
    mu = mean([left right]);
    SNRshift(ishift) = (x_spectrum(where) - mu)/sigma;
    fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
end


% non-parametric
Pval = sum(SNRshift >= Z)/Nshift;
fprintf('Frequency tagging effect, non-parametric test, P = %1.2f\n',Pval)

% control test
[~,where] = min(abs(foi - control));
xSNRshift = nan(1,Nshift);
% now compare SNR with the surrogate spectra
parfor ishift = 1:Nshift
    x_spectrum = abs(fft(mean(shiftedData(:,:,ishift)),N));
    left = x_spectrum(where-flank:where-2);
    right = x_spectrum(where+2:where+flank);
    sigma = std([left right]);
    mu = mean([left right]);
    xSNRshift(ishift) = (x_spectrum(where) - mu)/sigma;
    fprintf('FFT on surrogate data %1.2f%% complete\n',ishift/Nshift*100)
end

% non-parametric
xPval = sum(xSNRshift >= xZ)/Nshift;
fprintf('CONTROL: Frequency tagging effect, non-parametric test, P = %1.2f\n',xPval)


% save everything
save(sprintf('local_global_freqtag_newborn_MAXDATA_%s',date))

%% Output table

P = [p_one_tailed Pval Page Pwght Pbmi Pmage Psex Phrv];
stat = [Z nan Rm_age Rm_wght Rm_bmi Rm_mage Rm_sex Rm_hrv];
cilo = [nan nan Rci_age(1)  Rci_wght(1) Rci_bmi(1) Rci_mage(1) Rci_sex(1) Rci_hrv(1)];
cihi = [nan nan Rci_age(2)  Rci_wght(2) Rci_bmi(2) Rci_mage(2) Rci_sex(2) Rci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tout = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tout,'LG_newborn_stats_out.csv')

% Control test stats

P = [xp_one_tailed xPval xPage xPwght xPbmi xPmage xPsex xPhrv];
stat = [xZ nan xRm_age xRm_wght xRm_bmi xRm_mage xRm_sex xRm_hrv];
cilo = [nan nan xRci_age(1) xRci_wght(1) xRci_bmi(1) xRci_mage(1) xRci_sex(1) xRci_hrv(1)];
cihi = [nan nan xRci_age(2) xRci_wght(2) xRci_bmi(2) xRci_mage(2) xRci_sex(2) xRci_hrv(2)];
vars = {'parametric','non-parametric','age','birth_weight','BMI','maternal_age','Sex','HRV'};

Tctr = table(vars',stat',cilo',cihi',P','VariableNames',{'Test','stat','CI_lower_bound','CI_upper_bound','Pvalue'})

writetable(Tctr,'LG_newborn_CONTROL_stats_out.csv')


%% plot with smoothing 

% Wavelet parameters
cfg = [];
cfg.fsample = Fs;
cfg.foi_start = 0.5;
cfg.foi_end = 8;
cfg.window_shift=0.025;
cfg.oct_bw = 0.2;
cfg.oct_delta = cfg.oct_bw/4;
cfg.verbose = false;

[pow_all,foi] = ro_freq_wavelet_TFT2(avgdat,cfg);

myfigure2
plot(foi,log10(squeeze(mean(pow_all))),'linewidth',2)
xlabel('Frequency (Hz)')
ylabel('Power (Z2/Hz)')
title('Neonatal MEG power spectrum, Experiment 2','fontsize',18)
plot(ones(1,100).*testfreq,linspace(-3,-1,100),'k--')
legend({'MEG power spectrum','Stimulation frequency'},'box','off','fontsize',16,'autoupdate','off')
makefighandsome
xlim([0.5 8])
ylim([-3 -1])
print('-dpng','neo_power_spectrum_localglobal_smooth')
print('-dsvg','neo_power_spectrum_localglobal_smooth')

%% Figures for NYU talk

myfigure2
for irow = 1:size(alldat,1)
    plot(alldat(irow,:)+2.5*irow,'k')
end

plot((avgdat-1.25).*5,'linewidth',3,'color','k')
axis off
print('-dpng','stack_and_avg.png')
print('-dsvg','stack_and_avg.svg')

%
myfigure
plot(shiftedData(1,:,1),'k')
plot(shiftedData(1,:,2)+5,'k')
axis off
print('-dpng','surrogates.png')
print('-dsvg','surrogates.svg')

myfigure
[pow2,foi] = ro_freq_wavelet_TFT2(shiftedData(1,:,1),cfg);
fooof_r2 = fooof(foi', squeeze(mean(pow2)), f_range, settings,true);
amp_spectrum2 = sqrt(10.^(fooof_r2.power_spectrum - fooof_r2.ap_fit));
plot(foi,amp_spectrum2,'linewidth',2)
xlabel('Frequency (Hz)')
ylabel('Amplitude (Z-score)')
title('Surrogate spectrum','fontsize',18)
%plot(ones(1,100).*testfreq,linspace(0,3,100),'k--')
makefighandsome
xlim([0.5 8])
ylim([0.75 1.75])
print('-dpng','surrogate_spectrum.png')
print('-dsvg','surrogate_spectrum.svg')



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
