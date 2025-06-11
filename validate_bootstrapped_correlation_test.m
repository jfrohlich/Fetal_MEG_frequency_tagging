%% Compute Type I and Type II error rate for boot-strapped correlation test
clearvars
close all

rng(111)
tic
siglen = 6000;
Fs= 600; % pretend sampling rate
t = linspace(0,siglen/Fs,siglen);
testfreq = 3;
x = sin(t.*2*pi*testfreq);
nsub = 50;
alldat = nan(nsub,siglen);
coefs = rand(1,size(alldat,1))./50; % create random weights and scale them
N = round(Fs)*20;
foi = linspace(0,round(Fs)/2,N/2);


for irow = 1:nsub
    alldat(irow,:) = zscore(powernoise(2,siglen))' +  x.*coefs(irow);
end


ntrials = 2000;
rsvec = 10:10:500;
nrs = max(rsvec);
p1 = nan(ntrials,length(rsvec));
p2 = nan(ntrials,length(rsvec));

for itrial = 1:ntrials
    cx = nan(1,nrs);
    snrx = nan(1,nrs);
    zzx = nan(1,nrs);
    for irs = 1:nrs
        NDX = randi(nsub,1,nsub);
        tmp = abs(fft(mean(alldat(NDX,:)),N));
        cx(irs) = mean(coefs(NDX));
        amp_spectrum = tmp(1:N/2);
        [~,where] = min(abs(foi - testfreq));
        left = amp_spectrum(where-11:where-2);
        right = amp_spectrum(where+2:where+11);
        snrx(irs) = amp_spectrum(where)/mean([left right]);
        sigma = std([left right]);
        mu = mean([left right]);
        zzx(irs) = (amp_spectrum(where) - mu)/sigma;
    end
    for icnt = 1:length(rsvec)
        irs = rsvec(icnt);
        [~,p1(itrial,icnt)] = corr(cx(randperm(irs))',zzx(1:irs)');
        [~,p2(itrial,icnt)] = corr(cx(1:irs)',zzx(1:irs)');
    end
    fprintf('%1.1f%% complete\n',itrial/ntrials*100)
end

%% plotting

typeIerror = [];
typeIIerror = [];


for irs = 1:size(p1,2)
    typeIerror(irs) = sum(p1(:,irs) < 0.05)/ntrials;
    typeIIerror(irs) = sum(p2(:,irs) >= 0.05)/ntrials;
end

%% plot error rates
myfigure2
plot(rsvec,typeIerror,'linewidth',2)
title('parametric Type I error, alpha = 0.05','FontSize',20)
xlabel('Number of resamples per correlation')
ylabel('False positive rate')
ylim([0 0.15])
makefighandsome
print('-dpng','Parametric_TypeIerror_resample_test.png')
print('-dsvg','Parametric_TypeIerror_resample_test.svg')

myfigure2
plot(rsvec,typeIIerror,'linewidth',2)
title('parametric Type II error, alpha = 0.05','FontSize',20)
xlabel('Number of resamples per correlation')
ylabel('False negative rate')
ylim([-0.01 1])
makefighandsome
print('-dpng','Parametric_TypeIIerror_resample_test.png')
print('-dsvg','Parametric_TypeIIerror_resample_test.svg')


%% Type I error report

a=round(min(typeIerror),2,'significant');
b=round(max(typeIerror),2,'significant');
c=round(median(typeIerror),2,'significant');

fprintf('Type I error, min = %1.3f, max = %1.3f, median = %1.3f\n',a,b,c)

%% Type II error report

a=round(min(typeIIerror),2,'significant');
b=round(max(typeIIerror),2,'significant');

fprintf('Type II error, min = %1.3g, max = %1.3f\n',a,b)
save boostrap_validation_improved
