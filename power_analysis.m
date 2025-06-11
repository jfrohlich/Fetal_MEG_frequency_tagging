% Performs a Monte Carlo subsampling-based statistical power analysis for
% frequency tagging data, testing how effect sizes (Z-scores) change with 
% different sample sizes by repeatedly drawing random subsets of subjects 
% and calculating frequency tagging statistics.
%
% Joel Frohlich
% Last updated 11.06.2025

clear all
try
    parpool
catch
    fprintf('Parpool is probably already open\n')
end

%%% Load the most recent output %%%
% Uncomment below with file name, include date of most recent output
% load('stat_learning_freqtag_fetal_[DATE].mat','testfreq','Fs','control','alldat')

% stat test of amplitude at f Hz
N = round(Fs)*20;
foi = linspace(0,round(Fs)/2,N/2);
nit = 1000;
Z = nan(size(alldat,1),nit);


for i = 1:size(alldat,1)
    parfor j = 1:nit
        rng(j,'twister')
        drawme = randperm(size(alldat,1),i);
        
        tmp = abs(fft(mean(alldat(drawme,:)),N));
        amp_spectrum = tmp(1:N/2);
        ratio = round(Fs)/Fs;
        
        [~,where] = min(abs(foi - testfreq));
        left = amp_spectrum(where-11:where-2);
        right = amp_spectrum(where+2:where+11);
        sigma = std([left right]);
        mu = mean([left right]);
        Z(i,j) = (amp_spectrum(where) - mu)/sigma;
    end
end

%% Plotting

Zmed = median(Z,2);
Pmed =  1 - normcdf(Zmed);
Z10 = prctile(Z',10);
Z20 = prctile(Z',20);
Z50 = prctile(Z',50);

close all

myfigure

plot(1:size(alldat,1),Z50,'linewidth',2)
plot(1:size(alldat,1),Z20,'linewidth',2)
plot(1:size(alldat,1),Z10,'linewidth',2)
xlabel('Sample size (N)')
ylabel('Z(SNR)')
plot(0:size(alldat,1),ones(1,size(alldat,1)+1).*norminv(1 - 0.05),'k')
plot(0:size(alldat,1),ones(1,size(alldat,1)+1).*norminv(1 - 0.005),'k--')
plot(0:size(alldat,1),ones(1,size(alldat,1)+1).*norminv(1 - 0.0005),'k:')
legend({'50% power','80% power','90% power','P = 0.05','P = 0.005','P = 0.0005'},'fontsize',16,'box','off','location','southeast')
title('Subsampling-based power analysis (Experiment 1)','fontsize',18)
makefighandsome
print('-dsvg','PowerAnalysis.svg')
print('-dpng','PowerAnalysis.png')



