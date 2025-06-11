% Performs False Discovery Rate (FDR) correction on frequency tagging 
% statistical results from three datasets (local-global fetal, 
% local-global newborn, and statistical learning fetal), comparing 
% experimental versus control conditions and identifying significant 
% correlations that differ between conditions.

% Joel Frohlich
% Last updated 11.06.2025

clearvars
Nshift = 1000; % so we can repair empirical P-value = 0 

% FDR correct frequency tagging stats

T1 = readtable('LG_fetal_stats_out.csv'); 
T2 = readtable('LG_newborn_stats_out.csv'); 
T3 = readtable('SL_fetal_stats_out.csv');

% Only take the variables common amongst all data, so we can see if they
% replicate across experiments and datasets

vars = intersect(T1.Test,intersect(T2.Test,T3.Test));
T1(~contains(T1.Test,vars),:) = [];
T2(~contains(T2.Test,vars),:) = [];
T3(~contains(T3.Test,vars),:) = [];

% experimental family 
P1e = [T3.Pvalue(1:2); T1.Pvalue(1:2); T2.Pvalue(1:2)];
P2e = [T3.Pvalue(3:end); T1.Pvalue(3:end); T2.Pvalue(3:end)];

Tmain = [T3(1:2,[1 2 5]); T1(1:2,[1 2 5]); T2(1:2,[1 2 5])];
Tsupp = [T3(3:end,:); T1(3:end,:); T2(3:end,:)];

clear T1 T2 T3

T1 = readtable('LG_fetal_CONTROL_stats_out.csv'); 
T2 = readtable('LG_newborn_CONTROL_stats_out.csv'); 
T3 = readtable('SL_fetal_CONTROL_stats_out.csv');

% Only take the variables common amongst all data, so we can see if they
% replicate across experiments and datasets

vars = intersect(T1.Test,intersect(T2.Test,T3.Test));
T1(~contains(T1.Test,vars),:) = [];
T2(~contains(T2.Test,vars),:) = [];
T3(~contains(T3.Test,vars),:) = [];

% control family
P1c = [T3.Pvalue(1:2); T1.Pvalue(1:2); T2.Pvalue(1:2)];
P2c = [T3.Pvalue(3:end); T1.Pvalue(3:end); T2.Pvalue(3:end)];

for ip = 1:length(P1e)
    if P1e(ip) == 0
        P1e(ip) = 1/(Nshift+1);
    end
    if P1c(ip) == 0
        P1c(ip) = 1/(Nshift+1);
    end
end

for ip = 1:length(P2e)
    if P2e(ip) == 0
        P2e(ip) = 1/(Nshift+1);
    end
    if P2c(ip) == 0
        P2c(ip) = 1/(Nshift+1);
    end
end

Q1e = mafdr(P1e,'bhfdr',1);
Q2e = mafdr(P2e,'bhfdr',1);
Q1c = mafdr(P1c,'bhfdr',1);
Q2c = mafdr(P2c,'bhfdr',1);

Tmain = [Tmain; T3(1:2,[1 2 5]); T1(1:2,[1 2 5]); T2(1:2,[1 2 5])];
Tsupp = [Tsupp; T3(3:end,:); T1(3:end,:); T2(3:end,:)];

%%% Finish main table

Tmain.Qvalue = [Q1e; Q1c];

Tmain.data  = repmat([repmat({'Exp. 1 fetal'},2,1); repmat({'Exp. 2 fetal'},2,1); ...
    repmat({'Exp. 2 newborn'},2,1)],2,1);

Tmain.frequency  = [repmat({'Experimental'},6,1); repmat({'Control'},6,1)];
 
Tmain.stat = round(Tmain.stat,3,'significant');
Tmain.Pvalue = round(Tmain.Pvalue,3,'significant');
Tmain.Qvalue = round(Tmain.Qvalue,3,'significant');

Tmain.Significance = cell(size(Tmain,1),1);

for irow = 1:size(Tmain,1)
    if Tmain.Qvalue(irow) < 0.0005
        Tmain.Significance{irow} = '***';
    elseif Tmain.Qvalue(irow) < 0.005
        Tmain.Significance{irow} = '**';
    elseif Tmain.Qvalue(irow) < 0.05
        Tmain.Significance{irow} = '*';
    else
        main.Significance{irow} = ' ';
    end
end

% Round results and export latex table
Tmain{:,2:4} = round(Tmain{:,2:4},3,'significant');
table2latex(Tmain, 'Main_table.tex')

%%% Finish supplemental table

Tsupp.Qvalue = [Q2e; Q2c];

Tsupp.data  = repmat([repmat({'Exp. 1 fetal'},6,1); repmat({'Exp. 2 fetal'},6,1); ...
    repmat({'Exp. 2 newborn'},6,1)],2,1);

Tsupp.frequency  = [repmat({'Experimental'},18,1); repmat({'Control'},18,1)];
 
Tsupp.stat = round(Tsupp.stat,3,'significant');
Tsupp.Pvalue = round(Tsupp.Pvalue,3,'significant');
Tsupp.Qvalue = round(Tsupp.Qvalue,3,'significant');

Tsupp.Significance = cell(size(Tsupp,1),1);

for irow = 1:size(Tsupp,1)
    if Tsupp.Qvalue(irow) < 0.0005
        Tsupp.Significance{irow} = '***';
    elseif Tsupp.Qvalue(irow) < 0.005
        Tsupp.Significance{irow} = '**';
    elseif Tsupp.Qvalue(irow) < 0.05
        Tsupp.Significance{irow} = '*';
    else
        main.Significance{irow} = ' ';
    end
end

% Round results and export latex table
Tsupp{:,2:6} = round(Tsupp{:,2:6},3,'significant');
table2latex(Tsupp, 'Supplemental_table.tex')
