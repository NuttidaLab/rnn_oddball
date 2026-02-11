% Name: Robert Kim, Nuttida Rungratsameetaweemana
% Date: 07-21-2025
% Email: robert.f.kim@gmail.com, nr2869@columbia.edu
% For plotting the lesion results

clear; clc;

fig_dir = fullfile(pwd, '..', '..', 'figures');
result_dir = fullfile(pwd, '..', '..', 'results');

intact_results = load(fullfile(result_dir, 'svm_results.mat'));

% IE
ie_results = load(fullfile(result_dir, 'ie_svm_results_v2.mat'));

% EE
ee_results = load(fullfile(result_dir, 'ee_svm_results_v2.mat'));

% EI
ei_results = load(fullfile(result_dir, 'ei_svm_results_v2.mat'));

% II
ii_results = load(fullfile(result_dir, 'ii_svm_results_v2.mat'));


tone_accs = [intact_results.all_tone_acc, ie_results.all_tone_acc, ee_results.all_tone_acc, ei_results.all_tone_acc, ii_results.all_tone_acc];
ctx_accs  = [intact_results.all_ctx_acc, ie_results.all_ctx_acc, ee_results.all_ctx_acc, ei_results.all_ctx_acc, ii_results.all_ctx_acc];

labs = [ones(1, 10)*0, ones(1, 10)*1, ones(1, 10)*2, ones(1, 10)*3, ones(1, 10)*4];

figure;
boxplot(tone_accs, labs, 'symbol', '')
xticks(1:5)
xticklabels({'Intact', 'E->I', 'E->E', 'I->E', 'I->I'})

% Kruskal-Wallis test
[p_kw, tbl, stats] = kruskalwallis(tone_accs, labs, 'off');

% Post-hoc Dunn-Sidak correction
c = multcompare(stats, 'CType', 'dunn-sidak', 'Display', 'off');

%print(fullfile(fig_dir, 'tone_svm_lesion.eps'), '-painters', '-depsc');

