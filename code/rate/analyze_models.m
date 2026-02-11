% Name: Robert Kim, Nuttida Rungratsameetaweemana
% Date: 3/5/2025
% Email: robert.f.kim@gmail.com, nr2869@columbia.edu
% Script to analyze the models to plot the average performance and to get some numbers
% (# of neurons that switch coding)

clear; clc;
rng(820)

fig_dir = fullfile(pwd, '..', '..', 'figures', 'manuscript');
result_dir = fullfile(pwd, '..', '..', 'results');
model_path = fullfile(pwd, '..', '..', 'models', 'go-nogo2', 'P_rec_0.2_Taus_4.0_25.0');

% N = 200 models
model_ids = {'194212', '153216', '201555', '183934', '234618', '193127', '195849', '153730', '195329', '184637'}; % trained model IDs. Replace with the IDs that you trained
model_ind = [1:10]; % just include all 10 trained models

svm_results = load(fullfile(result_dir, 'svm_results.mat'));
shuffled_svm_results = load(fullfile(result_dir, 'shuffled_svm_results_v2.mat'));
tone_dec_ctx_inc_N = svm_results.tone_dec_ctx_inc_N(model_ind)/200*100;
tone_dec_ctx_dec_N = svm_results.tone_dec_ctx_dec_N(model_ind)/200*100;
tone_inc_ctx_dec_N = svm_results.tone_inc_ctx_dec_N(model_ind)/200*100;
tone_inc_ctx_inc_N = svm_results.tone_inc_ctx_inc_N(model_ind)/200*100;

tone_inc_N = svm_results.tone_inc_N(model_ind)/200*100;
tone_dec_N = svm_results.tone_dec_N(model_ind)/200*100;
ctx_inc_N = svm_results.ctx_inc_N(model_ind)/200*100;
ctx_dec_N = svm_results.ctx_dec_N(model_ind)/200*100;

unchanged_N = svm_results.unchanged_N(model_ind)/200*100;

disp(['% of Units whose decoding did not change: ' num2str(mean(unchanged_N)) '+/-' num2str(std(unchanged_N'))]);

disp(['% of Units whose tone decoding decreased: ' num2str(mean(tone_dec_N)) '+/-' num2str(std(tone_dec_N'))]);
disp(['% of Units whose tone decoding increased: ' num2str(mean(tone_inc_N)) '+/-' num2str(std(tone_inc_N'))]);
disp(['% of Units whose context decoding decreased: ' num2str(mean(ctx_dec_N)) '+/-' num2str(std(ctx_dec_N'))]);
disp(['% of Units whose context decoding increased: ' num2str(mean(ctx_inc_N)) '+/-' num2str(std(ctx_inc_N'))]);

disp(['% of Units whose tone decoding decreased AND ctx decoding decreased: ' num2str(mean(tone_dec_ctx_dec_N)) '+/-' num2str(std(tone_dec_ctx_dec_N'))]);
disp(['% of Units whose tone decoding decreased AND ctx decoding increased: ' num2str(mean(tone_dec_ctx_inc_N)) '+/-' num2str(std(tone_dec_ctx_inc_N'))]);
disp(['% of Units whose tone decoding increased AND ctx decoding decreased: ' num2str(mean(tone_inc_ctx_dec_N)) '+/-' num2str(std(tone_inc_ctx_dec_N'))]);
disp(['% of Units whose tone decoding increased AND ctx decoding increased: ' num2str(mean(tone_inc_ctx_inc_N)) '+/-' num2str(std(tone_inc_ctx_inc_N'))]);
disp('paused!')

[p, h] = signrank(tone_dec_N, tone_inc_N, 'tail', 'right')
[p, h] = signrank(ctx_dec_N, ctx_inc_N, 'tail', 'left')
pause;

A = round(sum(tone_dec_ctx_inc_N))
B = round(sum(tone_dec_ctx_dec_N))
C = round(sum(tone_inc_ctx_inc_N))
D = round(sum(tone_inc_ctx_dec_N))
contingencyTable = [A, B; C, D];
[pValue,~,stats] = fishertest(contingencyTable);

%% Boxplot of overall accuracy (tone and context)
figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
hold on;
boxplot([svm_results.all_tone_acc(model_ind); svm_results.all_ctx_acc(model_ind)]', 'symbol', '');
ylim([0.3 0.8]);
xlim([0.5, 2.5]);
plot(xlim, [0.5, 0.5], 'k--');
xticks([1 2]);
xticklabels({'Tone', 'Context'});
ylabel('Decoding accuracy');
%print(fullfile(fig_dir, 'overall_decoding_boxplot.eps'), '-painters', '-depsc');
[p, h] = ranksum(svm_results.all_tone_acc(model_ind), svm_results.all_ctx_acc(model_ind))

% Comparing against shuffled condition
[p, h] = signrank(svm_results.all_tone_acc(model_ind), shuffled_svm_results.all_tone_acc)
[p, h] = signrank(svm_results.all_ctx_acc(model_ind), shuffled_svm_results.all_ctx_acc)


%%
figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
hold on;
boxplot(svm_results.all_ntpc_ctx_acc', 'symbol', '');
ylim([0.4 0.70])
xlim([0.5, 2.5])
plot(xlim, [0.5, 0.5], 'k--');
xticks([1 2]);
xticklabels({'Early Trials', 'Late Trials'});
ylabel('Decoding accuracy');
title('Decoding context in ')
[p, h] = signrank(svm_results.all_ntpc_ctx_acc(1, :), svm_results.all_ntpc_ctx_acc(2, :))
%print(fullfile(fig_dir, 'ctx_decoding_boxplot.eps'), '-painters', '-depsc');

figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
hold on;
boxplot(svm_results.all_ntpc_tone_acc', 'symbol', '');
ylim([0.40 0.90])
xlim([0.5, 2.5])
plot(xlim, [0.5, 0.5], 'k--');
xticks([1 2]);
xticklabels({'Early Trials', 'Late Trials'});
ylabel('Decoding accuracy');
[p, h] = signrank(svm_results.all_ntpc_tone_acc(1, :), svm_results.all_ntpc_tone_acc(2, :))
%print(fullfile(fig_dir, 'tone_decoding_boxplot.eps'), '-painters', '-depsc');


%%
all_data = svm_results.all_data;
all_tone_inc_decoding = [];
all_tone_dec_decoding = [];
all_ctx_inc_decoding = [];
all_ctx_dec_decoding = [];
all_rest_tone_decoding = [];
all_rest_ctx_decoding = [];

for i = 1:length(all_data)
    fdr_right_tone_ps = mafdr(all_data(i).right_tone_ps, 'BHFDR', true);
    fdr_left_tone_ps = mafdr(all_data(i).left_tone_ps, 'BHFDR', true);
    tone_dec_units = find(fdr_right_tone_ps < 0.05);
    tone_inc_units = find(fdr_left_tone_ps < 0.05);

    tone_dec_decoding = mean(all_data(i).all_tone_decoding(tone_dec_units, :));
    tone_inc_decoding = mean(all_data(i).all_tone_decoding(tone_inc_units, :));
    all_tone_dec_decoding = [all_tone_dec_decoding; tone_dec_decoding];
    all_tone_inc_decoding = [all_tone_inc_decoding; tone_inc_decoding];

    fdr_left_ctx_ps = mafdr(all_data(i).left_ctx_ps, 'BHFDR', true);
    fdr_right_ctx_ps = mafdr(all_data(i).right_ctx_ps, 'BHFDR', true);
    ctx_inc_units = find(fdr_left_ctx_ps < 0.05);
    ctx_dec_units = find(fdr_right_ctx_ps < 0.05);

    if length(ctx_dec_units) > 1
        ctx_dec_decoding = mean(all_data(i).all_ctx_decoding(ctx_dec_units, :));
    else
        ctx_dec_decoding = all_data(i).all_ctx_decoding(ctx_dec_units, :);
    end
    if length(ctx_inc_units) > 1
        ctx_inc_decoding = mean(all_data(i).all_ctx_decoding(ctx_inc_units, :));
    else
        ctx_inc_decoding = all_data(i).all_ctx_decoding(ctx_inc_units, :);
    end
    all_ctx_dec_decoding = [all_ctx_dec_decoding; ctx_dec_decoding];
    all_ctx_inc_decoding = [all_ctx_inc_decoding; ctx_inc_decoding];

    temp_uni = union(tone_dec_units, tone_inc_units);
    temp_uni = union(temp_uni, ctx_inc_units);
    temp_uni = union(temp_uni, ctx_dec_units);
    unchanged_units = setdiff([1:200], temp_uni);

    if length(unchanged_units) > 1
        rest_tone_decoding = mean(all_data(i).all_tone_decoding(unchanged_units, :));
        rest_ctx_decoding = mean(all_data(i).all_ctx_decoding(unchanged_units, :));
    else
        rest_tone_decoding = all_data(i).all_tone_decoding(unchanged_units, :);
        rest_ctx_decoding = all_data(i).all_ctx_decoding(unchanged_units, :);
    end

    all_rest_tone_decoding = [all_rest_tone_decoding; rest_tone_decoding];
    all_rest_ctx_decoding = [all_rest_ctx_decoding; rest_ctx_decoding];
end

figure('Units', 'Inch', 'Outerposition', [0 0 3 4]);
hold on;
purple = [0.42, 0.39, 0.72];
green  = [0.25, 0.58, 0.29];
xVals = [1 2];
errorbar(xVals, nanmean(all_ctx_inc_decoding), sem(all_ctx_inc_decoding), '-o', ...
         'Color', green, 'LineWidth',2, 'CapSize',10, 'DisplayName', '', 'HandleVisibility', 'off');
plot(nan, nan, '-', 'Color', green, 'linewidth', 2, 'DisplayName', 'Oddball');
errorbar(xVals, nanmean(all_tone_dec_decoding), sem(all_tone_dec_decoding), '-o', ...
         'Color', purple, 'LineWidth',2, 'CapSize',10, 'DisplayName', '', 'HandleVisibility', 'off');
plot(nan, nan, '-', 'Color', purple, 'linewidth', 2, 'DisplayName', 'Frequency');

legend('Location','northeast');
xlim([0.5, 2.5])
ylim([0.6, 0.80])
xticks([1 2]);
xticklabels({'Early Trials', 'Late Trials'});
%print(fullfile(fig_dir, 'ctx_inc_tone_dec_lineplot.eps'), '-painters', '-depsc');

figure('Units', 'Inch', 'Outerposition', [0 0 3 4]);
hold on;
purple = [0.42, 0.39, 0.72];
green  = [0.25, 0.58, 0.29];
light_purple = 0.5*purple + 0.5*[1, 1, 1];
light_green = 0.5*green  + 0.5*[1, 1, 1];
xVals = [1 2];
errorbar(xVals, nanmean(all_ctx_dec_decoding), sem(all_ctx_dec_decoding), '-o', ...
         'Color', green, 'LineWidth',2, 'CapSize',10, 'DisplayName', '', 'HandleVisibility', 'off');
errorbar(xVals, nanmean(all_rest_ctx_decoding), sem(all_rest_ctx_decoding), '--o', ...
         'Color', light_green, 'LineWidth',2, 'CapSize',10, 'DisplayName', '', 'HandleVisibility', 'off');
plot(nan, nan, '-', 'Color', green, 'linewidth', 2, 'DisplayName', 'Oddball');
plot(nan, nan, '--', 'Color', green, 'linewidth', 2, 'DisplayName', 'Oddball, unchanged');

errorbar(xVals, nanmean(all_tone_inc_decoding), sem(all_tone_inc_decoding), '-o', ...
         'Color', purple, 'LineWidth',2, 'CapSize',10, 'DisplayName', '', 'HandleVisibility', 'off');
errorbar(xVals, nanmean(all_rest_tone_decoding), sem(all_rest_tone_decoding), '--o', ...
         'Color', light_purple, 'LineWidth',2, 'CapSize',10, 'DisplayName', '', 'HandleVisibility', 'off');
plot(nan, nan, '-', 'Color', purple, 'linewidth', 2, 'DisplayName', 'Frequency');
plot(nan, nan, '--', 'Color', purple, 'linewidth', 2, 'DisplayName', 'Frequency, unchanged');

legend('Location','northeast');
xlim([0.5, 2.5])
ylim([0.45, 1.10])
xticks([1 2]);
xticklabels({'Early Trials', 'Late Trials'});
%print(fullfile(fig_dir, 'Supplementary_ctx_dec_tone_inc_lineplot.eps'), '-painters', '-depsc');

