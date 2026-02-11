% Name: Robert Kim, Nuttida Rungratsameetaweemana
% Date: 05-28-2023
% Email: robert.f.kim@gmail.com, nr2869@columbia.edu
% Description: GLM script

clear; clc;
rng(820);

fig_dir = <DIRECTORY WHERE FIGURES WILL BE SAVED>
model_path = <WHERE TRAINED RNNS ARE SAVED>
model_id = '153216';

all_fnames = dir(fullfile(model_path, ['Stage_*' model_id '*']));
all_fnames = {all_fnames.name};
trs = [101:100:length(all_fnames)*100];

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 200; % trial duration (timesteps)
task_info.stim_on = 50;
task_info.stim_dur = 25;
task_info.prob = 0.50;
stim_off = task_info.stim_on + task_info.stim_dur + 1;

all_coeffs = zeros(length(trs), 200, 3);
all_ps = zeros(length(trs), 200, 3);

curr_neus = 1:200;
for i = 1:length(trs)
  %curr_fi = ['Stage_0_Tr_' num2str(trs(i)) '_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_153216.mat'];
  curr_fi = find(contains(all_fnames, ['_' num2str(trs(i)) '_']));

  full_path = fullfile(model_path, all_fnames{curr_fi});
  temp = load(full_path);

  % Tone B is oddball
  task_info.prob = 0.80;
  num_trials = 100;
  labs = zeros(1, num_trials);
  outs = zeros(num_trials, task_info.trial_dur);
  rrs = zeros(num_trials, temp.N, task_info.trial_dur);
  for trial = 1:num_trials
    [u, lab] = fnc_generate_trials('go-nogo2', task_info);
    labs(trial) = lab;
    out = fnc_eval_model(full_path, u, 'sigmoid', []);

    outs(trial, :) = out('O');
    rrs(trial, :, :) = out('X')';
  end
  rrs_A = squeeze(nanmean(rrs(labs == 0, curr_neus, stim_off:end), 3)); % 20%
  rrs_B = squeeze(nanmean(rrs(labs == 1, curr_neus, stim_off:end), 3)); % 80%
  tone_ids = [ones(1, size(rrs_A, 1)), ones(1, size(rrs_B, 1))*2];

  % Tone A is oddball
  task_info.prob = 0.20;
  num_trials = 100;
  labs = zeros(1, num_trials);
  outs = zeros(num_trials, task_info.trial_dur);
  rrs = zeros(num_trials, temp.N, task_info.trial_dur);
  for trial = 1:num_trials
    [u, lab] = fnc_generate_trials('go-nogo2', task_info);
    labs(trial) = lab;
    out = fnc_eval_model(full_path, u, 'sigmoid', []);

    outs(trial, :) = out('O');
    rrs2(trial, :, :) = out('X')';
  end
  rrs_A2 = squeeze(nanmean(rrs2(labs == 0, curr_neus, stim_off:end), 3)); % 20%
  rrs_B2 = squeeze(nanmean(rrs2(labs == 1, curr_neus, stim_off:end), 3)); % 80%
  tone_ids = [tone_ids, ones(1, size(rrs_A2, 1)), ones(1, size(rrs_B2, 1))*2];
  ctx_ids = [ones(1, num_trials), ones(1, num_trials)*2];

  rrs_combined = [rrs_A; rrs_B; rrs_A2; rrs_B2];

  for ii = 1:200 % for each neu
    T = table(tone_ids', ctx_ids', rrs_combined(:, ii));
    mdl = fitglm(T,'Var3~1+Var1*Var2','DummyVarCoding','effects');
    all_coeffs(i, ii, :) = mdl.Coefficients.Estimate(2:4)';
    all_ps(i, ii, :) = mdl.Coefficients.pValue(2:4)';
  end
end

init_coeffs = squeeze(all_coeffs(1, :, :));
init_ps = squeeze(all_ps(1, :, :));

last_coeffs = squeeze(all_coeffs(end, :, :));
last_ps = squeeze(all_ps(end, :, :));

init_ratios = init_coeffs(:, 2)./init_coeffs(:, 1);
last_ratios = last_coeffs(:, 2)./last_coeffs(:, 1);
combined_ratios = [init_ratios, last_ratios];

figure; hold on;
plot(combined_ratios', 'k');
plot(mean(combined_ratios), 'r', 'linewidth', 2)
set(gca, 'YScale', 'log')
xlim([0.5, 2.5])
ylabel('beta ratio (context/tone)')
xticks([1, 2])
xticklabels({'Early', 'Late'});
%print(fullfile(fig_dir, 'RNN_glm_results_all.eps'), '-painters', '-depsc');


diff_ratios = combined_ratios(:, 2) - combined_ratios(:, 1);
pos_ratios = find(diff_ratios > 0);

figure; hold on;
plot(combined_ratios(pos_ratios, :)', 'k');
plot(mean(combined_ratios(pos_ratios, :)), 'r', 'linewidth', 2)
xlim([0.5, 2.5])
ylabel('beta ratio (context/tone)')
xticks([1, 2])
xticklabels({'Early', 'Late'});
%print(fullfile(fig_dir, 'RNN_glm_results.eps'), '-painters', '-depsc');


sens_w_in = sum(abs(temp.w_in), 2);
sens_ratios = combined_ratios(sens_w_in < 1, :);
figure; hold on;
plot(sens_ratios', 'k');
set(gca, 'YScale', 'log')
plot(mean(sens_ratios), 'r', 'linewidth', 2)
xlim([0.5, 2.5])
ylabel('beta ratio (context/tone)')
xticks([1, 2])
xticklabels({'Early', 'Late'});



