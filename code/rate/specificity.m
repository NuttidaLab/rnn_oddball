% Name: Robert Kim, Nuttida Rungratsameetaweemana
% Date: 11-18-2024
% Email: robert.f.kim@gmail.com, nr2869@columbia.edu
% specificity.m
% Description: To find units with selectivity

clear; clc;
rng(821);
model_path = '/home/shared/oddball/models/go-nogo2/P_rec_0.2_Taus_4.0_25.0';
model_ids = {'005148', '005732', '010345', '153216', '152648', '153730', '162159', '161647', '161047', '155832'};

all_ps = zeros(length(model_ids), 3);
stim_ps = zeros(length(model_ids), 3);
ctx_ps = zeros(length(model_ids), 3);

num_trials = 100;

for ii = 1:length(model_ids)
  model_id = model_ids{ii};
  fname = dir(fullfile(model_path, ['END_Stage_1_*' model_id '*']));
  fname = fname.name;

  % Task structure
  task_info = struct();
  task_info.trial_dur = 200; % trial duration (timesteps)
  task_info.stim_on = 50;
  task_info.stim_dur = 25;
  %task_info.prob = 0.50;
  stim_on = task_info.stim_on;
  stim_off = task_info.stim_on + task_info.stim_dur + 1;

  full_path = fullfile(model_path, fname);
  temp = load(full_path);

  % Tone B is oddball
  task_info.prob = 0.80;
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
  rrs1 = rrs;
  labs1 = labs;

  % Tone A is oddball
  task_info.prob = 0.20;
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
  rrs2 = rrs;
  labs2 = labs;

  combined_rrs = [rrs1; rrs2];
  ctx_ids = [ones(1, num_trials), ones(1, num_trials)*2];

  % Any response (baseline vs. stim window)
  pre_rrs = combined_rrs(:, :, 1:stim_on);
  pre_rrs = squeeze(nanmean(pre_rrs, 3));

  stim_rrs = combined_rrs(:, :, (stim_on+5):(stim_on+task_info.stim_dur));
  stim_rrs = squeeze(nanmean(stim_rrs, 3));

  ps = zeros(1, size(pre_rrs, 2));
  for n = 1:size(pre_rrs, 2) % for each neuron
    [p, h] = signrank(pre_rrs(:, n), stim_rrs(:, n));
    ps(n) = p;
  end

  exc_ps = ps(find(temp.exc == 1));
  inh_ps = ps(find(temp.exc == 0));

  all_ps(ii, 1) = length(find(exc_ps < 0.001))/length(find(temp.exc == 1));
  all_ps(ii, 2) = length(find(inh_ps < 0.001))/length(find(temp.exc == 0));
  all_ps(ii, 3) = length(find(ps < 0.001))/200;

  % Stim 1 vs stim 2
  stim1_rrs = [rrs1(labs1 == 0, :, :); rrs2(labs2 == 0, :, :)];
  stim2_rrs = [rrs1(labs1 == 1, :, :); rrs2(labs2 == 1, :, :)];

  stim1_rrs = squeeze(nanmean(stim1_rrs(:, :, (stim_on+5):(stim_on+task_info.stim_dur)), 3));
  stim2_rrs = squeeze(nanmean(stim2_rrs(:, :, (stim_on+5):(stim_on+task_info.stim_dur)), 3));

  ps = zeros(1, size(stim1_rrs, 2));
  for n = 1:size(stim1_rrs, 2) % for each neuron
    [p, h] = ranksum(stim1_rrs(:, n), stim2_rrs(:, n));
    ps(n) = p;
  end

  exc_ps = ps(find(temp.exc == 1));
  inh_ps = ps(find(temp.exc == 0));
  stim_ps(ii, 1) = length(find(exc_ps < 0.001))/length(find(temp.exc == 1));
  stim_ps(ii, 2) = length(find(inh_ps < 0.001))/length(find(temp.exc == 0));
  stim_ps(ii, 3) = length(find(ps < 0.001))/200;

  % Context
  ctx1_rrs = rrs1;
  ctx2_rrs = rrs2;

  ctx1_rrs = squeeze(nanmean(ctx1_rrs(:, :, (stim_on+5):(stim_on+task_info.stim_dur)), 3));
  ctx2_rrs = squeeze(nanmean(ctx2_rrs(:, :, (stim_on+5):(stim_on+task_info.stim_dur)), 3));

  ps = zeros(1, size(ctx1_rrs, 2));
  for n = 1:size(ctx1_rrs, 2) % for each neuron
    [p, h] = ranksum(ctx1_rrs(:, n), ctx2_rrs(:, n));
    ps(n) = p;
  end

  exc_ps = ps(find(temp.exc == 1));
  inh_ps = ps(find(temp.exc == 0));
  ctx_ps(ii, 1) = length(find(exc_ps < 0.001))/length(find(temp.exc == 1));
  ctx_ps(ii, 2) = length(find(inh_ps < 0.001))/length(find(temp.exc == 0));
  ctx_ps(ii, 3) = length(find(ps < 0.001))/200;
end

%% Chi square test
% Data: proportions for two groups
group1 = ctx_ps(:, 1)';
group2 = ctx_ps(:, 2)';

% Assuming 100 trials per subject
n_trials = 100;
successes_group1 = round(group1 * n_trials);
failures_group1 = n_trials - successes_group1;

successes_group2 = round(group2 * n_trials);
failures_group2 = n_trials - successes_group2;

% Combine into a contingency table
observed = [
    sum(successes_group1), sum(successes_group2);
    sum(failures_group1), sum(failures_group2)
];

% Compute row and column totals
row_totals = sum(observed, 2);
col_totals = sum(observed, 1);
total = sum(observed(:));

% Compute expected values
expected = (row_totals * col_totals) / total;

% Compute the chi-square statistic
chi2_stat = sum((observed - expected).^2 ./ expected, 'all');

% Degrees of freedom
df = (size(observed, 1) - 1) * (size(observed, 2) - 1);

% Compute p-value
p = 1 - chi2cdf(chi2_stat, df);

% Display results
disp('Chi-Square Test Results:');
disp(['Chi-Square Statistic: ', num2str(chi2_stat)]);
disp(['Degrees of Freedom: ', num2str(df)]);
disp(['P-value: ', num2str(p)]);



