% Name: Robert Kim, Nuttida Rungratsameetaweemana
% Date: 09-02-2024
% Email: rkim@salk.edu
% svm_decoding.m
% Description:

clear; clc;

rng(825)

fig_dir = fullfile(pwd, '..', '..', 'figures');
result_dir = fullfile(pwd, '..', '..', 'results');
model_path = fullfile(pwd, '..', '..', 'models', 'go-nogo2', 'P_rec_0.2_Taus_4.0_25.0');

model_ids = {'194212', '153216', '201555', '183934', '234618', '193127', '195849', '153730', '195329', '184637'}; % trained model IDs. Replace with the IDs that you trained

shuffle_cond = false;

svm_results = load(fullfile(result_dir, 'svm_results.mat'));
temp_all_data = svm_results.all_data;

num_neus = 200;
disp(['Number of neurons is set to ' num2str(num_neus)]);

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 200; % trial duration (timesteps)
task_info.stim_on = 50;
task_info.stim_dur = 25;
task_info.prob = 0.50;

stim_on = task_info.stim_on;
stim_dur = task_info.stim_dur;
stim_off = task_info.stim_on + task_info.stim_dur + 1;

all_tone_acc = [];
all_ctx_acc = [];
all_ntpc_ctx_acc = [];
all_ntpc_tone_acc = [];
all_trs = [];

all_data = struct(); % variable containing all the SVM decoding accuracy and p-values

for fi = 1:length(model_ids)
  fi
  model_id = model_ids{fi};
  curr_stage0_fnames = dir(fullfile(model_path, ['END_Stage_0*' model_id '*']));
  %curr_stage0_fnames = dir(fullfile(model_path, ['Stage_0_Tr_101_*' model_id '*']));
  curr_stage2_fnames = dir(fullfile(model_path, ['END_Stage_2*' model_id '*']));

  fnames = {curr_stage0_fnames.name, curr_stage2_fnames.name};

  fdr_left_ctx_ps = mafdr(temp_all_data(fi).left_ctx_ps, 'BHFDR', true);
  ctx_inc_units = find(fdr_left_ctx_ps < 0.05);

  fdr_right_tone_ps = mafdr(temp_all_data(fi).right_tone_ps, 'BHFDR', true);
  fdr_left_tone_ps = mafdr(temp_all_data(fi).left_tone_ps, 'BHFDR', true);
  tone_dec_units = find(fdr_right_tone_ps < 0.05);
  tone_inc_units = find(fdr_left_tone_ps < 0.05);

  if isempty(ctx_inc_units)
    disp('Skipping')
    continue;
  end

  %% make sure this is set correctly
  reps = 100; % repetition for SVM

  pref_accs = zeros(length(fnames), reps);
  pref_accs2 = zeros(length(fnames), reps);

  for ii = 1:length(fnames)
    curr_fi = fnames{ii};

    full_path = fullfile(model_path, curr_fi);
    load(full_path);

    if ii == 1
      all_trs = [all_trs, tr]
    end

    % Tone B is oddball
    task_info.prob = 0.80;
    num_trials = 100;
    labs = zeros(1, num_trials);
    outs = zeros(num_trials, task_info.trial_dur);
    rrs = zeros(num_trials, N, task_info.trial_dur);
    for trial = 1:num_trials
      [u, lab] = fnc_generate_trials('go-nogo2', task_info);
      %u = u + randn(size(u));
      labs(trial) = lab;
      out = fnc_eval_model(full_path, u, 'sigmoid', []);

      outs(trial, :) = out('O');
      rrs(trial, :, :) = out('R')'; % trials x neurons x time
    end

    % Tone A is oddball
    task_info.prob = 0.20;
    num_trials = 100;
    labs2 = zeros(1, num_trials);
    outs = zeros(num_trials, task_info.trial_dur);
    rrs2 = zeros(num_trials, N, task_info.trial_dur);
    for trial = 1:num_trials
      [u, lab] = fnc_generate_trials('go-nogo2', task_info);
      %u = u + randn(size(u));
      labs2(trial) = lab;
      out = fnc_eval_model(full_path, u, 'sigmoid', []);

      outs(trial, :) = out('O');
      rrs2(trial, :, :) = out('R')'; % trials x neurons x time
    end

    combined_data = squeeze(nanmean(rrs(:, :, :), 3));
    combined_data = squeeze(nanmean(combined_data, 2));
    combined_data2 = squeeze(nanmean(rrs2(:, :, :), 3));
    combined_data2 = squeeze(nanmean(combined_data2, 2));
    ctx_labs = [zeros(1, num_trials), ones(1, num_trials)];

    combined_data = [combined_data; combined_data2];
    combined_labs = [labs, labs2];

    for rp = 1:reps
      % First shuffle the data
      shuffle_rows = randperm(size(combined_data, 1));
      sh_combined_data = combined_data(shuffle_rows, :);
      sh_labs = combined_labs(shuffle_rows);
      sh_ctx_labs = ctx_labs(shuffle_rows);
      if shuffle_cond == true
        shuffle_rows2 = randperm(size(combined_data, 1));
        sh_labs = sh_labs(shuffle_rows2);
        sh_ctx_labs = sh_ctx_labs(shuffle_rows2);
      end

      train_n = floor(size(sh_combined_data, 1)*0.70);
      train_x = sh_combined_data(1:train_n, :);
      train_y = sh_labs(1:train_n);
      train_y2 = sh_ctx_labs(1:train_n);

      test_x = sh_combined_data(train_n+1:end, :);
      test_y = sh_labs(train_n+1:end)';
      test_y2 = sh_ctx_labs(train_n+1:end)';

      % Train an SVM classifier
      svmModel = fitcsvm(train_x, train_y, 'KernelFunction', 'linear');
      svmModel2 = fitcsvm(train_x, train_y2, 'KernelFunction', 'linear');

      % Predict and evaluate accuracy
      pred_y = predict(svmModel, test_x);
      pred_y2 = predict(svmModel2, test_x);
      accuracy = sum(pred_y == test_y) / length(test_y);
      accuracy2 = sum(pred_y2 == test_y2) / length(test_y2);
      pref_accs(ii, rp) = accuracy;
      pref_accs2(ii, rp) = accuracy2;
      end
  end

  all_data(fi).pref_accs = pref_accs;
  all_data(fi).pref_accs2 = pref_accs2;
  all_data(fi).inh = inh;

  tone_acc = squeeze(nanmean(pref_accs, 2)); % tr x neuron
  ctx_acc = squeeze(nanmean(pref_accs2, 2));


end

all_mean = [];
for i = 1:length(all_data)
  curr = mean(all_data(i).pref_accs2, 2);
  all_mean = [all_mean; curr'];
end
