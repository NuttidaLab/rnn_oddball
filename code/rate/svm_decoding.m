% Name: Robert Kim
% Date: 09-02-2024
% Email: rkim@salk.edu
% svm_decoding.m
% Description:
%
% INPUT
%
%
% OUTPUT

clear; clc;

fig_dir = '/home/shared/oddball/figures';

model_path = '/home/shared/oddball/models/go-nogo2/P_rec_0.2_Taus_4.0_25.0';
%model_id = '153216';
%model_id = '153730';
model_id = '152648';

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 200; % trial duration (timesteps)
task_info.stim_on = 50;
task_info.stim_dur = 25;
task_info.prob = 0.50;

stim_on = task_info.stim_on;
stim_off = task_info.stim_on + task_info.stim_dur + 1;

%trs = [101:100:1801];
trs = [101:100:1701];

%% make sure this is set correctly
nboot = 100;
reps = 100; % repetition for SVM
rand_neus = 10; % random subsampling neurons

pref_accs = zeros(length(trs), nboot, reps);
pref_accs2 = zeros(length(trs), nboot, reps);

for ii = 1:length(trs)
  ii
  %curr_fi = ['Stage_0_Tr_' num2str(trs(ii)) '_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_153216.mat']; % 80 - 20
  curr_fi = ['Stage_0_Tr_' num2str(trs(ii)) '_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat']; % 80 - 20

  full_path = fullfile(model_path, curr_fi);
  load(full_path);

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

  for nb = 1:nboot
    rand_sub = randperm(200);
    combined_data = squeeze(mean(rrs(:, rand_sub(1:rand_neus), :), 3));
    combined_data2 = squeeze(mean(rrs2(:, rand_sub(1:rand_neus), :), 3));
    ctx_labs = [zeros(1, num_trials), ones(1, num_trials)];

    combined_data = [combined_data; combined_data2];
    combined_labs = [labs, labs2];

    parfor rp = 1:reps
      % First shuffle the data
      shuffle_rows = randperm(size(combined_data, 1));
      sh_combined_data = combined_data(shuffle_rows, :);
      sh_labs = combined_labs(shuffle_rows);
      sh_ctx_labs = ctx_labs(shuffle_rows);

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
      pref_accs(ii, nb, rp) = accuracy;
      pref_accs2(ii, nb, rp) = accuracy2;
    end
  end
end

aa = squeeze(nanmean(pref_accs, 2));   
bb = squeeze(nanmean(pref_accs2, 2));   

%figure; 
hold on;
plot(aa, 'ro-'); % stimulus decoding
plot(bb, 'bo-'); % oddball/statistics decoding
ylim([0.4 1]); 
plot(mean(aa, 2), 'ko-', 'linewidth', 2);
plot(mean(bb, 2), 'ko-', 'linewidth', 2) 

%print(fullfile(fig_dir, 'svm_decoding.eps'), '-painters', '-depsc');


