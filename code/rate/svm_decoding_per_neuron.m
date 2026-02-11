% Name: Robert Kim, Nuttida Rungratsameetaweemana
% Date: 09-02-2024
% Email: robert.f.kim@gmail.com, nr2869@columbia.edu
% svm_decoding.m

clear; clc;

% Close any existing pool
delete(gcp('nocreate'));

% Start a local pool with multiple *process* workers (not threads):
parpool('Processes');  % or parpool('local', 'SpmdEnabled', true)

% Now you can run:
pctRunOnAll('rng(8330)');


lesion_type = 'ie';
disp(['Lesioning ' lesion_type])

fig_dir = fullfile(pwd, '..', '..', 'figures');
result_dir = fullfile(pwd, '..', '..', 'results');
model_path = fullfile(pwd, '..', '..', 'models', 'go-nogo2', 'P_rec_0.2_Taus_4.0_25.0');

model_ids = {'194212', '153216', '201555', '183934', '234618', '193127', '195849', '153730', '195329', '184637'}; % trained model IDs. Replace with the IDs that you trained

shuffle_cond = false;

num_neus = 200;
disp(['Number of neurons is set to ' num2str(num_neus)]);

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 200; % trial duration (timesteps)
task_info.stim_on = 50;
task_info.stim_dur = 25;
task_info.prob = 0.50;

stim_on = task_info.stim_on;
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

  %% make sure this is set correctly
  reps = 100; % repetition for SVM

  pref_accs = zeros(length(fnames), num_neus, reps);
  pref_accs2 = zeros(length(fnames), num_neus, reps);

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
      out = fnc_eval_model(full_path, u, 'sigmoid', [], lesion_type);

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
      out = fnc_eval_model(full_path, u, 'sigmoid', [], lesion_type);

      outs(trial, :) = out('O');
      rrs2(trial, :, :) = out('R')'; % trials x neurons x time
    end

    parfor neu = 1:double(N) % for each neuron
      combined_data = squeeze(nanmean(rrs(:, neu, :), 3));
      combined_data2 = squeeze(mean(rrs2(:, neu, :), 3));
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
        pref_accs(ii, neu, rp) = accuracy;
        pref_accs2(ii, neu, rp) = accuracy2;
      end
    end
  end

  all_data(fi).pref_accs = pref_accs;
  all_data(fi).pref_accs2 = pref_accs2;
  all_data(fi).inh = inh;

  tone_acc = squeeze(nanmean(pref_accs, 3)); % tr x neuron
  ctx_acc = squeeze(nanmean(pref_accs2, 3));

  % Find the neurons with tone difference
  left_tone_ps = zeros(1, double(N));
  right_tone_ps = zeros(1, double(N));
  all_tone_decoding = zeros(double(N), 2);
  for neu = 1:double(N)
    [p, h] = ranksum(squeeze(pref_accs(1, neu, :)), squeeze(pref_accs(2, neu, :)), 'tail', 'left');
    left_tone_ps(neu) = p;

    [p, h] = ranksum(squeeze(pref_accs(1, neu, :)), squeeze(pref_accs(2, neu, :)), 'tail', 'right');
    right_tone_ps(neu) = p;

    all_tone_decoding(neu, :) = [squeeze(mean(pref_accs(1, neu, :), 3)), squeeze(mean(pref_accs(2, neu, :), 3))];
  end

  all_data(fi).left_tone_ps = left_tone_ps;
  all_data(fi).right_tone_ps = right_tone_ps;
  all_data(fi).all_tone_decoding = all_tone_decoding;

  % Find the neurons with ctx difference
  left_ctx_ps = zeros(1, double(N));
  right_ctx_ps = zeros(1, double(N));
  all_ctx_decoding = zeros(double(N), 2);
  for neu = 1:double(N)
    [p, h] = ranksum(squeeze(pref_accs2(1, neu, :)), squeeze(pref_accs2(2, neu, :)), 'tail', 'left');
    left_ctx_ps(neu) = p;

    [p, h] = ranksum(squeeze(pref_accs2(1, neu, :)), squeeze(pref_accs2(2, neu, :)), 'tail', 'right');
    right_ctx_ps(neu) = p;

    all_ctx_decoding(neu, :) = [squeeze(mean(pref_accs2(1, neu, :), 3)), squeeze(mean(pref_accs2(2, neu, :), 3))];
  end
  all_data(fi).left_ctx_ps = left_ctx_ps;
  all_data(fi).right_ctx_ps = right_ctx_ps;
  all_data(fi).all_ctx_decoding = all_ctx_decoding;

  ov_tone_acc = tone_acc(2, :);
  ov_ctx_acc = ctx_acc(2, :);

  all_tone_acc = [all_tone_acc, mean(ov_tone_acc)];
  all_ctx_acc = [all_ctx_acc, mean(ov_ctx_acc)];

  % Tone neurons
  tone_neus = tone_acc(2, :) - tone_acc(1, :); % only use the first two time points for now
  pos_tone_neus = find(tone_neus > 0);
  neg_tone_neus = find(tone_neus < 0);

  % CTX neurons
  ctx_neus = ctx_acc(2, :) - ctx_acc(1, :); % only use the first two time points for now
  pos_ctx_neus = find(ctx_neus > 0);
  neg_ctx_neus = find(ctx_neus < 0);

  % CTX and Tone interaction
  pos_tone_neg_ctx = intersect(pos_tone_neus, neg_ctx_neus);
  neg_tone_pos_ctx = intersect(neg_tone_neus, pos_ctx_neus);

  ntpc_tone = mean(tone_acc(:, neg_tone_pos_ctx), 2);
  ntpc_ctx = mean(ctx_acc(:, neg_tone_pos_ctx), 2);

  all_ntpc_tone_acc = [all_ntpc_tone_acc, ntpc_tone];
  all_ntpc_ctx_acc = [all_ntpc_ctx_acc, ntpc_ctx];

end

tone_inc_N = [];
tone_dec_N = [];
ctx_inc_N = [];
ctx_dec_N = [];
unchanged_N = [];
tone_dec_ctx_inc_N = [];
tone_dec_ctx_dec_N = [];
tone_inc_ctx_inc_N = [];
tone_inc_ctx_dec_N = [];
tone_dec_ctx_inc_inh = [];
for i = 1:length(all_data)
  fdr_right_tone_ps = mafdr(all_data(i).right_tone_ps, 'BHFDR', true);
  fdr_left_tone_ps = mafdr(all_data(i).left_tone_ps, 'BHFDR', true);
  tone_dec_units = find(fdr_right_tone_ps < 0.05);
  tone_inc_units = find(fdr_left_tone_ps < 0.05);

  fdr_left_ctx_ps = mafdr(all_data(i).left_ctx_ps, 'BHFDR', true);
  fdr_right_ctx_ps = mafdr(all_data(i).right_ctx_ps, 'BHFDR', true);
  ctx_inc_units = find(fdr_left_ctx_ps < 0.05);
  ctx_dec_units = find(fdr_right_ctx_ps < 0.05);

  temp_uni = union(tone_dec_units, tone_inc_units);
  temp_uni = union(temp_uni, ctx_inc_units);
  temp_uni = union(temp_uni, ctx_dec_units);
  unchanged_units = setdiff([1:num_neus], temp_uni);
  unchanged_N = [unchanged_N, length(unchanged_units)];

  tone_inc_N = [tone_inc_N, length(tone_inc_units)];
  tone_dec_N = [tone_dec_N, length(tone_dec_units)];
  ctx_inc_N = [ctx_inc_N, length(ctx_inc_units)];
  ctx_dec_N = [ctx_dec_N, length(ctx_dec_units)];

  % Neurons with decreased tone decoding over time and with increased ctx decoding over time
  tone_dec_ctx_inc = intersect(tone_dec_units, ctx_inc_units);
  tone_dec_ctx_dec = intersect(tone_dec_units, ctx_dec_units);
  tone_inc_ctx_inc = intersect(tone_inc_units, ctx_inc_units);
  tone_inc_ctx_dec = intersect(tone_inc_units, ctx_dec_units);

  if ~isempty(tone_dec_ctx_inc)
    tone_dec_ctx_inc_inh = [tone_dec_ctx_inc_inh, sum(all_data(i).inh(tone_dec_ctx_inc))/length(tone_dec_ctx_inc)];
  else
    tone_dec_ctx_inc_inh = [tone_dec_ctx_inc_inh, nan];
  end
  all_data(i).tone_dec_ctx_inc_inh = all_data(i).inh(tone_dec_ctx_inc);

  tone_accs = all_data(i).pref_accs(:, tone_dec_ctx_inc);
  ctx_accs = all_data(i).pref_accs2(:, tone_dec_ctx_inc);

  tone_dec_ctx_inc_N = [tone_dec_ctx_inc_N, length(tone_dec_ctx_inc)];
  tone_dec_ctx_dec_N = [tone_dec_ctx_dec_N, length(tone_dec_ctx_dec)];
  tone_inc_ctx_inc_N = [tone_inc_ctx_inc_N, length(tone_inc_ctx_inc)];
  tone_inc_ctx_dec_N = [tone_inc_ctx_dec_N, length(tone_inc_ctx_dec)];
end

if shuffle_cond == false
  save(fullfile(result_dir, [lesion_type '_svm_results_v2.mat']), 'all_data', 'model_ids', 'all_trs', "tone_dec_ctx_inc_N", "tone_dec_ctx_dec_N", "tone_inc_ctx_dec_N", "tone_inc_ctx_inc_N", "tone_dec_ctx_inc_inh", "tone_inc_N", "unchanged_N", "tone_dec_N", "ctx_inc_N", "ctx_dec_N", "all_tone_acc", "all_ctx_acc", "all_ntpc_tone_acc", "all_ntpc_ctx_acc");
else
  save(fullfile(result_dir, 'shuffled_svm_results_v2.mat'), 'all_data', 'model_ids', 'all_trs', "tone_dec_ctx_inc_N", "tone_dec_ctx_dec_N", "tone_inc_ctx_dec_N", "tone_inc_ctx_inc_N", "tone_dec_ctx_inc_inh", "tone_inc_N", "unchanged_N", "tone_dec_N", "ctx_inc_N", "ctx_dec_N", "all_tone_acc", "all_ctx_acc", "all_ntpc_tone_acc", "all_ntpc_ctx_acc");
end

%% Boxplot of overall accuracy (tone and context)
figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
hold on;
boxplot([all_tone_acc; all_ctx_acc]', 'symbol', '');
ylim([0.3 0.8]);
xlim([0.5, 2.5]);
plot(xlim, [0.5, 0.5], 'k--');
%print(fullfile(fig_dir, 'overall_decoding_boxplot.eps'), '-painters', '-depsc');
pause;


%{
figure; hold on;
h1 = plot(tone_acc(:, pos_tone_neus)*100, 'r');
h2 = plot(tone_acc(:, neg_tone_neus)*100, 'b');
xlim([0.5, 2.5])
ylim([40 100])
xticks([1 2]);
xticklabels({'early', 'late'});
legend([h1(1), h2(1)], {'"Positive" Tone', '"Negative" Tone'}, 'Location', 'northwest');
ylabel('SVM Decoding Performance (%)');
title('Tone Identity Decodability');
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');

figure; hold on;
h1 = plot(ctx_acc(:, pos_ctx_neus)*100, 'r');
h2 = plot(ctx_acc(:, neg_ctx_neus)*100, 'b');
xlim([0.5, 2.5])
ylim([40 100])
xticks([1 2]);
xticklabels({'early', 'late'});
legend([h1(1), h2(1)], {'"Positive" Context', '"Negative" Context'}, 'Location', 'northwest');
ylabel('SVM Decoding Performance (%)')
title('Context Decodability');
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
%}


figure; hold on;
h1 = plot(all_ntpc_ctx_acc, 'r');
h2 = plot(all_ntpc_tone_acc, 'b');
xlim([0.5, 2.5])
ylim([0.4 0.6])
xticks([1 2]);
xticklabels({'early', 'late'});
legend([h1(1), h2(1)], {'"Positive" Context', '"Negative" Tone'}, 'Location', 'northwest');
ylabel('SVM Decoding Performance (%)')
title('Neurons that switch encoding');
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
%print(fullfile(fig_dir, 'code_switching_neus.eps'), '-painters', '-depsc');
 

figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
hold on;
boxplot(all_ntpc_ctx_acc', 'symbol', '');
ylim([0.4 0.70])
xlim([0.5, 2.5])
plot(xlim, [0.5, 0.5], 'k--');
[p, h] = signrank(all_ntpc_ctx_acc(1, :), all_ntpc_ctx_acc(2, :))
%print(fullfile(fig_dir, 'ctx_decoding_boxplot.eps'), '-painters', '-depsc');

figure('Units', 'Inch', 'Outerposition', [0 0 4 4]);
hold on;
boxplot(all_ntpc_tone_acc', 'symbol', '');
ylim([0.40 0.85])
xlim([0.5, 2.5])
plot(xlim, [0.5, 0.5], 'k--');
[p, h] = signrank(all_ntpc_tone_acc(1, :), all_ntpc_tone_acc(2, :))
%print(fullfile(fig_dir, 'tone_decoding_boxplot.eps'), '-painters', '-depsc');

%{
%% CTX and Tone difference
figure; hold on;
plot((ctx_acc(:, neg_tone_pos_ctx) - tone_acc(:, neg_tone_pos_ctx))*100, 'r');
xlim([0.5, 2.5])
xticks([1 2]);
xticklabels({'early', 'late'});
ylabel('SVM Decoding Performance (%)')
title('Difference in SVM decodability (Context - Tone)')
set(gcf, 'Color', 'w');
set(gca, 'Color', 'w');
% print(fullfile(fig_dir, 'code_switching_neus.eps'), '-painters', '-depsc');
%}


%{
aa = squeeze(nanmean(pref_accs, 2));   
bb = squeeze(nanmean(pref_accs2, 2));   

%figure; 
hold on;
plot(aa, 'ro-'); % stimulus decoding
plot(bb, 'bo-'); % oddball/statistics decoding
ylim([0.4 1]); 
plot(mean(aa, 2), 'ko-', 'linewidth', 2); hold on; 
plot(mean(bb, 2), 'ko-', 'linewidth', 2) 

%}


