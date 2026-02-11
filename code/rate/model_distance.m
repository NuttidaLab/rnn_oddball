% Name: Robert Kim
% Date: 05-28-2023
% Email: rkim@salk.edu
% eval_model.m
% Description:

clear; clc;
rng(821);

fig_dir = fullfile(pwd, '..', '..', 'figures', 'manuscript');
model_path = fullfile(pwd, '..', '..', 'models', 'go-nogo2', 'P_rec_0.2_Taus_4.0_25.0');
%model_id = '153216';
%model_id = '152648';
%model_id = '153730';
%model_ids = {'005148', '005732', '010345', '153216', '152648', '153730', '162159', '161647', '161047', '155832'};
%model_ids = {'152648', '153216', '153730', '224716', '225712', '234618', '232315', '232920', '234109', '233512'};
model_ids = {'194212', '153216', '201555', '183934', '234618', '193127', '195849', '153730', '195329', '184637'};


% First get the maximum number of trials across all the files
max_trs = 0;
for fi = 1:length(model_ids)
  model_id = model_ids{fi};
  all_fnames = dir(fullfile(model_path, ['Stage_*' model_id '*']));
  all_fnames = {all_fnames.name};
  trs = [1:500:length(all_fnames)*100];
  max_trs = max(max_trs, length(trs));
end

resample_N = 100;
vecDs = nan(length(model_ids), resample_N, max_trs);
angs = nan(length(model_ids), resample_N, max_trs);
num_neus = 10;

for fi = 1:length(model_ids)
  model_id = model_ids{fi};
  disp(model_id)
  all_fnames = dir(fullfile(model_path, ['Stage_*' model_id '*']));
  all_fnames = {all_fnames.name};
  trs = [1:500:length(all_fnames)*100];


  parfor ii = 1:resample_N
    ii

    % Test the trained model on randomly generated trials
    task_info = struct();
    task_info.trial_dur = 200; % trial duration (timesteps)
    task_info.stim_on = 50;
    task_info.stim_dur = 25;
    %task_info.prob = 0.50;
    stim_on = task_info.stim_on;
    stim_off = task_info.stim_on + task_info.stim_dur + 1;

    curr_neus = randsample(200, num_neus, 'false');
    temp_vec = zeros(1, length(trs));
    temp_ang = zeros(1, length(trs));
    for i = 1:length(trs)
      curr_fi = find(contains(all_fnames, ['Tr_' num2str(trs(i)) '_']));
      %curr_fi = ['Stage_0_Tr_' num2str(trs(i)) '_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_153216.mat'];

      full_path = fullfile(model_path, all_fnames{curr_fi});
      temp = load(full_path);

      % Tone A is oddball
      task_info.prob = 0.20
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
      rrs_A = squeeze(nanmean(rrs(labs == 0, curr_neus, stim_off:end))); % 20%
      rrs_B = squeeze(nanmean(rrs(labs == 1, curr_neus, stim_off:end))); % 80%
      %rrs_A = squeeze(nanmean(rrs(labs == 0, curr_neus, stim_on:stim_off))); % 20%
      %rrs_B = squeeze(nanmean(rrs(labs == 1, curr_neus, stim_on:stim_off))); % 80%

      rrs_A = squeeze(nanmean(rrs_A, 2));
      rrs_B = squeeze(nanmean(rrs_B, 2));

      %rrs_A = rrs_A'/norm(rrs_A);
      %rrs_B = rrs_B'/norm(rrs_B);
      rrs_A = rrs_A';
      rrs_B = rrs_B';

      vecD = pdist2(rrs_A, rrs_B, "euclidean");
      ang = Vangle(rrs_A, rrs_B);
      temp_vec(i) = vecD;
      temp_ang(i) = ang;
    end
    temp_nan = nan(1, max_trs);
    temp_nan(1:length(temp_vec)) = temp_vec;
    vecDs(fi, ii, :) = temp_nan;

    temp_nan = nan(1, max_trs);
    temp_nan(1:length(temp_ang)) = temp_ang;
    angs(fi, ii, :) = temp_nan;
  end
end

meanVec1 = squeeze(nanmean(vecDs(:, :, :), 2));
meanVec = nanmean(meanVec1);
semVec = sem(meanVec1);

meanAng1 = squeeze(nanmean(angs(:, :, :), 2));
meanAng = nanmean(meanAng1);
semAng = sem(meanAng1);

%save(fullfile(fig_dir, 'dist_angle_nonnorm_results_v2.mat'), 'vecDs', 'angs', 'model_ids', 'resample_N', 'num_neus', 'meanVec', 'semVec', 'meanAng', 'semAng');

%{

figure; hold on;
yyaxis left

% Create x values for the plot
x = 1:length(meanVec);

% Plot the shaded area for SEM
fill([x, fliplr(x)], [meanVec + semVec, fliplr(meanVec - semVec)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % Adjust FaceAlpha for transparency

% Plot the mean line on top
plot(x, meanVec, 'ro-', 'MarkerFaceColor', 'r');
ylabel('distance (a.u.)')

yyaxis right
% Create x values for the plot
x = 1:length(meanAng);

% Plot the shaded area for SEM
fill([x, fliplr(x)], [meanAng + semAng, fliplr(meanAng - semAng)], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none');  % Adjust FaceAlpha for transparency

% Plot the mean line on top
plot(x, meanAng, 'bo-', 'MarkerFaceColor', 'b');
xlim([1, 6.2])
ylabel('cosine angle')
xlabel('trials (x500)')
xticks([1 2 3 4 5 6]);
xticklabels({'0', '1', '2', '3', '4', '5'});

%print(fullfile(fig_dir, 'RNN_dist_angle_across_RNNs_non_norm_v2.eps'), '-painters', '-depsc');

%errorbar(nanmean(vecDs), sem(vecDs), 'ro-')
%}








