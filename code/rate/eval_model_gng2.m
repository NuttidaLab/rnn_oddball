% Name: Robert Kim
% Date: 05-28-2023
% Email: rkim@salk.edu
% eval_model.m
% Description:

clear; clc;

model_path = '/home/shared/oddball/models/go-nogo2/P_rec_0.2_Taus_4.0_25.0';
%model_name = 'END_Stage_0_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_153216.mat';
%model_name = 'END_Stage_1_Prob_0.5_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_153216.mat';
model_name = 'END_Stage_2_Prob_0.2_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_153216.mat';

%model_name = 'END_Stage_0_Prob_0.5_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_21_010345_2stages.mat';
%model_name = 'END_Stage_1_Prob_0.2_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_21_010345_2stages.mat';

%mat_file = dir(fullfile(model_path, '*.mat'));
%model_name = mat_file(1).name;
model_path = fullfile(model_path, model_name);
load(model_path);

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 200; % trial duration (timesteps)
task_info.stim_on = 50;
task_info.stim_dur = 25;
task_info.prob = 0.50;

num_trials = 100;
labs = zeros(1, num_trials);
outs = zeros(num_trials, task_info.trial_dur);
rrs = zeros(num_trials, N, task_info.trial_dur);
for trial = 1:num_trials
  [u, lab] = fnc_generate_trials('go-nogo2', task_info);
  labs(trial) = lab;
  out = fnc_eval_model(model_path, u, 'sigmoid', []);

  outs(trial, :) = out('O');
  rrs(trial, :, :) = out('X')';
end

hold on;
plot(outs(labs== 0, :)', 'r');
plot(outs(labs== 1, :)', 'b');
plot([task_info.stim_on, task_info.stim_on], ylim, 'k--');
plot([task_info.stim_on + task_info.stim_dur, task_info.stim_on + task_info.stim_dur], ylim, 'k--');


