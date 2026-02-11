% Name: Robert Kim
% Date: 05-28-2023
% Email: robert.f.kim@gmail.com, nr2869@columbia.edu
% Description:

clear; clc;

model_path = <WHERE TRAINED RNNS ARE SAVED>
model_name = <TRAINED RNN MODEL FILENAME>
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


