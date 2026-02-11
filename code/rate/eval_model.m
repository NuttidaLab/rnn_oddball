% Name: Robert Kim
% Date: 05-28-2023
% Email: rkim@salk.edu
% eval_model.m
% Description:

clear; clc;

model_path = '../../models/go-nogo/P_rec_0.2_Taus_4.0_25.0'; 
mat_file = dir(fullfile(model_path, '*.mat'));
model_name = mat_file(1).name;
model_path = fullfile(model_path, model_name);
load(model_path);

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 450; % trial duration (timesteps)
task_info.stim_on = 300;
task_info.stim_dur = 25;

num_trials = 500;
outs = zeros(num_trials, task_info.trial_dur);
rrs = zeros(num_trials, N, task_info.trial_dur);
gos = [];
go_rrs = zeros(N, task_info.trial_dur);
nogos = [];
nogo_rrs = zeros(N, task_info.trial_dur);
for trial = 1:num_trials
  u = zeros(1, task_info.trial_dur) + randn(1, task_info.trial_dur);
  u(task_info.stim_on:task_info.stim_on+task_info.stim_dur) = 0.4;
  out = fnc_eval_model(model_path, u, 'sigmoid', []);
  curr_out = out('O');
  if max(curr_out(task_info.stim_on:end)) > 0.7 & max(curr_out(1:task_info.stim_on)) < 0.7
    gos = [gos; curr_out];
    go_rrs = go_rrs + out('R')';
  elseif max(curr_out) < 0.7;
    nogos = [nogos; curr_out];
    nogo_rrs = nogo_rrs + out('R')';
  end
  outs(trial, :) = curr_out;
  rrs(trial, :, :) = out('X')';
end
go_rrs = go_rrs/size(gos, 1);
for i = 1:N
  curr_rr = go_rrs(i, :);
  baseline = mean(curr_rr(30:100));
  curr_rr = curr_rr - baseline;
  go_rrs(i, :) = curr_rr;
end

nogo_rrs = nogo_rrs/size(nogos, 1);
for i = 1:N
  curr_rr = nogo_rrs(i, :);
  baseline = mean(curr_rr(30:100));
  curr_rr = curr_rr - baseline;
  nogo_rrs(i, :) = curr_rr;
end

hold on;
plot(mean(go_rrs), 'r');
plot(mean(nogo_rrs), 'b');
ylim([-2.5, -2])




