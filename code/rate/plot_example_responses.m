% Name: Robert Kim
% Date: 09-02-2024
% Email: rkim@salk.edu
% svm_decoding.m
% Description:

clear; clc;

fig_dir = '/home/shared/oddball/figures';

model_path = '/home/shared/oddball/models/go-nogo2/P_rec_0.2_Taus_4.0_25.0';
model_id = '152648';

% Test the trained model on randomly generated trials
task_info = struct();
task_info.trial_dur = 200; % trial duration (timesteps)
task_info.stim_on = 50;
task_info.stim_dur = 25;
task_info.prob = 0.50;

stim_on = task_info.stim_on;
stim_off = task_info.stim_on + task_info.stim_dur + 1;

% fnames = {'END_Stage_0_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat', 'Stage_0_Tr_101_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat', 'END_Stage_1_Prob_0.5_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat', 'END_Stage_2_Prob_0.2_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat'};

fnames = {'END_Stage_0_Prob_0.8_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat', 'END_Stage_1_Prob_0.5_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat', 'END_Stage_2_Prob_0.2_Task_go-nogo2_N_200_Taus_4.0_25.0_Act_sigmoid_2024_08_20_152648.mat'};
time_vec = [1:199]*5/1000;
for i = 1:length(fnames)
    load(fullfile(model_path, fnames{i}));
    figure; hold on;
    plot(time_vec, eval_os(eval_labels==0, :)', 'r')
    plot(time_vec, eval_os(eval_labels==1, :)', 'b')
    if i == 1
        fname = 'Stage0_output';
    elseif i == 2
        fname = 'Stage1_output';
    elseif i == 3
        fname = 'Stage2_output';
    end
    print(fullfile(fig_dir, [fname '.eps']), '-painters', '-depsc');
    pause(2)
    close;
end