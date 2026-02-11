% Name: Robert Kim
% Date: 3/4/2025

clear; clc;
rng(820)

fig_dir = fullfile(pwd, '..', '..', 'figures');
result_dir = fullfile(pwd, '..', '..', 'results');
model_path = fullfile(pwd, '..', '..', 'models', 'go-nogo2', 'P_rec_0.2_Taus_4.0_25.0');

% N = 200 models
%model_ids = {'152648', '153216', '153730', '224716', '225712', '234618', '232315', '232920', '234109', '233512'};
%model_ids = {'152648', '153216', '153730', '224716', '225712', '234618', '232315', '232920', '234109', '233512', '175820', '180438', '181032', '181637', '182143', '182720', '183343', '183934', '184637', '185315', '193043', '193127', '193633', '194212', '194125', '194819', '194910', '195329', '195624', '195849'};

model_ids = {'152648', '153216', '153730', '224716', '225712', '234618', '232315', '232920', '234109', '233512', '175820', '180438', '181032', '181637', '182143', '182720', '183343', '183934', '184637', '185315', '193043', '193127', '193633', '194212', '194125', '194819', '194910', '195329', '195624', '195849', '200232', '200522', '200806', '201033', '201342', '201555', '201905', '202049', '202426', '222224'};

% These are obtained from compare_models.m
%model_ids = {'153216', '153730', '183934', '185315', '193127', '194212', '195329', '195849', '201033', '201555'};

% SVM results (from svm_decoding_per_neuron.m)
%svm_results = load(fullfile(result_dir, 'svm_results_N10.mat'));
svm_results = load(fullfile(result_dir, 'svm_results.mat'));
tone_dec_ctx_inc_sig_N = svm_results.tone_dec_ctx_inc_N;
tone_dec_ctx_dec_sig_N = svm_results.tone_dec_ctx_dec_N;
tone_inc_ctx_dec_sig_N = svm_results.tone_inc_ctx_dec_N;
tone_inc_ctx_inc_sig_N = svm_results.tone_inc_ctx_inc_N;

many_ids = find(tone_dec_ctx_inc_sig_N > prctile(tone_dec_ctx_inc_sig_N, 75));
few_ids = setdiff([1:length(model_ids)], many_ids);

many_iis = [];
for i = 1:length(many_ids)
    first_model_id = model_ids{many_ids(i)};
    first_model_fname = dir(fullfile(model_path, ['END_Stage_0*' first_model_id '*']));
    second_model_fname = dir(fullfile(model_path, ['END_Stage_2*' first_model_id '*']));

    first_model = load(fullfile(model_path, first_model_fname.name));
    second_model = load(fullfile(model_path, second_model_fname.name));

    first_ww = first_model.w*first_model.m;
    second_ww = second_model.w*second_model.m;

    first_ii = first_ww(:, :);
    second_ii = second_ww(:, :);

    diff_ii = abs(second_ii - first_ii);

    %first_ii = first_ii(:); first_ii = first_ii(first_ii ~= 0);
    first_ii = diff_ii(:); diff_ii = diff_ii(diff_ii ~= 0);
    many_iis = [many_iis, mean(first_ii)];
end

few_iis = [];
for i = 1:length(few_ids)
    first_model_id = model_ids{few_ids(i)};
    first_model_fname = dir(fullfile(model_path, ['END_Stage_0*' first_model_id '*']));
    second_model_fname = dir(fullfile(model_path, ['END_Stage_2*' first_model_id '*']));

    first_model = load(fullfile(model_path, first_model_fname.name));
    second_model = load(fullfile(model_path, second_model_fname.name));

    first_ww = first_model.w*first_model.m;
    second_ww = second_model.w*second_model.m;

    first_ii = first_ww(:, :);
    second_ii = second_ww(:, :);

    diff_ii = abs(second_ii - first_ii);

    %first_ii = first_ii(:); first_ii = first_ii(first_ii ~= 0);
    first_ii = diff_ii(:); diff_ii = diff_ii(diff_ii ~= 0);
    few_iis = [few_iis, mean(first_ii)];
end

[p, h] = ranksum(many_iis, few_iis, 'tail', 'left')

all_iis = [many_iis, few_iis];
sig_N = svm_results.tone_dec_ctx_inc_N([many_ids, few_ids]);
[r, p] = corr(all_iis', sig_N', 'type', 'Pearson')

nonzero_ind = find(all_iis ~= 0);
[r, p] = corr(all_iis(nonzero_ind)', sig_N(nonzero_ind)', 'type', 'Pearson')

%{

figure;
scatter(all_iis, sig_N, 'filled');

%}