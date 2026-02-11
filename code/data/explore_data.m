% Name: Robert Kim
% Date: 02-01-2024
% Email: rkim@salk.edu
% explore_data.m
% Description: Script to explore the pixel data

clear; clc;

fig_dir = '/home/shared/oddball/figures';

data_dir = '/home/shared/oddball/data';
data_fname = 'spikeData05_oddball.mat';
load(fullfile(data_dir, data_fname));

%% COLORMAP
hex = ['#06D2AC'; '#206975'; '#6F3AA4'; '#2B1644'];
vec = [100, 50, 25, 0];
raw = sscanf(hex','#%2x%2x%2x',[3,size(hex,1)]).' / 255;
N = 128;
NR_RK_map = interp1(vec,raw,linspace(100,0,N),'pchip');

num_neus = length(goodData);

% Window in ms
window_st = -100;
window_end = 1000;

% times are in ms
all_bins = nan(num_neus, length(astamps));
for n = 1:num_neus
  curr_neu = double(goodData(n).spikeTimes);
  curr_bin = nan(1, length(astamps));
  for i = 1:length(astamps) % for each auditory onset
    curr_onset = astamps(i);
    curr_win = [curr_onset+window_st, curr_onset+window_end];
    curr_bin(i) = length(find(curr_neu > curr_win(1) & curr_neu < curr_win(2)));
  end
  all_bins(n, :) = curr_bin;
end

% "Go" is 200Hz tones (freqLabel = 1)
% "NoGo" is 5kHz tones (freqLabel = 2)
up_ps_gos = zeros(1, num_neus);
down_ps_gos = zeros(1, num_neus);
for n = 1:num_neus
  [p, h] = ranksum(all_bins(n, taskLabel == 1 & freqLabel==1), all_bins(n, taskLabel == 3 & freqLabel == 1), 'tail', 'left');
  up_ps_gos(n) = p;

  [p, h] = ranksum(all_bins(n, taskLabel == 1 & freqLabel==1), all_bins(n, taskLabel == 3 & freqLabel == 1), 'tail', 'right');
  down_ps_gos(n) = p;
end
sig_up_ps_gos = find(up_ps_gos < 0.001);
sig_down_ps_gos = find(down_ps_gos < 0.001);

up_ps_nogos = zeros(1, num_neus);
down_ps_nogos = zeros(1, num_neus);
for n = 1:num_neus
  [p, h] = ranksum(all_bins(n, taskLabel == 1 & freqLabel==2), all_bins(n, taskLabel == 3 & freqLabel == 2), 'tail', 'left');
  up_ps_nogos(n) = p;

  [p, h] = ranksum(all_bins(n, taskLabel == 1 & freqLabel==2), all_bins(n, taskLabel == 3 & freqLabel == 2), 'tail', 'right');
  down_ps_nogos(n) = p;
end
sig_up_ps_nogos = find(up_ps_nogos < 0.001);
sig_down_ps_nogos = find(down_ps_nogos < 0.001);

% Comparing channel depth
%nsig_nogo_ps = find(nogo_ps > 0.001);
%sig_nogo_CD = [];
%for i = 1:length(sig_nogo_ps)
  %sig_nogo_CD = [sig_nogo_CD, goodData(sig_nogo_ps(i)).chanDepth];
%end
%nsig_nogo_CD = [];
%for i = 1:length(nsig_nogo_ps)
  %nsig_nogo_CD = [nsig_nogo_CD, goodData(nsig_nogo_ps(i)).chanDepth];
%end

% Get mean firing rates
gos0_fr = mean(all_bins(:, taskLabel == 1 & freqLabel == 1), 2);
gos2_fr = mean(all_bins(:, taskLabel == 3 & freqLabel == 1), 2);
nogos0_fr = mean(all_bins(:, taskLabel == 1 & freqLabel == 2), 2);
nogos2_fr = mean(all_bins(:, taskLabel == 3 & freqLabel == 2), 2);

%--------------------------------------------------------------
% Raster example
% ex_neu = 69 for upregulation gos
%--------------------------------------------------------------
ex_neu = 69;
ex_spks = double(goodData(ex_neu).spikeTimes);

gos0 = struct(); gos2 = struct();
counter0 = 1; counter2 = 1;
for i = 1:length(astamps)
  curr_onset = astamps(i);
  curr_win = [curr_onset+window_st, curr_onset+window_end];
  if taskLabel(i) == 1 & freqLabel(i) == 1 % gos0
    gos0(counter0).times = ex_spks(find(ex_spks > curr_win(1) & ex_spks < curr_win(2))) - curr_onset;
    counter0 = counter0 + 1;
  elseif taskLabel(i) == 3 & freqLabel(i) == 1 % gos2
    gos2(counter2).times = ex_spks(find(ex_spks > curr_win(1) & ex_spks < curr_win(2))) - curr_onset;
    counter2 = counter2 + 1;
  end
end

figure; hold on;
for i = 1:length(gos2)
  plot(gos2(i).times, ones(1, length(gos2(i).times))*i, '.', 'markers', 8, 'color', hex(3, :));
end
for i = 1:length(gos0)
  plot(gos0(i).times, ones(1, length(gos0(i).times))*(i+length(gos2)), '.', 'markers', 8, 'color', hex(3, :));
end
axis tight;
xlim([window_st, window_end]);
plot(xlim, [length(gos2), length(gos2)], 'k--')
plot([0, 0], ylim, 'k--')
%print(fullfile(fig_dir, 'data_gos_up.eps'), '-painters', '-depsc');


ex_neu = 115;
ex_spks = double(goodData(ex_neu).spikeTimes);

nogos0 = struct(); nogos2 = struct();
counter0 = 1; counter2 = 1;
for i = 1:length(astamps)
  curr_onset = astamps(i);
  curr_win = [curr_onset+window_st, curr_onset+window_end];
  if taskLabel(i) == 1 & freqLabel(i) == 2 % nogos0
    nogos0(counter0).times = ex_spks(find(ex_spks > curr_win(1) & ex_spks < curr_win(2))) - curr_onset;
    counter0 = counter0 + 1;
  elseif taskLabel(i) == 3 & freqLabel(i) == 2 % nogos2
    nogos2(counter2).times = ex_spks(find(ex_spks > curr_win(1) & ex_spks < curr_win(2))) - curr_onset;
    counter2 = counter2 + 1;
  end
end

figure; hold on;
for i = 1:length(nogos2)
  plot(nogos2(i).times, ones(1, length(nogos2(i).times))*i, '.', 'markers', 8, 'color', hex(2, :));
end
for i = 1:length(nogos0)
  plot(nogos0(i).times, ones(1, length(nogos0(i).times))*(i+length(nogos2)), '.', 'markers', 8, 'color', hex(2, :));
end
axis tight;
xlim([window_st, window_end]);
plot(xlim, [length(nogos2), length(nogos2)], 'k--')
plot([0, 0], ylim, 'k--')
%print(fullfile(fig_dir, 'data_nogos_down.eps'), '-painters', '-depsc');





