% Name: Robert Kim
% Date: 03-29-2023
% Email: robert.f.kim@gmail.com 
% fnc_eval_model.m
% Description:

function out = fnc_eval_model(model_dir, uu, activation_fnc, mod_ww, lesion_type)
% INPUT
%   model_dir: full path (including .mat filename) for a trained RNN model
%   uu: input matrix (channel x time)
%   activation_fnc: activation function (sigmoid or clipped ReLU)
%   mod_ww: modified recurrent connectivity matrix. No modification if left blank ([])
%   lesion_type: 'ii', 'ee', 'ei', 'ie'
% OUTPUT
%   out: container that contains
%     out('X'): synpatic current matrix (time x neurons)
%     out('R'): firing rate estimate matrix (time x neurons)
%     out('O'): RNN output signal (time x 1)

% Load the data first
load(model_dir);

% Time points-by-neurons
X = zeros(size(x)); % x is the synaptic current variable
R = zeros(size(r)); % r is the firing rate estimate

X(1, :) = x(1, :);
if strcmpi(activation_fnc, 'sigmoid') == 1
  R(1, :) = 1./(1+exp(-X(1, :)));
elseif strcmpi(activation_fnc, 'clipped_relu') == 1
  R(1, :) = max(min(X(1, :), 1), 0);
end

% Get the trained recurrent connectivity weight
% w is the non-negative weight matrix
% m is the mask matrix (diagonal matrix) that indicates which
% neurons are inhibitory/excitatory units.
% To get the weight matrix that abides by Dale's principle,
% we need to multiple w by m.
intact_ww = w*m.*som_m;
if isempty(mod_ww)
  ww = intact_ww;
else
  ww = mod_ww;
end

inh_ind = find(inh);
exc_ind = find(exc);
if strcmpi(lesion_type, 'ii')
  ww(inh_ind, inh_ind) = ww(inh_ind, inh_ind)*0;
elseif strcmpi(lesion_type, 'ie')
  ww(inh_ind, exc_ind) = ww(inh_ind, exc_ind)*0;
elseif strcmpi(lesion_type, 'ee')
  ww(exc_ind, exc_ind) = ww(exc_ind, exc_ind)*0;
elseif strcmpi(lesion_type, 'ei')
  ww(exc_ind, inh_ind) = ww(exc_ind, inh_ind)*0;
end

% Synaptic decay constants constrained by min and max limit
% taus = [min_tau, max_tau]
taus_sig = (1./(1+exp(-taus_gaus)))*(taus(2) - taus(1)) + taus(1);

% Total number of time-points
T = size(uu, 2);
DeltaT = 1;

O = zeros(1, T);
for t = 2:T+1

  next_x = (1-DeltaT./taus_sig).*transpose(X(t-1, :)) + ...
      (DeltaT./taus_sig).*(ww*transpose(R(t-1, :)) + ...
      w_in*uu(:, t-1)) + randn(N, 1)/10;

  if strcmpi(activation_fnc, 'sigmoid') == 1
    next_r = 1./(1+exp(-next_x));
  elseif strcmpi(activation_fnc, 'clipped_relu') == 1
    next_r = max(min(next_x, 1), 0);
  end
  
  X(t, :) = transpose(next_x);
  R(t, :) = transpose(next_r);
  O(t) = w_out*transpose(R(t, :)) + b_out;
end

out = containers.Map;

out('X') = X(2:end, :);
out('R') = R(2:end, :);
out('O') = O(2:end);


