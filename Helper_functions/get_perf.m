function [ tpr, fpr, mcc, auc, ...
    tpr_diff, fpr_diff, mcc_diff, auc_diff, ...
    fl ] = ...
    get_perf( truth, ppi, ppi_diff, true_C, est_C )
% 1. Given truth array p x p x K of true adjacency matrices and results array
% of p x p x K edge selections, calculate true positive rate, false
% positive rate, MCC, and AUC
% 2. Given p * K by p * K matrix of posterior difference probabilities,
% calculate true positive rate, false positive rate, MCC and AUC for
% differential edge selection
% 3. Given true precision matrices + estimated precision matrices,
% compute Frobenius loss

[~, ~, K] = size(truth);

% Treat median model as selections
median_model = ppi > .5;

% 1. ---------------------------------------------------------------------

[tpr, fpr, mcc] = get_tpr_fpr_mcc(truth, median_model);

opts = (0:1000) / 1000;
tpr_roc = zeros(1001, 1);
fpr_roc = zeros(1001, 1);

for i = 1:1001
    cur_threshold = opts(i);
    
    % Get estimated adjacency matrix based on ppi matrices and current threshold
    cur_adj = ppi > cur_threshold;
    
    [tpr_roc(i), fpr_roc(i), ~] = get_tpr_fpr_mcc(truth, cur_adj);
end

auc = sum((fpr_roc(1:1000) - fpr_roc(2:1001)) .* ...
    (tpr_roc(2:1001) + tpr_roc(1:1000)) / 2);

% 2. --------------------------------------------------------------------

% Median model for differential edge selection
median_model_diff = ppi_diff > .5;

[tpr_diff, fpr_diff, mcc_diff] = get_tpr_fpr_mcc_diff(truth, median_model_diff);

tpr_diff_roc = zeros(1001, 1);
fpr_diff_roc = zeros(1001, 1);

for i = 1:1001
    cur_threshold = opts(i);
    
    [tpr_diff_roc(i), fpr_diff_roc(i), ~] = get_tpr_fpr_mcc_diff(truth, ...
        ppi_diff > cur_threshold);
end

auc_diff = sum((fpr_diff_roc(1:1000) - fpr_diff_roc(2:1001)) .* ...
    (tpr_diff_roc(2:1001) + tpr_diff_roc(1:1000)) / 2);

% 3. --------------------------------------------------------------------
fl = 0;

for k = 1:K
    fl = fl + ...
        norm(squeeze(true_C(:, :, k)) - squeeze(est_C(:, :, k)), 'fro')^2 / ...
        norm(squeeze(true_C(:, :, k)), 'fro')^2 / K;
end

end

