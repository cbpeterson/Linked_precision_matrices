function [ tpr, fpr, mcc ] = get_tpr_fpr_mcc( truth, results )
% Given truth array p x p x K of true adjacency matrices and results array
% of p x p x K edge selections, calculate true positive and false positive
% rates across the graphs

[~, p, K] = size(truth);

% Get indices of upper triangular entries i.e. unique edges
indmx = reshape([1:p^2], p, p);
upperind = indmx(triu(indmx, 1) > 0);

% Get number of fp, tp, fn and tn
fp = 0;
tp = 0;
fn = 0;
tn = 0;
for k = 1:K
    cur_truth = squeeze(truth(:, :, k));
    true_up = cur_truth(upperind);
    cur_res = squeeze(results(:, :, k));
    res_up = cur_res(upperind);
    tp = tp + sum(true_up & res_up);
    fp = fp + sum(~true_up & res_up);
    tn = tn + sum(~true_up & ~res_up);
    fn = fn + sum(true_up & ~res_up);
end

tpr = tp / (tp + fn);
fpr = fp / (fp + tn);
mcc = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));

end
