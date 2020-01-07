function [ tpr_diff, fpr_diff, mcc_diff ] = get_tpr_fpr_mcc_diff( truth, results_diff )
% Given truth array p x p x K of true adjacency matrices and results array
% of p x K by p x K differential edge selections, compute tpr, fpr and mcc
% of differential edge selection

[~, p, K] = size(truth);

% Get number of fp, tp, fn and tn for differential edge selection
fp = 0;
tp = 0;
fn = 0;
tn = 0;
for k1 = 1:(K - 1)
    for k2 = (k1 + 1):K
        cur_truth1 = truth(:, :, k1);
        cur_truth2 = truth(:, :, k2);
        
        start_row = (k1 - 1) * p + 1;
        end_row = k1 * p;
        start_col = (k2 - 1) * p + 1;
        end_col = k2 * p;
        cur_res_diff = results_diff(start_row:end_row, start_col:end_col);
        
        for i = 1:(p - 1)
            for j = (i + 1):p
                true_diff = abs(cur_truth1(i, j) - cur_truth2(i, j));
                res_diff = cur_res_diff(i, j);
                tp = tp + sum(true_diff & res_diff);
                fp = fp + sum(~true_diff & res_diff);
                tn = tn + sum(~true_diff & ~res_diff);
                fn = fn + sum(true_diff & ~res_diff);
            end
        end
    end
end

tpr_diff = tp / (tp + fn);
fpr_diff = fp / (fp + tn);
mcc_diff = (tp * tn - fp * fn) / sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn));

end
