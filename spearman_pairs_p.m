function [r, p] = spearman_pairs_p(varargin)

[r, p] = corr_pairs_p('spearman', varargin{:});

end