function [r, p] = pearson_pairs_p(varargin)

[r, p] = corr_pairs_p('pearson', varargin{:});

end