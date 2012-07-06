function [r, p] = corr_pairs_p(meth, x, y, z, n, pairs, circular)

x = x(:);
y = y(:);
z = z(:);

% Define correlation functions
    function rho = pearson(x1, x2)
        rho = corrcoef(x1, x2);
        rho = rho(1,2);
    end

    function rho = spearman(x1, x2)
        r1 = ranks(x1);
        r2 = ranks(x2);
        rho = corrcoef(r1, r2);
        rho = rho(1,2);
    end

% Define feature space distance functions
    function d = fdz_pr(in)
        d = abs(in(pairs(:,1)) - in(pairs(:,2)));
    end

    function d = fdz_pr_circ(in)
        d = abs(circ_dist(in(pairs(:,1)), in(pairs(:,2))));
    end

    function d = fdz_npr(in)
        d = abs(bsxfun(@minus, in, in'));
        d = d(mask);
    end

    function d = fdz_npr_circ(in)
        d = abs(bsxfun(@circ_dist, in, in'));
        d = d(mask);
    end

% Assign correlation function
switch meth
    case 'pearson'
        fcorr = @pearson;
    case 'spearman'
        fcorr = @spearman;
    otherwise
        error('Unknown correlation method: %s', meth)
end

% Compute map space distances and assign feature space function
if ~isempty(pairs)
    dx = x(pairs(:,1)) - x(pairs(:,2));
    dy = y(pairs(:,1)) - y(pairs(:,2));
    dCort = sqrt(dx.^2 + dy.^2);
    
    if circular
        fdz = @fdz_pr_circ;
    else
        fdz = @fdz_pr;
    end
else
    mask = triu(true(length(x)), 1);
    dx = bsxfun(@minus, x, x');
    dy = bsxfun(@minus, y, y');
    dCort = sqrt(dx(mask).^2 + dy(mask).^2);
    
    if circular
        fdz = @fdz_npr_circ;
    else
        fdz = @fdz_npr;
    end
end

% Compute feature space diatances
dz = fdz(z);

% Compute correlation coefficient
r = fcorr(dCort, dz);

% Do MC permutation
rArr = zeros(n,1);

for i = 1 : n
    % Shuffle z values
    zShuffled = z(randperm(length(z)));
    
    % Compute new feature space distances
    dzShuffled = fdz(zShuffled);
    
    % Compute rho for sample
    rArr(i) = fcorr(dCort, dzShuffled);
end

% Compute p value
p = sum(rArr > r) ./ n;

end