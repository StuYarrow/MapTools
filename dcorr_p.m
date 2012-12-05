function [dcorr, p] = dcorr_p(x, y, z, nMC, pairs, circular)

x = x(:);
y = y(:);
z = z(:);
n = length(x);

% Assign feature space function
if circular
    fdz = @(in) abs(bsxfun(@circ_dist, in, in'));
else
    fdz = @(in) abs(bsxfun(@minus, in, in'));
end

% Setup mask to mark valid pairs
if ~isempty(pairs)
    pairs2 = [pairs ; fliplr(pairs)];
    mask = false([n n]);
    mask(sub2ind([n n], pairs2(:,1), pairs2(:,2))) = true;
else
    mask = true([n n]);
end

% Counts of valid pairs
colCnt = sum(mask, 1);
rowCnt = sum(mask, 2);
gndCnt = sum(mask(:));

% Compute recentered map space distances
dx = bsxfun(@minus, x, x');
dy = bsxfun(@minus, y, y');
dCort = sqrt(dx.^2 + dy.^2);
dCort(~mask) = 0;
dCortR = bsxfun(@minus, dCort, sum(dCort, 1) ./ colCnt);
dCortR = bsxfun(@minus, dCortR, sum(dCort, 2) ./ rowCnt);
dCortR = dCortR + sum(dCort(:)) ./ gndCnt;

% Compute recentered feature space diatances
dz = fdz(z);
dz(~mask) = 0;
dzR = bsxfun(@minus, dz, sum(dz, 1) ./ colCnt);
dzR = bsxfun(@minus, dzR, sum(dz, 2) ./ rowCnt);
dzR = dzR + sum(dz(:)) ./ gndCnt;

% Compute dcorr
dcov = sqrt(mean(dCortR(:) .* dzR(:)));
dvarCort = sqrt(mean(dCortR(:) .* dCortR(:)));
dvarZ = sqrt(mean(dzR(:) .* dzR(:)));
dcorr = dcov / sqrt(dvarCort * dvarZ);

% Is n small enough to do exact permutation?
nExact = factorial(n);
exact = nExact < nMC;

if exact
    nPerm = nExact;
    zPerms = perms(z);
else
    nPerm = nMC;
end

% Do permutation analysis
rArr = zeros(nPerm,1);

for i = 1 : nPerm
    if exact
        % Get next permutation
        zShuffled = zPerms(i,:);
    else
        % Shuffle z values randomly
        zShuffled = z(randperm(length(z)));
    end
    
    % Compute new feature space distances
    dzShuffled = fdz(zShuffled);
    dzShuffled(~mask) = 0;
    dzRshuf = bsxfun(@minus, dzShuffled, sum(dzShuffled, 1) ./ colCnt);
    dzRshuf = bsxfun(@minus, dzRshuf, sum(dzShuffled, 2) ./ rowCnt);
    dzRshuf = dzRshuf + sum(dzShuffled(:)) ./ gndCnt;
    
    % Compute dcorr for sample
    dcov = sqrt(mean(dCortR(:) .* dzRshuf(:)));
    dvarZ = sqrt(mean(dzRshuf(:) .* dzRshuf(:)));
    rArr(i) = dcov / sqrt(dvarCort * dvarZ);
end

% Compute p-value
if exact
    p = sum(rArr >= dcorr) ./ nExact;
else
    p = (sum(rArr >= dcorr) + 1) ./ (nMC + 1);
end
    
end

