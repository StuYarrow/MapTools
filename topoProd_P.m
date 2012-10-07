function [tp, p] = topoProd_P(x, y, z, nMC, pairs, circular)

% pairs is not used, but is explicitly in the args for compatibility with R2008b

nMC2 = 100;

if length(y) ~= length(x) || length(z) ~= length(x)
    error('All vectors must be same length')
end

x = x(:);
y = y(:);
z = z(:);

n = length(x);
In = ~~eye(n);
sortIndices = repmat((1:n)', [1 n]);

% Define feature space distance functions
if circular
    fdz = @(in) abs(circ_dist2(in));
else
    fdz = @(in) abs(bsxfun(@minus, in, in'));
end

mapDist = sqrt(bsxfun(@minus, x, x').^2 + bsxfun(@minus, y, y').^2);
mapDist(In) = -1; % make sure self distances sort first
[mapDistSort mapInd] = sort(mapDist, 2);

% Allocate memory for samples
tpSamples = zeros(nMC,1);
tpNoise = zeros(nMC2,1);

k = repmat(1 : n-1, [n 1]);
kk = 1 ./ (2 .* k);

% Do initial short-run MC to resolve identical values in feature space
for i = 1 : nMC2
    zNoise = z + 0.001 .* randn(size(z));
    featDist = fdz(zNoise);
    featDist(In) = -1;
    [featDistSort featInd] = sort(featDist, 2);
    
    mapDistSortFeat = mapDist(sub2ind([n n], sortIndices, featInd));
    featDistSortMap = featDist(sub2ind([n n], sortIndices, mapInd));
    
    Q1 = featDistSortMap ./ featDistSort;
    Q2 = mapDistSort ./ mapDistSortFeat;
    
    % Discard zeroth-order neighbours (self)
    Q1(:,1) = [];
    Q2(:,1) = [];
    
    logP3 = cumsum(log(Q1) + log(Q2), 2) .* kk;
    tpNoise(i) = 1 / (n^2 - n) .* sum(logP3(:));
end

% Do MC permutation
for i = 1 : nMC
    zNoise = z + 0.001 .* randn(size(z));
    shuffled = zNoise(randperm(n));
    featDist = fdz(shuffled);
    featDist(In) = -1;
    [featDistSort featInd] = sort(featDist, 2);
    
    mapDistSortFeat = mapDist(sub2ind([n n], sortIndices, featInd));
    featDistSortMap = featDist(sub2ind([n n], sortIndices, mapInd));
    
    Q1 = featDistSortMap ./ featDistSort;
    Q2 = mapDistSort ./ mapDistSortFeat;

    % Discard zeroth-order neighbours (self)
    Q1(:,1) = [];
    Q2(:,1) = [];
    
    logP3 = cumsum(log(Q1) + log(Q2), 2) .* kk;    
    tpSamples(i) = 1 / (n^2 - n) .* sum(logP3(:));
end

tp = mean(tpNoise);
p = (sum(abs(tpSamples) <= abs(tp)) + 1) ./ (nMC + 1);

end