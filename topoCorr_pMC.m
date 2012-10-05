function [tc, p] = topoCorr_pMC(x, y, z, nMC, pairs, circular)

% pairs is not used, but is explicitly in the args for compatibility with R2008b

x = x(:);
y = y(:);
z = z(:);

n = length(x);
mask = ~~tril(ones(n), -1);

% Define feature space distance functions
    function d = fdz(r)
        d = abs(bsxfun(@minus, r, r'));
        d = d(mask);
    end

    function d = fdz_circ(r)
        d = abs(bsxfun(@minus, r, r'));
        d(d > n) = n - d(d > n); % Resolve 'long way round' distances
        d = d(mask);
    end

% Assign distance function
if circular
    fDist = @fdz_circ;
else
    fDist = @fdz;
end

% map space triangulation
dt = delaunay(x, y);
edges = delaunayEdges(dt);
mapGraph = sparse(edges(:,2), edges(:,1), ones(size(edges(:,2))), n, n, length(edges));

mapDists = zeros(n);
for i = 1 : n
    mapDists(i,:) = graphshortestpath(mapGraph, i, 'directed', false, 'method', 'BFS');
end

mapDists = mapDists(mask);
dMapDists = mapDists - mean(mapDists);

featRank = ranks(z);
featDist = fDist(featRank);
dFeatDists = featDist - mean(featDist);

tc = sum(dFeatDists .* dMapDists) ./ sqrt(sum(dFeatDists .^ 2) .* sum(dMapDists .^ 2));

% Is n small enough to do exact permutation?
if factorial(n) > nMC
    % Do MC permutation
    Tp = zeros(nMC, 1);
    for i = 1 : nMC
        featRankShuf = featRank(randperm(n));
        featDist = fDist(featRankShuf);
        dFeatDists = featDist - mean(featDist);
        
        Tp(i) = sum(dFeatDists .* dMapDists) ./ sqrt(sum(dFeatDists .^ 2) .* sum(dMapDists .^ 2));
    end
    
    p = (sum(Tp > tc) + 1) ./ (nMC + 1);
else
    % Do exact permutation
    nMC = factorial(n);
    Tp = zeros(nMC, 1);
    featRankShuf = perms(featRank);
    
    for i = 1 : nMC
        featDist = fDist(featRankShuf(i,:));
        dFeatDists = featDist - mean(featDist);
        
        Tp(i) = sum(dFeatDists .* dMapDists) ./ sqrt(sum(dFeatDists .^ 2) .* sum(dMapDists .^ 2));
    end
    
    p = sum(Tp > tc) ./ nMC;
end

end