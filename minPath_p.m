function [mp, p] = minPath_p(x, y, z, nMC, pairs, circular)

% pairs is not used, but is explicitly in the args for compatibility with R2008b

x = x(:);
y = y(:);
z = z(:);

n = length(x);

% Compute map space triangulation
dt = delaunay(x, y);
edges = delaunayEdges(dt);

% Define path length functions
if ~circular
    fPath = @(in) mean( (in(edges(:,1)) - in(edges(:,2))).^2 );
else
    fPath = @(in) mean( circ_dist(in(edges(:,1)), in(edges(:,2))).^2 );
end

% Compute measure
mp = fPath(z);

% Is n small enough to do exact permutation?
nExact = factorial(n);

if nExact > nMC
    % Do MC permutation
    mcSamps = zeros(nMC,1);
    
    for i = 1 : nMC
        zShuf = z(randperm(n));
        mcSamps(i) = fPath(zShuf);
    end

    p = (sum(mcSamps <= mp) + 1) ./ (nMC + 1);
else
    % Do exact permutation
    mcSamps = zeros(nExact,1);
    zShuf = perms(z);
    
    for i = 1 : nExact
        mcSamps(i) = fPath(zShuf(i,:));
    end
    
    p = sum(mcSamps <= mp) ./ nExact;
end

end
