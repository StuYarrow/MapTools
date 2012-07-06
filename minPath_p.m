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
if circular
    fPath = @(in) mean( (in(edges(:,1)) - in(edges(:,2))).^2 );
else
    fPath = @(in) mean( circ_dist(in(edges(:,1)), in(edges(:,2))).^2 );
end

% Compute measure
mp = fPath(z);

% Do MC permutation
mcSamps = zeros(nMC,1);
for i = 1 : nMC
    zShuf = z(randperm(n));
    mcSamps(i) = fPath(zShuf);
end

p = sum(mcSamps < mp) ./ nMC;

end
