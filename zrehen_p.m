function [zrehen, p] = zrehen_p(x, y, z, nMC, pairs, circular)

% pairs is not used, but is explicitly in the args for compatibility with R2008b

x = x(:);
y = y(:);
z = z(:);

n = length(x);

% Define intruder counting functions
    function d = fInt(r)
        d = abs(r(edges(:,1)) - r(edges(:,2))) - 1;
    end

    function d = fInt_circ(r)
        d = abs(r(edges(:,1)) - r(edges(:,2)));
        d(d > n) = n - d(d > n); % Resolve 'long way round' distances
        d = d - 1;
    end

% Find neighbours in map space
dt = delaunay(x, y);
edges = delaunayEdges(dt);

% Count intruders
if circular
    % Make sure we don't have -pi AND pi values as this will mess up the ranks
    z(z == pi) = -pi;
    zr = ranks(z);
    
    fIntruders = @fInt_circ;
else
    zr = ranks(z);
    
    fIntruders = @fInt;
end

% Compute measure
zrehen = mean(fIntruders(zr)) ./ n;

% Do MC permutation
samps = zeros(nMC,1);
for i = 1 : nMC
    zrShuf = zr(randperm(n));
    samps(i) = mean(fIntruders(zrShuf)) ./ n;
end

p = sum(samps < zrehen) ./ nMC;

end