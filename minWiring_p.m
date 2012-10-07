function [mw, p] = minWiring_p(x, y, z, nMC, pairs, circular)

% pairs is not used, but is explicitly in the args for compatibility with R2008b

nResolve = 10^4;

x = x(:);
y = y(:);
z = z(:);
n = length(x);
zJitterSD = 0.001 .* std(z);

% Make path circular if necessary
if circular
    fCirc = @(in) [in ; in(1)];
else
    fCirc = @(in) in;
end 

% MC to resolve feature space sort order ambiguity
mwSamps = zeros(nResolve,1);

for i = 1 : nResolve
    % Add a small amount of noise to resolve feature space sort order
    [dummy, order] = sort(z + zJitterSD .* randn(size(z)));
    order = fCirc(order);
    
    % Compute map space distances
    dx = diff(x(order));
    dy = diff(y(order));
    
    % Compute measure sample
    mwSamps(i) = mean(dx.^2 + dy.^2);
end

mw = mean(mwSamps);

% Is n small enough to do exact permutation?
nExact = factorial(n);

if nExact > nMC
    % Do MC permutation
    mcSamps = zeros(nMC,1);
    
    for i = 1 : nMC
        order = randperm(n)';
        order = fCirc(order);
        dx = diff(x(order));
        dy = diff(y(order));
        mcSamps(i) = mean(dx.^2 + dy.^2);
    end
    
    p = (sum(mcSamps <= mw) + 1) ./ (nMC + 1);    
else
    % Do exact permutation
    mcSamps = zeros(nExact,1);
    orders = perms(1:n);
    
    for i = 1 : nExact
        order = fCirc(orders(i,:)');
        dx = diff(x(order));
        dy = diff(y(order));
        mcSamps(i) = mean(dx.^2 + dy.^2);
    end
    
    p = sum(mcSamps <= mw) ./ nExact;
end

end