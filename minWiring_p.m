function [mw, p] = minWiring_p(x, y, z, nMC, pairs, circular)

% pairs is not used, but is explicitly in the args for compatibility with R2008b

x = x(:);
y = y(:);
z = z(:);

n = length(x);

% Make path circular if necessary
if circular
    fCirc = @(in) [in ; in(1)];
else
    fCirc = @(in) in;
end 

mwSamps = zeros(nMC,1);
mcSamps = zeros(nMC,1);

for i = 1 : nMC
    % Add a small amount of noise to resolve feature space sort order
    [dummy, order] = sort(z + 0.001 .* randn(size(z)));
    order = fCirc(order);
    
    % Compute map space distances
    dx = diff(x(order));
    dy = diff(y(order));
    
    % Compute actual measure sample
    mwSamps(i) = mean(dx.^2 + dy.^2);
    
    % Do MC permutation
    order = randperm(n);
    order = fCirc(order);
    
    dx = diff(x(order));
    dy = diff(y(order));
    mcSamps(i) = mean(dx.^2 + dy.^2);
end

mw = mean(mwSamps);
p = sum(mcSamps < mw) ./ nMC;

end