function r = ranks(in)

% Sort values
[sorted, order] = sort(in);

% Locate values that are identical to their neighbours
id = diff(sorted) == 0;
id = id(:)';

% Locate the beginnings and ends of blocks of identical values
blockStart = find([id false] & ~[false id]);
blockEnd = find(~[id false] & [false id]);
nBlocks = length(blockStart);

% Generate uncorrected ranks
r = reshape(1 : length(in), size(in));

% Correct for rank ties by averaging rank within identical blocks
for i = 1 : nBlocks
    s = blockStart(i);
    e = blockEnd(i);
    r(s:e) = sum(r(s:e)) ./ (e - s + 1);
end

% Rearrange into correct order
r(order) = r;
end