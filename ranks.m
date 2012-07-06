function r = ranks(in)

[sorted, order] = sort(in);

id = diff(sorted) == 0;
id = id(:)';
blockStart = find([id false] & ~[false id]);
blockEnd = find(~[id false] & [false id]);
nBlocks = length(blockStart);

r = reshape(1 : length(in), size(in));

for i = 1 : nBlocks
    s = blockStart(i);
    e = blockEnd(i);
    r(s:e) = sum(r(s:e)) ./ (e - s + 1);
end

r(order) = r;
end