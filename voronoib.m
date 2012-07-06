function [verts, cells, allpoints] = voronoib(points)

n = size(points, 1);
[v, c] = voronoin(points);

unbounded = false(n,1);
areas = zeros(n,1);

for i = 1 : n
    % test for unbounded cells
    if any(c{i} == 1)
        unbounded(i) = true;
    else
        areas(i,1) = polyarea(v(c{i},1), v(c{i},2));
    end
end

meanArea = median(areas(~unbounded));
radius = sqrt(meanArea ./ pi);

hull = convhull(points(:,1), points(:,2));
hullp = [hull(end-1) ; hull];
priPoints = [];
secPoints = [];

for i = 2 : size(hullp,1) - 1
    vec1 = points(hullp(i-1),:) - points(hullp(i),:);
    vec2 = points(hullp(i+1),:) - points(hullp(i),:);
    vec1 = vec1 ./ norm(vec1);
    vec2 = vec2 ./ norm(vec2);
    
    vecIn = (vec1 + vec2) ./ norm(vec1 + vec2);
    theta = acos(dot(vec1, vec2) ./ (norm(vec1) .* norm(vec2)));
    offset = 2 * radius ./ sin(theta ./ 2);
    
    if ~any(isnan(vecIn))
        priPoints(end+1,:) = points(hullp(i),:) - (offset .* vecIn);
    end
end

for i = 1 : size(priPoints,1)
    if i == 1
        vec = priPoints(end,:) - priPoints(i,:);
    else
        vec = priPoints(i-1,:) - priPoints(i,:);
    end
    
    nSeg = floor(norm(vec) ./ (2 * radius));
    
    for j = 1 : nSeg - 1
        secPoints(end+1,:) = priPoints(i,:) + j .* vec ./ nSeg;
    end
end

allpoints = [points ; priPoints ; secPoints];
[v2, c2] = voronoin(allpoints);
verts = v2;
cells = c2(1:n,:);

%figure

%for i = 1 : length(cells)
%    patch(verts(cells{i},1), verts(cells{i},2), [0.9 0.9 0.9])
%end

%hold on
%box on
%plot(points(hull,1), points(hull,2), 'k--')
%plot(points(:,1), points(:,2), 'k+')
%plot(priPoints(:,1), priPoints(:,2), 'rx')
%plot(secPoints(:,1), secPoints(:,2), 'rx')

end