function edges = delaunayEdges(tri)

% Function to extract list of edges from old-style MATLAB Delaunay triangulation

% Build list of every side of every triangle
edges = [tri(:,1:2) ; tri(:,2:3) ; tri(:,[3 1])];

% Sort
edges = sort(edges, 2);
edges = sortrows(edges);

% Remove duplicates
dups = ~any(diff(edges, 1), 2);
edges(dups,:) = [];