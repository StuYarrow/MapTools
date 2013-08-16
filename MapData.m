classdef MapData
    %MapData Dataset representing a neural map
    %   map = MapData(mapSpaceCoords, featureSpaceCoords [,mapSpaceAngFlags [,featSpaceAngFlags [,pairs]]])
    %   
    
    properties (SetAccess=protected)
        map = [];
        ftr = [];
        mapAngular = [false false];
        ftrAngular = [];
        pairs = false;
        n = 0;
        mapDims = 2;
        ftrDims = 1;
    end
    
    methods
        
        function obj = MapData(varargin)
            % MapData( mapSpaceCoords, featureSpaceCoords [, mapSpaceAngFlags [, featSpaceAngFlags [, pairs] ] ] )
            assert(nargin ~= 1 && nargin <= 5, 'wrong number of arguments')
            
            if nargin >= 2
                mapIn = varargin{1};
                ftrIn = varargin{2};
                assert(size(mapIn, 1) == size(ftrIn, 1), 'mismatch in number of map space and feature space points')

                obj.mapDims = size(mapIn, 2);
                obj.ftrDims = size(ftrIn, 2);
                obj.n = size(mapIn, 1);
                obj.map = mapIn;
                obj.ftr = ftrIn;
            end
            
            if nargin >= 3
                circIn = varargin{3};
                
                assert(numel(circIn) == obj.mapDims, 'mismatch between map space dimensionality and length of circular flag vector')
                obj.mapAngular = circIn(:)';
                
                if any(obj.mapAngular)
                    minAng = min(min(obj.map(:,obj.mapAngular)));
                    maxAng = max(max(obj.map(:,obj.mapAngular)));
                    assert(minAng > -pi && maxAng <= pi, 'map space angular coordinates must be in the interval (-pi pi]')
                end
            end
            
            if nargin >= 4
                circIn = varargin{4};
                
                assert(numel(circIn) == obj.ftrDims, 'mismatch between feature space dimensionality and length of circular flag vector')
                obj.ftrAngular = circIn(:)';
                
                if any(obj.ftrAngular)
                    minAng = min(min(obj.ftr(:,obj.ftrAngular)));
                    maxAng = max(max(obj.ftr(:,obj.ftrAngular)));
                    assert(minAng > -pi && maxAng <= pi, 'feature space angular coordinates must be in the interval (-pi pi]')
                end
            end
            
            if nargin >= 5
                pairsIn = varargin{5};
                assert(size(pairsIn,2) == 2, 'pairs matrix must be m by 2')
                assert(sum(pairsIn(:) - round(pairsIn(:))) == 0, 'pairs indices must be integers')
                assert(max(pairsIn(:)) <= obj.n, 'all elements of pairs must be valid indices')
                obj.pairs = pairsIn;
            end
        end
        
        function [c, p] = measure(obj, m, nMC)
            fNames = obj.measures('funcs');
            [c, p] = obj.(fNames{m})(nMC);
        end
        
        function [c, p] = pearson_pairs_p(obj, nMC)
            [c, p] = obj.corr_pairs_p('pearson', nMC);
        end
        
        function [c, p] = spearman_pairs_p(obj, nMC)
            [c, p] = obj.corr_pairs_p('spearman', nMC);
        end
        
        function [c, p] = zrehen_p(obj, nMC)
            assert(~obj.pairs, 'does not support specified pairs/merged maps')
            
            % Find neighbours in map space
            assert(obj.mapDims > 1, 'zrehen measure does not support 1-D map spaces')
            edges = obj.neighbours(obj.map, obj.mapAngular);
            
            % Set up for intruder counting
            % Compute radii for intruder test
            ftrRadii = 0.5 .* obj.distMatrix(obj.ftr, obj.ftrAngular);
            % Compute midpoints
            ftrRe1 = repmat(permute(obj.ftr, [1 3 2]), [1 obj.n 1]);
            ftrRe2 = repmat(permute(obj.ftr, [3 1 2]), [obj.n 1 1]);
            % linear dimension midpoints
            ftrMidpoints(:,:,~obj.ftrAngular) = mean(cat(4, ftrRe1(:,:,~obj.ftrAngular), ftrRe2(:,:,~obj.ftrAngular)), 4);
            % circular dimension midpoints
            ftrMidpoints(:,:,obj.ftrAngular) = circ_mean(cat(4, ftrRe1(:,:,obj.ftrAngular), ftrRe2(:,:,obj.ftrAngular)), [], 4);
            
            intruders = zeros(obj.n);
            dMP = zeros(obj.n, obj.n, obj.ftrDims);
            for i = 1 : obj.n
                % distance from midpoint
                dMP(:,:,~obj.ftrAngular) = bsxfun(@minus, ftrMidpoints(:,:,~obj.ftrAngular), ftrRe1(i,1,~obj.ftrAngular));
                dMP(:,:,obj.ftrAngular) = bsxfun(@circ_dist, ftrMidpoints(:,:,obj.ftrAngular), ftrRe1(i,1,obj.ftrAngular));
                dMP = sqrt(sum(dMP.^2, 3));
                % test for intruders
                intrd = dMP < ftrRadii;
                intrd(i,:) = false;
                intrd(:,i) = false;
                intruders = intruders + intrd;
            end
            
            % Compute measure
            c = mean(intruders(sub2ind([obj.n obj.n], edges(:,1), edges(:,2)))) ./ obj.n;
            
            % permutation count for exact test
            nExact = factorial(obj.n);
            
            if nExact > nMC
                % Do MC permutation
                samps = zeros(nMC,1);
                for i = 1 : nMC
                    % shuffle intruder count matrix
                    perm = randperm(obj.n);
                    intruderShuf = intruders(perm,perm);
                    samps(i) = mean(intruderShuf(sub2ind([obj.n obj.n], edges(:,1), edges(:,2)))) ./ obj.n;
                end
            else
                % Do exact permutation
                samps = zeros(nExact,1);
                permInds = perms(1:obj.n);
                for i = 1 : nExact
                    % shuffle intruder count matrix
                    perm = permInds(i,:);
                    intruderShuf = intruders(perm,perm);
                    samps(i) = mean(intruderShuf(sub2ind([obj.n obj.n], edges(:,1), edges(:,2)))) ./ obj.n;
                end
            end
            
            % Compute p value
            if nExact > nMC
                p = (sum(samps <= c) + 1) ./ (nMC + 1);
            else
                p = sum(samps <= c) ./ nExact;
            end
        end
        
        function [c, p] = wiringLength_p(obj, nMC)
            assert(~obj.pairs, 'does not support specified pairs/merged maps')
            
            switch obj.ftrDims
                case 1 % 1-D feature space case
                    nResolve = 10^4;
                    jitterSD = 0.001 .* obj.std(obj.ftr, obj.ftrAngular);

                    % MC to resolve feature space sort order ambiguity
                    resSamps = zeros(nResolve,1);
                    
                    if obj.ftrAngular
                        fOrder2Pairs = @(ord) [ord circshift(ord, 1)];
                    else
                        fOrder2Pairs = @(ord) [ord(1:end-1) ord(2:end)];
                    end
                    
                    for i = 1 : nResolve
                        % Add a small amount of noise to resolve feature space sort order
                        [~, order] = sort(obj.ftr + jitterSD .* randn([obj.n 1]));
                        prs = fOrder2Pairs(order(:));
                        % Compute mean map space distances
                        resSamps(i) = mean(obj.distEuclidPairs2(obj.map, prs, obj.mapAngular));
                    end

                    cNonNorm = mean(resSamps);
                    normTerm = 1 ./ mean(obj.distEuclidPairs2(obj.map, nchoosek(1:obj.n, 2), obj.mapAngular));
                    c = cNonNorm .* normTerm;
            
                    % Is n small enough to do exact permutation?
                    nExact = factorial(obj.n);

                    if nExact > nMC
                        % Do MC permutation
                        permSamps = zeros(nMC,1);

                        for i = 1 : nMC
                            order = randperm(obj.n);
                            prs = fOrder2Pairs(order(:));
                            permSamps(i) = mean(obj.distEuclidPairs2(obj.map, prs, obj.mapAngular));
                        end

                        p = (sum(permSamps <= cNonNorm) + 1) ./ (nMC + 1);
                    else
                        % Do exact permutation
                        permSamps = zeros(nExact,1);
                        orders = perms(1:obj.n);
                        
                        for i = 1 : nExact
                            order = orders(i,:);
                            prs = fOrder2Pairs(order(:));
                            permSamps(i) = mean(obj.distEuclidPairs2(obj.map, prs, obj.mapAngular));
                        end

                        p = sum(permSamps <= c) ./ nExact;
                    end
                    
                case {2, 3} % 2-D or 3-D feature space cases
                    % identify feature space neighbours
                    g = obj.constructGraph(obj.ftr, obj.ftrAngular);
                    [e1, e2] = find(g);
                    edges = [e1 e2];
                    
                    % compute measure
                    cNonNorm = mean(obj.distEuclidPairs2(obj.map, edges, obj.mapAngular));
                    normTerm = 1 ./ mean(obj.distEuclidPairs2(obj.map, nchoosek(1:obj.n, 2), obj.mapAngular));
                    c = cNonNorm .* normTerm;
                    
                    % Permutation analysis.  Slightly different to 1-D case,
                    % but equivalent
                    % Is n small enough to do exact permutation?
                    nExact = factorial(obj.n);

                    if nExact > nMC
                        % Do MC permutation
                        permSamps = zeros(nMC,1);
                        
                        for i = 1 : nMC
                            % Shuffle map space positions
                            mapShuffled = obj.map(randperm(obj.n),:);
                            % Compute new map space distances
                            permSamps(i) = mean(obj.distEuclidPairs2(mapShuffled, edges, obj.mapAngular));
                        end
                        
                        p = (sum(permSamps <= cNonNorm) + 1) ./ (nMC + 1);
                    else
                        % Do exact permutation
                        permSamps = zeros(nExact,1);
                        permInds = perms(1:obj.n);
                        
                        for i = 1 : nExact
                            % Shuffle map space positions
                            mapShuffled = obj.map(permInds(i,:),:);
                            % Compute new map space distances
                            permSamps(i) = mean(obj.distEuclidPairs2(mapShuffled, edges, obj.mapAngular));
                        end

                        p = sum(permSamps <= cNonNorm) ./ nExact;
                    end
                    
                otherwise
                    error('wiring length measure only supports feature spaces with 1-3 dimensions')
            end
        end
        
        function [c, p] = pathLength_p(obj, nMC)
            assert(~obj.pairs, 'does not support specified pairs/merged maps')
            
            % Find neighbours in map space
            assert(obj.mapDims > 1, 'min path measure does not support 1-D map spaces')
            edges = obj.neighbours(obj.map, obj.mapAngular);
            
            % Compute measure
            cNonNorm = mean(obj.distEuclidPairs2(obj.ftr, edges, obj.ftrAngular));
            normTerm = 1 ./ mean(obj.distEuclidPairs2(obj.ftr, nchoosek(1:obj.n, 2), obj.ftrAngular));
            c = cNonNorm .* normTerm;
            
            % Is n small enough to do exact permutation?
            nExact = factorial(obj.n);

            if nExact > nMC
                % Do MC permutation
                samps = zeros(nMC,1);

                for i = 1 : nMC
                    % Shuffle feature space positions
                    ftrShuffled = obj.ftr(randperm(obj.n),:);
                    % Compute sample
                    samps(i) = mean(obj.distEuclidPairs2(ftrShuffled, edges, obj.ftrAngular));
                end

                p = (sum(samps <= cNonNorm) + 1) ./ (nMC + 1);
            else
                % Do exact permutation
                samps = zeros(nExact,1);
                permInds = perms(1:obj.n);

                for i = 1 : nExact
                    % Shuffle feature space positions
                    ftrShuffled = obj.ftr(permInds(i,:),:);
                    % Compute sample
                    samps(i) = mean(obj.distEuclidPairs2(ftrShuffled, edges, obj.ftrAngular));
                end

                p = sum(samps <= cNonNorm) ./ nExact;
            end
        end
           
        function [c, p] = topoProd_p(obj, nMC)
            assert(~obj.pairs, 'does not support specified pairs/merged maps')
            In = ~~eye(obj.n);
            sortIndices = repmat((1:obj.n)', [1 obj.n]);
            
            nResolve = 10^4;
            sampsResolve = zeros(nResolve,1);
            
            k = repmat(1 : obj.n-1, [obj.n 1]);
            kk = 1 ./ (2 .* k);
            
            % Compute distances
            mapDist = obj.distMatrix(obj.map, obj.mapAngular);
            ftrDist = obj.distMatrix(obj.ftr, obj.ftrAngular);
            % make sure self distances sort first
            mapDist(In) = -1;
            ftrDist(In) = -1;

            % Do initial short-run MC to resolve identical values
            % i.e. multiple cells with identical characteristic stimuli or 
            % from same electrode penetration
            for i = 1 : nResolve
                % shuffle both distance matrices identically
                perm = randperm(obj.n);
                mapShuf = mapDist(perm,perm);
                ftrShuf = ftrDist(perm,perm);
                % do sorting
                [mapDistSort, mapInd] = sort(mapShuf, 2);
                [featDistSort, featInd] = sort(ftrShuf, 2);
                % construct cross-sorted matrices
                mapDistSortFeat = mapShuf(sub2ind([obj.n obj.n], sortIndices, featInd));
                featDistSortMap = ftrShuf(sub2ind([obj.n obj.n], sortIndices, mapInd));
                % compute ratios, discarding zeroth-order neighbours (self)
                Q1 = featDistSortMap(:,2:end) ./ featDistSort(:,2:end);
                Q2 = mapDistSort(:,2:end) ./ mapDistSortFeat(:,2:end);
                % compute measure
                logP3 = cumsum(log(Q1) + log(Q2), 2) .* kk;
                sampsResolve(i) = mean(abs(logP3(:)));
            end

            % Allocate memory for samples
            sampsMC = zeros(nMC,1);
            
            % Do MC permutation
            % No exact option in this case as we need to jitter the data to
            % resolve ties in the sort orders
            for i = 1 : nMC
                % shuffle distance matrices independently
                mapPerm = randperm(obj.n);
                ftrPerm = randperm(obj.n);
                mapShuf = mapDist(mapPerm,mapPerm);
                ftrShuf = ftrDist(ftrPerm,ftrPerm);
                % do sorting
                [mapDistSort, mapInd] = sort(mapShuf, 2);
                [featDistSort, featInd] = sort(ftrShuf, 2);
                % construct cross-sorted matrices
                mapDistSortFeat = mapShuf(sub2ind([obj.n obj.n], sortIndices, featInd));
                featDistSortMap = ftrShuf(sub2ind([obj.n obj.n], sortIndices, mapInd));
                % compute ratios, discarding zeroth-order neighbours (self)
                Q1 = featDistSortMap(:,2:end) ./ featDistSort(:,2:end);
                Q2 = mapDistSort(:,2:end) ./ mapDistSortFeat(:,2:end);
                % compute measure
                logP3 = cumsum(log(Q1) + log(Q2), 2) .* kk;
                sampsMC(i) = mean(abs(logP3(:)));
            end

            c = mean(sampsResolve);
            p = (sum(sampsMC <= c) + 1) ./ (nMC + 1);
        end
        
        function [c, p] = topoCorr_p(obj, nMC)
            assert(~obj.pairs, 'does not support specified pairs/merged maps')            
            % Construct mask for pulling out lower triangle
            mask = ~~tril(ones(obj.n), -1);
            
            % compute map space distances
            mapGraph = obj.constructGraph(obj.map, obj.mapAngular);
            mapDistMat = zeros(obj.n);
            for i = 1 : obj.n
                mapDistMat(i,:) = graphshortestpath(mapGraph, i, 'directed', false, 'method', 'BFS');
            end
            
            % compute feature space distances
            if obj.ftrDims == 1 % 1-D case
                rnk = ranks(obj.ftr);
                ftrDistMat = abs(bsxfun(@minus, rnk, rnk'));
                if obj.ftrAngular
                    % Resolve 'long way round' distances
                    lwr = ftrDistMat > obj.n/2;
                    ftrDistMat(lwr) = obj.n - ftrDistMat(lwr);
                end
            else % 2-D or higher
                ftrGraph = obj.constructGraph(obj.ftr, obj.ftrAngular);
                ftrDistMat = zeros(obj.n);
                for i = 1 : obj.n
                    ftrDistMat(i,:) = graphshortestpath(ftrGraph, i, 'directed', false, 'method', 'BFS');
                end
            end
            
            % subtract out means, keep full feature distance matrix for
            % shuffling
            dMapDists = mapDistMat(mask) - mean(mapDistMat(mask));
            dFtrDistMat = ftrDistMat - mean(ftrDistMat(mask));
            % Compute measure
            normTerm = 1 ./ sqrt(sum(dFtrDistMat(mask) .^ 2) .* sum(dMapDists .^ 2));
            c = sum(dFtrDistMat(mask) .* dMapDists) .* normTerm;

            % Is n small enough to do exact permutation?
            nExact = factorial(obj.n);
            exact = nExact <= nMC;

            if exact
                nSamps = nExact;
                permInds = perms(1:obj.n);
            else
                nSamps = nMC;
            end
            
            samps = zeros(nSamps, 1);
            % Do permutation
            for i = 1 : nSamps
                % Shuffle feature space distance matrix
                if exact
                    perm = permInds(i,:);
                else
                    perm = randperm(obj.n);
                end
                % shuffle distance matrix
                dFtrDistMatShuf = dFtrDistMat(perm,perm);
                % compute sample
                samps(i) = sum(dFtrDistMatShuf(mask) .* dMapDists) .* normTerm;
            end
            
            % compute p-value
            if exact
                p = sum(samps >= c) ./ nExact;
            else
                p = (sum(samps >= c) + 1) ./ (nMC + 1);
            end
        end        
        
        function [c, p] = dcorr_p(obj, nMC)
            % Setup mask to mark valid pairs
            if obj.pairs
                pairs2 = [obj.pairs ; fliplr(obj.pairs)];
                mask = false(onj.n);
                mask(sub2ind([obj.n obj.n], pairs2(:,1), pairs2(:,2))) = true;
            else
                mask = true(obj.n);
            end
            
            % Counts of valid pairs
            colCnt = sum(mask, 1);
            rowCnt = sum(mask, 2);
            gndCnt = sum(colCnt);
            
            % Compute recentered map space distances
            dMap = obj.distMatrix(obj.map, obj.mapAngular);
            dMap(~mask) = 0;
            dMapR = bsxfun(@minus, dMap, sum(dMap, 1) ./ colCnt);
            dMapR = bsxfun(@minus, dMapR, sum(dMap, 2) ./ rowCnt);
            dMapR = dMapR + sum(dMap(:)) ./ gndCnt;

            % Compute recentered feature space distances
            dFtr = obj.distMatrix(obj.ftr, obj.ftrAngular);
            dFtr(~mask) = 0;
            dFtrR = bsxfun(@minus, dFtr, sum(dFtr, 1) ./ colCnt);
            dFtrR = bsxfun(@minus, dFtrR, sum(dFtr, 2) ./ rowCnt);
            dFtrR = dFtrR + sum(dFtr(:)) ./ gndCnt;
            
            % Compute dcorr
            dcov = sqrt(mean(dMapR(:) .* dFtrR(:)));
            dvarMap = sqrt(mean(dMapR(:) .* dMapR(:)));
            dvarFtr = sqrt(mean(dFtrR(:) .* dFtrR(:)));
            c = dcov ./ sqrt(dvarMap .* dvarFtr);

            % Is n small enough to do exact permutation?
            nExact = factorial(obj.n);
            exact = nExact < nMC;

            if exact
                nPerm = nExact;
                permInds = perms(1:obj.n);
            else
                nPerm = nMC;
            end

            % Do permutation analysis
            samps = zeros(nPerm,1);

            for i = 1 : nPerm
                if exact
                    % Get next permutation
                    ftrShuffled = obj.ftr(permInds(i,:),:);
                else
                    % Shuffle z values randomly
                    ftrShuffled = obj.ftr(randperm(obj.n),:);
                end
                
                % Compute new feature space distances
                dFtrShuffled = obj.distMatrix(ftrShuffled, obj.ftrAngular);
                dFtrShuffled(~mask) = 0;
                dFtrShuffledR = bsxfun(@minus, dFtrShuffled, sum(dFtrShuffled, 1) ./ colCnt);
                dFtrShuffledR = bsxfun(@minus, dFtrShuffledR, sum(dFtrShuffled, 2) ./ rowCnt);
                dFtrShuffledR = dFtrShuffledR + sum(dFtrShuffled(:)) ./ gndCnt;
                
                % Compute dcorr for sample
                dcov = sqrt(mean(dMapR(:) .* dFtrShuffledR(:)));
                dvarFtr = sqrt(mean(dFtrShuffledR(:) .* dFtrShuffledR(:)));
                samps(i) = dcov ./ sqrt(dvarMap .* dvarFtr);
            end

            % Compute p-value
            if exact
                p = sum(samps >= c) ./ nExact;
            else
                p = (sum(samps >= c) + 1) ./ (nMC + 1);
            end
        end
        
    end
    
    methods (Static)
        
        function funcs = measures(kind)
            switch kind
                case 'funcs'
                    funcs = {'pearson_pairs_p',...
                             'spearman_pairs_p',...
                             'zrehen_p',...
                             'wiringLength_p',...
                             'pathLength_p',...
                             'topoProd_p',...
                             'topoCorr_p',...
                             'dcorr_p'};
                case 'abbrv'
                    funcs = {'PC',...
                             'SC',...
                             'ZM',...
                             'WL',...
                             'PL',...
                             'TP',...
                             'TC',...
                             'DC'};
                case 'names'
                    funcs = {'Pearson distance correlation',...
                             'Spearman distance correlation',...
                             'Zrehen measure',...
                             'Wiring length',...
                             'Path length',...
                             'Topographic product',...
                             'Topological correlation',...
                             'Szekely distance correlation'};
                otherwise
                    error('valid arguments are: ''funcs'' ''abbrv'' ''names''')
            end
        end
        
    end
    
    methods (Access=protected)
        
        function [r, p] = corr_pairs_p(obj, meth, nMC)
            if obj.pairs
                prs = obj.pairs;
            else
                prs = nchoosek(1:obj.n, 2);
            end
            
            % Define correlation functions
            function rho = pearson(x1, x2)
                rho = corrcoef(x1, x2);
                rho = rho(1,2);
            end

            function rho = spearman(x1, x2)
                r1 = ranks(x1);
                r2 = ranks(x2);
                rho = corrcoef(r1, r2);
                rho = rho(1,2);
            end

            % Assign correlation function
            switch meth
                case 'pearson'
                    fcorr = @pearson;
                case 'spearman'
                    fcorr = @spearman;
                otherwise
                    error('Unknown correlation method: %s', meth)
            end

            % Compute distances and true correlation coeff
            dMap = obj.distEuclidPairs(obj.map, prs, obj.mapAngular);
            dFtr = obj.distEuclidPairs(obj.ftr, prs, obj.ftrAngular);
            r = fcorr(dMap, dFtr);

            % Is n small enough to do exact permutation?
            nExact = factorial(obj.n);

            if nExact > nMC
                % Do MC permutation
                samps = zeros(nMC,1);

                for i = 1 : nMC
                    % Shuffle feature space positions
                    ftrShuffled = obj.ftr(randperm(obj.n),:);
                    % Compute new feature space distances
                    dFtrShuffled = obj.distEuclidPairs(ftrShuffled, prs, obj.ftrAngular);
                    % Compute rho for sample
                    samps(i) = fcorr(dMap, dFtrShuffled);
                end

                % Compute p value
                p = (sum(samps >= r) + 1) ./ (nMC + 1);
            else
                % Do exact permutation
                samps = zeros(nExact, 1);
                permInds = perms(1:obj.n);
                
                for i = 1 : nExact
                    % Shuffle feature space positions
                    ftrShuffled = obj.ftr(permInds(i,:),:);
                    % Compute new feature space distances
                    dFtrShuffled = obj.distEuclidPairs(ftrShuffled, prs, obj.ftrAngular);
                    % Compute rho for sample
                    samps(i) = fcorr(dMap, dFtrShuffled);        
                end

                p = sum(samps >= r) ./ nExact;
            end
        end
        
    end
    
    methods (Access=protected,Static)
        
        function d = distEuclidPairs2(points, prs, circ)
            % Compute distances in linear dimensions
            d(:,~circ) = points(prs(:,1),~circ) - points(prs(:,2),~circ);
            % Compute distances in circular dimensions
            d(:,circ) = circ_dist(points(prs(:,1),circ), points(prs(:,2),circ));
            % Euclidean norm
            d = sum(d.^2, 2);
        end
        
        function d = distEuclidPairs(points, prs, circ)
            d = sqrt(MapData.distEuclidPairs2(points, prs, circ));
        end
        
        function d = distMatrix2(points, circ)
            % Compute distances in linear dimensions
            d(:,~circ,:) = bsxfun(@minus, points(:,~circ), permute(points(:,~circ), [3 2 1]));
            % Compute distances in circular dimensions
            d(:,circ,:) = bsxfun(@circ_dist, points(:,circ), permute(points(:,circ), [3 2 1]));
            % Euclidean norm
            d = permute(sum(d .* d, 2), [1 3 2]);
        end
        
        function d = distMatrix(points, circ)
            d = realsqrt(MapData.distMatrix2(points, circ));
        end
        
        function g = constructGraph(points, circ)
            nn = size(points, 1);
            dim = size(points, 2);
            pts = points;
            
            for d = 1 : dim
                if circ(d)
                    % for periodic dimensions, replicate point set in
                    % either direction
                    offset = zeros(1, dim);
                    offset(d) = 2 * pi;
                    upperRep = bsxfun(@plus, pts, offset);
                    lowerRep = bsxfun(@minus, pts, offset);
                    pts = [pts ; lowerRep ; upperRep]; %#ok<AGROW>
                end
            end
            
            % Construct triangulation
            dt = delaunayn(pts);
            edges = delaunayEdges(dt);
            
            % Resolve edges across periodic boundaries
            if any(circ)
                % Discard any edges that do not contain an original point
                nonOrig = all(edges > nn, 2);
                edges(nonOrig,:) = [];
                % Convert indices of replicate points back to original
                edges = mod(edges - 1, nn) + 1;
            end
            
            % Create graph representation
            g = sparse(edges(:,2), edges(:,1), ones(size(edges(:,2))), nn, nn, length(edges));
            % Make symmetric
            g = g + g';
            % Remove duplicate edges
            g = min(g, 1);
            % Convert to lower tri form, discarding self-edges
            g = tril(g, -1);
        end
        
        function e = neighbours(points, circ)
            g = MapData.constructGraph(points, circ);
            [e1, e2] = find(g);
            e = [e1 e2];
        end
        
        function sd = std(points, circ)
            sd = zeros(1, size(points, 2));
            % SD in lin dims
            sd(~circ) = std(points(:,~circ));
            % SD in circular dims
            sd(circ) = circ_std(points(:,circ));
        end
    end
end
