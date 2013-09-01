classdef TopoMap < handle
    
    properties
        debug = true;
        mapsize = 500;
        circular = false;
        featureDims = 1;
        map = [];
        type = '';
    end
    
    methods
        
        function obj = TopoMap(type, varargin)
            switch type
                case 'linear'
                    obj.type = type;
                    obj.circular = false;
                    obj.featureDims = 1;
                    
                    % Sample gradient
                    %theta = rand * 2 * pi - pi;
                    theta = 0;
                    a = cos(theta);
                    b = sin(theta);
                    
                    % Generate map space
                    x = -0.5 : 1/obj.mapsize : 0.5;
                    
                    % Generate map
                    fmap = @(x1, x2) a .* x1 + b .* x2;
                    obj.map = bsxfun(fmap, x, x');
                                        
                case 'orient'
                    obj.type = type;
                    obj.circular = true;
                    obj.featureDims = 1;
                    assert(nargin == 2, 'args should be of the form: (''orient'', scale)')
                    scale = varargin{1};
                    
                    % Define filter coordinate space
                    spacing = 1/400;
                    fX = -3*scale : spacing : 3*scale;
                    
                    % Generate difference of gaussians filter
                    fXnorm2 = bsxfun(@(x1, x2) x1.^2 + x2.^2, fX, fX');
                    h = 1 / (2 * pi * scale^2) * exp(-fXnorm2 / (2 * scale^2)) - 2 / (pi * scale^2) * exp(-2 * fXnorm2 / (scale^2));
                    
                    % Generate noise, allowing additional area for filtering
                    noiseSize = obj.mapsize + length(fX);
                    
                    % Filter noise and compute Arg()                    
                    xi = randn(noiseSize);
                    xr = randn(noiseSize);
                    
                    obj.map = angle( TopoMap.fftconv2(h, xr) + TopoMap.fftconv2(h, xi) * 1i );
                    
                case 'clust'
                    obj.type = type;
                    obj.circular = true;
                    obj.featureDims = 1;
                    assert(nargin == 2, 'args should be of the form: (''clust'', scale)')
                    scale = varargin{1};
                    
                    % Max density (40 x 40 seed points)
                    seedDens = 40;
                    
                    % Generate seed points
                    qrng = haltonset(2, 'Skip', 1e3 + ceil(1e6 * rand), 'Leap', 1e2);
                    qrng = scramble(qrng, 'RR2');
                    seeds = net(qrng, seedDens^2);
                    
                    % Re-centre on origin and scale
                    seeds = seeds - 0.5;
                    %seeds = (scale / 0.031) * seeds;
                    seeds = scale * seedDens * seeds;
                    
                    % Sample z-values
                    zSeeds = rand(size(seeds, 1), 1) * 2 * pi - pi;
                    
                    % Construct triangulation
                    tri = delaunay(seeds(:,1), seeds(:,2));
                    
                    % Generate map
                    x = -0.5 : 1/obj.mapsize : 0.5;
                    x = repmat(x, [length(x) 1]);
                    y = x';
                    z = zSeeds(dsearchn(seeds, tri, [x(:) y(:)]));                    
                    obj.map = reshape(z, size(x));
                
                case 'linear2'
                    obj.type = type;
                    obj.circular = [false false];
                    obj.featureDims = 2;
                    assert(nargin == 2, 'args should be of the form: (''linear2'', theta)')
                    
                    % Set gradients
                    % scale arg is angle difference between feature dim gradients
                    theta = [0 varargin{1}];
                    a = cos(theta);
                    b = sin(theta);
                    
                    % Generate map space
                    x = -0.5 : 1/obj.mapsize : 0.5;
                    
                    % Generate map
                    fmap1 = @(x1, x2) a(1) .* x1 + b(1) .* x2;
                    fmap2 = @(x1, x2) a(2) .* x1 + b(2) .* x2;
                    obj.map = cat(3, bsxfun(fmap1, x, x'), bsxfun(fmap2, x, x'));
                    
                case 'composite'
                    obj.type = type;
                    obj.circular = [false false true];
                    obj.featureDims = 3;
                    assert(nargin == 3, 'args should be of the form: (''composite'', theta, scale)')
                    
                    theta = [0 varargin{1}];
                    a = cos(theta);
                    b = sin(theta);
                    scale = varargin{2};
                    
                    % Generate map space
                    x = -0.5 : 1/obj.mapsize : 0.5;
                    
                    % Generate linear map dims
                    fmap1 = @(x1, x2) a(1) .* x1 + b(1) .* x2;
                    fmap2 = @(x1, x2) a(2) .* x1 + b(2) .* x2;
                    
                    % Define filter coordinate space
                    spacing = 1/400;
                    fX = -3*scale : spacing : 3*scale;
                    
                    % Generate difference of gaussians filter
                    fXnorm2 = bsxfun(@(x1, x2) x1.^2 + x2.^2, fX, fX');
                    h = 1 / (2 * pi * scale^2) * exp(-fXnorm2 / (2 * scale^2)) - 2 / (pi * scale^2) * exp(-2 * fXnorm2 / (scale^2));
                    
                    % Generate noise, allowing additional area for filtering
                    noiseSize = obj.mapsize + length(fX);
                    
                    % Filter noise and compute Arg()                    
                    xi = randn(noiseSize);
                    xr = randn(noiseSize);
                    
                    obj.map = cat(3,...
                                  bsxfun(fmap1, x, x'),...
                                  bsxfun(fmap2, x, x'),...
                                  angle( TopoMap.fftconv2(h, xr) + TopoMap.fftconv2(h, xi) * 1i ));
                    
                    % renormalise linear dims to same range as circular dim
                    obj.map(obj.map == -1) = 1;
                    obj.map(:,:,[1 2]) = pi .* obj.map(:,:,[1 2]);
                    
                otherwise
                    error('Unsupported map type')
            end
        end
        
        function map = observe(obj, sampleSpacing, SNR, makePlots, publish, varargin)
            if length(varargin) >= 1
                option = varargin{1};
                
                if option == 'N'
                    nPoints = sampleSpacing;
                end
            end
            
            % number of points at which to observe map 
            if ~exist('nPoints', 'var')
                nPoints = ceil(1 / sampleSpacing^2);
            end
            
            % Draw quasi-random points on unit disc
            qrng = haltonset(2, 'Skip', 1e3 + ceil(1e6 * rand), 'Leap', 1e2);
            qrng = scramble(qrng, 'RR2');
            overSample = 2;
            
            while true
                obsPts = net(qrng, nPoints * overSample);
                onUnitDisc = sum((obsPts - 0.5).^2, 2) <= 0.25;
                obsPts = obsPts(onUnitDisc,:);

                if size(obsPts, 1) >= nPoints
                    obsPts = obsPts(1:nPoints,:);
                    break
                else
                    overSample = overSample + 1;
                end
            end
            
            % Snap to nearest pixel
            xi = ceil(obsPts(:,1) * size(obj.map,1));
            yi = ceil(obsPts(:,2) * size(obj.map,2));
            pos = 0 : 1/obj.mapsize : obj.mapsize;
            x = pos(xi)';
            y = pos(yi)';
            xMap = [x y];
            
            for fDim = 1 : obj.featureDims
                % Lookup z values
                zRaw = obj.map(sub2ind(size(obj.map), xi, yi, repmat(fDim, size(xi))))';
                
                % Add noise
                subMap = obj.map(:,:,fDim);
                
                if obj.circular(fDim)
                    stdNoise = circ_std(subMap(:)) ./ SNR;
                else
                    stdNoise = std(subMap(:)) ./ SNR;
                end
                
                xFeature(:,fDim) = zRaw + randn(size(zRaw)) * stdNoise; %#ok<AGROW>
                
                if obj.circular(fDim)
                    xFeature(:,fDim) = circ_dist(xFeature(:,fDim), 0); %#ok<AGROW> % wrap to interval [-pi pi]
                end
            end
            
            % Construct MapData object to return
            map = MapData(xMap, xFeature, [false false], obj.circular);
            
            if makePlots && obj.featureDims == 1
                if obj.circular
                    cl = [-pi pi];
                    cmap = 'HSV';
                else
                    minVal = min([obj.map(:) ; xFeature(:)]);
                    maxVal = max([obj.map(:) ; xFeature(:)]);
                    cl = [minVal maxVal];
                    cmap = 'jet';
                end
                
                figure
                
                if ~publish
                    subplot(1,3,1)
                end
                
                imagesc(obj.map', cl)
                colormap(cmap)
                hold on
                plot(xi, yi, 'w+', 'markersize', 15, 'linewidth', 2)
                set(gca, 'dataaspectratio', [1 1 1], 'ydir', 'normal')
                set(gca, 'xtick', [], 'ytick', [])
                title('Ground truth map')
                
                if publish
                    set(gca, 'ActivePositionProperty', 'OuterPosition')
                    set(gcf, 'Color', 'w')
                    set(gcf, 'Units', 'centimeters');
                    set(gcf, 'OuterPosition', [5 10 16 16]);
                    %shrinkfig(ax, 0.9)
                    
                    export_fig(sprintf('modelplots/mapModelA_%s', datestr(now, 30)), '-pdf', '-painters')
                    %system(sprintf('~/Scripts/pdf2eps modelplots/mapModelA_%s.pdf', datestr(now, 30)));
                    %system(sprintf('rm modelplots/mapModelA_%s.pdf', datestr(now, 30)));
                    figure
                else
                    subplot(1,3,2)
                end
                
                scatter(x, y, 50, zRaw, 'filled')
                caxis(cl)
                colormap(cmap)
                xlim([-0.05 1.05])
                ylim([-0.05 1.05])
                set(gca, 'dataaspectratio', [1 1 1])
                set(gca, 'xtick', [], 'ytick', [])
                box on
                title('Observations (no noise)')
                
                if publish
                    set(gca, 'ActivePositionProperty', 'OuterPosition')
                    set(gcf, 'Color', 'w')
                    set(gcf, 'Units', 'centimeters');
                    set(gcf, 'OuterPosition', [5 10 16 16]);
                    %shrinkfig(ax, 0.9)

                    export_fig(sprintf('modelplots/mapModelB_%s', datestr(now, 30)), '-pdf', '-painters')
                    figure
                else
                    subplot(1,3,3)
                end
                
                scatter(x, y, 50, xFeature, 'filled')
                caxis(cl)
                colormap(cmap)
                xlim([-0.05 1.05])
                ylim([-0.05 1.05])
                set(gca, 'dataaspectratio', [1 1 1])
                set(gca, 'xtick', [], 'ytick', [])
                box on
                title('Observations (incl. noise)')
                
                if publish
                    set(gca, 'ActivePositionProperty', 'OuterPosition')
                    set(gcf, 'Color', 'w')
                    set(gcf, 'Units', 'centimeters');
                    set(gcf, 'OuterPosition', [5 10 16 16]);
                    %shrinkfig(ax, 0.9)

                    export_fig(sprintf('modelplots/mapModelC_%s', datestr(now, 30)), '-pdf', '-painters')
                end
            end
        end
        
        function plot(obj)
            figure
            
            for d = 1 : obj.featureDims
                if obj.circular(d)
                    cl = [-pi pi];
                else
                    minVal = min(min(obj.map(:,:,d)));
                    maxVal = max(max(obj.map(:,:,d)));
                    cl = [minVal maxVal];
                    %ticks = minVal : diff(cl)/4 : maxVal;
                end
                
                subplot(1,obj.featureDims,d)
                imagesc(obj.map(:,:,d)', cl)
                set(gca, 'fontsize', 11)
                hold on
                set(gca, 'dataaspectratio', [1 1 1], 'ydir', 'normal')
                set(gca, 'xtick', [], 'ytick', [])
                %xlabel('x')
                %ylabel('y')
                %title('Ground truth map')

                cbax = colorbar('southoutside');

                if obj.circular(d)
                    colormap('HSV')
                    set(cbax, 'xtick', [-pi -pi/2 0 pi/2 pi])
                    set(cbax, 'xticklabel', {'-p' '-p/2' '0' 'p/2' 'p'}, 'fontname', 'Symbol', 'fontsize', 14)
                else
                    colormap('jet')
                    set(cbax, 'fontname', 'Symbol', 'fontsize', 14)
                end
                
                if exist('cbfreeze', 'file') == 2
                    cbfreeze
                end
                
                if exist('freezeColors', 'file') == 2
                    freezeColors
                end
            end
            
            set(gcf, 'Color', 'w')
            set(gca, 'ActivePositionProperty', 'OuterPosition')
            set(gcf, 'Units', 'centimeters');
            set(gcf, 'OuterPosition', [5 10 16 16]);
        end
        
        function plotInset(obj, pos)
            if obj.featureDims ~= 1, return; end
            
            if obj.circular
                cl = [-pi pi];
                %cmap = [0 : 1/63 : 1, 62/63 : -1/63 : 0]';
                %cmap = [cmap cmap cmap];
                cmap = 'HSV';
            else
                minVal = min(obj.map(:));
                maxVal = max(obj.map(:));
                cl = [minVal maxVal];
                cmap = 'gray';
            end
            
            ax = axes('position', pos);
            imagesc(obj.map', cl)
            colormap(cmap)
            box on
            set(ax, 'ydir', 'normal', 'xtick', [], 'ytick', [])
        end
        
    end
    
    methods(Static, Access=protected)
        
        function ret = fftconv2(h, x)
            % Work out dimensions and padding
            fftSize = max([size(h) ; size(x)]);
            outSize = abs(size(h) - size(x)) + 1;
            pad = (fftSize - outSize) ./ 2;
            % Forward transform filter and data
            H = fft2(h, fftSize(1), fftSize(2));
            X = fft2(x, fftSize(1), fftSize(2));
            % Do convolution
            HX = H .* X;
            % Inverse transform result
            hx = real(ifft2(HX));
            % Return useful region of result
            ret = hx(pad(1)+1 : end-pad(1), pad(2)+1 : end-pad(2));
        end
        
    end
end