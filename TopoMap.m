classdef TopoMap < handle
    
    properties
        debug = true;
        mapsize = 500;
        circular = false;
        map = [];
        type = '';
    end
    
    methods
        
        function obj = TopoMap(type, scale)
            switch type
                case 'linear'
                    obj.type = type;
                    obj.circular = false;
                    
                    % Sample gradient
                    theta = rand * 2 * pi - pi;
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
                    
                    obj.map = angle( fftconv2(h, xr) + fftconv2(h, xi) * 1i );
                    
                case 'clust'
                    obj.type = type;
                    obj.circular = true;
                    
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
                    
                otherwise
                    error('Unsupported map type')
            end
        end
        
        
        function [x, y, z] = observe(obj, sampleSpacing, SNR, makePlots, publish)
            % number of points at which to observe map 
            nPoints = ceil(1 / sampleSpacing^2);
            
            % Draw quasi-random points
            qrng = haltonset(2, 'Skip', 1e3 + ceil(1e6 * rand), 'Leap', 1e2);
            qrng = scramble(qrng, 'RR2');
            obsPts = net(qrng, nPoints);
            
            % Snap to nearest pixel
            xi = ceil(obsPts(:,1) * size(obj.map,1));
            yi = ceil(obsPts(:,2) * size(obj.map,2));
            pos = 0 : 1/obj.mapsize : obj.mapsize;
            x = pos(xi);
            y = pos(yi);
            
            % Lookup z values
            zRaw = obj.map(sub2ind(size(obj.map), xi, yi))';
            
            % Add noise
            if obj.circular
                stdNoise = circ_std(obj.map(:)) ./ SNR;
            else
                stdNoise = std(obj.map(:)) ./ SNR;
            end
            
            z = zRaw + randn(size(zRaw)) * stdNoise;
            
            if obj.circular
                z = circ_dist(z, 0); % wrap to interval [-pi pi]
            end
            
            if makePlots 
                if obj.circular
                    cl = [-pi pi];
                    cmap = 'HSV';
                else
                    minVal = min([obj.map(:) ; z(:)]);
                    maxVal = max([obj.map(:) ; z(:)]);
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
                plot(xi, yi, 'w+', 'markersize', 10, 'linewidth', 1)
                set(gca, 'dataaspectratio', [1 1 1], 'ydir', 'normal')
                set(gca, 'xtick', [], 'ytick', [])
                title('Ground truth map')
                
                if publish
                    set(gca, 'ActivePositionProperty', 'OuterPosition')
                    set(gcf, 'Color', 'w')
                    set(gcf, 'Units', 'centimeters');
                    set(gcf, 'OuterPosition', [5 10 10 10]);
                    %shrinkfig(ax, 0.9)
                    
                    export_fig(sprintf('modelplots/mapModelA_%s', datestr(now, 30)), '-pdf')
                    system(sprintf('~/Scripts/pdf2eps modelplots/mapModelA_%s.pdf', datestr(now, 30)));
                    system(sprintf('rm modelplots/mapModelA_%s.pdf', datestr(now, 30)));
                    figure
                else
                    subplot(1,3,2)
                end
                
                scatter(x, y, 20, zRaw, 'filled')
                caxis(cl)
                colormap(cmap)
                xlim([-0.02 1.02])
                ylim([-0.02 1.02])
                set(gca, 'dataaspectratio', [1 1 1])
                set(gca, 'xtick', [], 'ytick', [])
                box on
                title('Observations (no noise)')
                
                if publish
                    set(gca, 'ActivePositionProperty', 'OuterPosition')
                    set(gcf, 'Color', 'w')
                    set(gcf, 'Units', 'centimeters');
                    set(gcf, 'OuterPosition', [5 10 10 10]);
                    %shrinkfig(ax, 0.9)

                    export_fig(sprintf('modelplots/mapModelB_%s', datestr(now, 30)), '-eps')
                    figure
                else
                    subplot(1,3,3)
                end
                
                scatter(x,y,20,z,'filled')
                caxis(cl)
                colormap(cmap)
                xlim([-0.02 1.02])
                ylim([-0.02 1.02])
                set(gca, 'dataaspectratio', [1 1 1])
                set(gca, 'xtick', [], 'ytick', [])
                box on
                title('Observations (incl. noise)')
                
                if publish
                    set(gca, 'ActivePositionProperty', 'OuterPosition')
                    set(gcf, 'Color', 'w')
                    set(gcf, 'Units', 'centimeters');
                    set(gcf, 'OuterPosition', [5 10 10 10]);
                    %shrinkfig(ax, 0.9)

                    export_fig(sprintf('modelplots/mapModelC_%s', datestr(now, 30)), '-eps')
                end
            end
        end
        
        
        function plot(obj)
            if obj.circular
                cl = [-pi pi];
                ticks = -pi : pi/2 : pi;
                cmap = 'HSV';
            else
                minVal = min(obj.map(:));
                maxVal = max(obj.map(:));
                cl = [minVal maxVal];
                %ticks = minVal : diff(cl)/4 : maxVal;
                cmap = 'jet';
            end
                
            figure                
            imagesc(obj.map', cl)
            set(gca, 'fontsize', 11)
            colormap(cmap)
            hold on
            set(gca, 'dataaspectratio', [1 1 1], 'ydir', 'normal')
            set(gca, 'xtick', [], 'ytick', [])
            %xlabel('x')
            %ylabel('y')
            %title('Ground truth map')
            
            cbax = colorbar('southoutside');
            
            if obj.circular
                set(cbax, 'xtick', ticks)
                set(cbax, 'xticklabel', {'-p' '-p/2' '0' 'p/2' 'p'}, 'fontname', 'Symbol', 'fontsize', 14)
            else
                set(cbax, 'fontname', 'Symbol', 'fontsize', 14)
            end
            
            set(gcf, 'Color', 'w')
            set(gca, 'ActivePositionProperty', 'OuterPosition')
            set(gcf, 'Units', 'centimeters');
            set(gcf, 'OuterPosition', [5 10 16 16]);


        end
                
    end
end