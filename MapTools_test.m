function MapTools_test(arg)

nMC = 50000;

switch arg
    case 'generate'
        maps{1} = TopoMap('orient', 0.4);
        maps{2} = TopoMap('clust', 0.6);
        maps{3} = TopoMap('linear');

        Ns = [8 30 60];
        SNRs = [0.5 1 2];

        for m = 1:3
            for n = 1:3
                for s = 1:3
                    data{m,n,s} = maps{m}.observe(Ns(n), SNRs(s), false, false, 'N');
                    x{m,n,s} = data{m,n,s}.map(:,1);
                    y{m,n,s} = data{m,n,s}.map(:,2);
                    z{m,n,s} = data{m,n,s}.ftr;
                    circ{m,n,s} = data{m,n,s}.ftrAngular;
                end
            end
        end
        
        filename = sprintf('testdata_%s', datestr(now, 'HHMM-ddmmYYYY'));
        save(filename)
        
    case 'test'
        load('testdata')
        
        oldMeasures = {@pearson_pairs_p,...
                       @spearman_pairs_p,...
                       @zrehen_p,...
                       @minWiring_p,...
                       @minPath_p,...
                       @topoProd_P,...
                       @topoCorr_pMC,...
                       @dcorr_p};
        
        newMeasures = MapData.measures('funcs');
        mName = MapData.measures('abbrv');
        
        for m = 1:3
            for n = 1:3
                for s = 1:3
                    %for q = 1:8
                    for q = 6
                        measure = oldMeasures{q};
                        fprintf('Map%d N=%d SNR=%.1f %s...', m, Ns(n), SNRs(s), mName{q})
                        
                        [cOld(m,n,s,q) pOld(m,n,s,q)] = measure(x{m,n,s},...
                                                                y{m,n,s},...
                                                                z{m,n,s},...
                                                                nMC,...
                                                                [],...
                                                                circ{m,n,s});
                        [cNew(m,n,s,q) pNew(m,n,s,q)] = data{m,n,s}.measure(q, nMC);
                        
                        cDiff(m,n,s,q) = cNew(m,n,s,q) - cOld(m,n,s,q);
                        pDiff(m,n,s,q) = pNew(m,n,s,q) - pOld(m,n,s,q);
                        cDiffR(m,n,s,q) = cDiff(m,n,s,q) ./ cOld(m,n,s,q);
                        pDiffR(m,n,s,q) = pDiff(m,n,s,q) ./ pOld(m,n,s,q);
                        
                        fprintf('c=%.2f/%.2f p=%.2g/%.2g\n', cNew(m,n,s,q), cOld(m,n,s,q), pNew(m,n,s,q), pOld(m,n,s,q))
                    end
                end
            end
        end
        
        fprintf('\n\nMeasures relative difference: mean=%f SD=%f\n', mean(cDiffR(:)), std(cDiffR(:)))
        fprintf('\n\nP-values relative difference: mean=%f SD=%f\n', mean(pDiffR(:)), std(pDiffR(:)))
        save('testresults');
        
    case 'testnew'
        maps{1} = TopoMap('linear2', pi/2);
        maps{2} = TopoMap('composite', pi/2, 0.6);
        Ns = [8 30];
        SNR = 1.5;
        
        mName = MapData.measures('abbrv');
        
        for m = 1:2
            for n = 1 : length(Ns)
                for q = 1:8
                    fprintf('Map%d N=%d SNR=%.1f %s...', m, Ns(n), SNR, mName{q})
                    data = maps{m}.observe(Ns(n), SNR, false, false, 'N');
                    
                    [cNew(m,q) pNew(m,q)] = data.measure(q, nMC);
                    
                    fprintf('c=%.2f p=%.2g\n', cNew(m,q), pNew(m,q))
                end
            end
        end
        
        save('newtestresult')
        
    otherwise
        error('valid arguments are: ''generate'' ''test'' ''testnew''')
end