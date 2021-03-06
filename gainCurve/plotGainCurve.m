function plotGainCurve(inputFn, ext, plotSubthreshold, plotInput, plotDendV, sizeSize)
    textFontSize = 6;
    if nargin < 6
        sizeSize = 'int64';
        if nargin < 5
            plotDendV = false;
            if nargin < 4
                plotInput = true;
                if nargin < 3
                    plotSubthreshold = true;
                    if nargin < 2
                        ext = '';
                    end
                end
            end
        end
    end
    format = ext;
    if isequal(ext, 'psc')
        format = 'epsc';
    end
    if isequal(ext, 'jpg')
        format = 'jpeg';
    end
    p = read_cfg(inputFn);
    dimsFn = 'readoutDimension.bin';
    if isfield(p,'reformatInputFn')
        dimsFn = p.reformatInputFn;
    end
    if ~isstruct(p)
        return
    end
    fdr = struct2cell(dir());
    % name, fdr, date, bytes, isdir, datenum
    filter = [1,5];
    fdr = fdr(filter,:);
    fnstr = strjoin(fdr(1,:))
    load(p.libFile,'tstep','n','gList','dtRange','sPSP','loc');
    ldur = size(sPSP,1) - p.ignoreT;

    nE = 0;
    nI = 0;
    for i=1:n
        if gList(i) > 0
            nE = nE + 1;
        else
            nI = nI + 1;
        end
    end
    datafilePattern = ['(Data|jND|Raster|tIn|cpuTime)-',p.theme,'\-s\d*\.bin'];
    datafile = regexp(fnstr, datafilePattern,'match')
    assert(length(datafile)==5);
    dimsFid = fopen(dimsFn,'r');
    if dimsFid
        nTrial = fread(dimsFid, 1, 'int');
        nDimSim = fread(dimsFid, nTrial, 'uint64')
        nDim = fread(dimsFid, nTrial, 'uint64')
        dt = fread(dimsFid, nTrial, 'double')
        inputLevel = fread(dimsFid, nTrial, 'double')
        runTime = fread(dimsFid, nTrial, 'double')
        fclose(dimsFid);
    end
    RasterFn = datafile{2};
    RasterFid = fopen(RasterFn,'r');
    ss = zeros(nTrial,1);
    ts = cell(nTrial,1);
    sb = zeros(nTrial,1);
    tb = cell(nTrial,1);
    sl = zeros(nTrial,1);
    tl = cell(nTrial,1);
    sjb = zeros(nTrial,1);
    tjb = cell(nTrial,1);
    sjl = zeros(nTrial,1);
    tjl = cell(nTrial,1);
    sb0 = zeros(nTrial,1);
    tb0 = cell(nTrial,1);
    sl0 = zeros(nTrial,1);
    tl0 = cell(nTrial,1);
    cpuFn = datafile{3};
    cpuFid = fopen(cpuFn,'r');
    if p.one
        cpuTime = fread(cpuFid,[nTrial,7],'double');
        cpuTime = cpuTime';
    else
        cpuTime = fread(cpuFid,[7,nTrial],'double');
    end
    figure;
    if p.one
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            ss(i) = fread(RasterFid,1,sizeSize);
            if ss(i) >0
                ts{i} = fread(RasterFid,ss(i),'double');
            end
            plot(ts{i},zeros(ss(i),1)+1,'.k','MarkerSize',10);
        end
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            sb(i) = fread(RasterFid,1,sizeSize);
            if sb(i) >0
                tb{i} = fread(RasterFid,sb(i),'double');
            end
            plot(tb{i},zeros(sb(i),1)+2,'.b','MarkerSize',10);
        end
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            sl(i) = fread(RasterFid,1,sizeSize);
            if sl(i) >0
                tl{i} = fread(RasterFid,sl(i),'double');
            end
            plot(tl{i},zeros(sl(i),1)+3,'.r','MarkerSize',10);
        end
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            sjb(i) = fread(RasterFid,1,sizeSize);
            if sjb(i) >0
                tjb{i} = fread(RasterFid,sjb(i),'double');
            end
            plot(tjb{i},zeros(sjb(i),1)+4,'ob','MarkerSize',3);
        end
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            sjl(i) = fread(RasterFid,1,sizeSize);
            if sjl(i) >0
                tjl{i} = fread(RasterFid,sjl(i),'double');
            end
            plot(tjl{i},zeros(sjl(i),1)+5,'or','MarkerSize',3);
        end
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            sb0(i) = fread(RasterFid,1,sizeSize);
            if sb0(i) >0
                tb0{i} = fread(RasterFid,sb0(i),'double');
            end
            plot(tb0{i},zeros(sb0(i),1)+6,'*b','MarkerSize',3);
        end
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            sl0(i) = fread(RasterFid,1,sizeSize);
            if sl0(i) >0
                tl0{i} = fread(RasterFid,sl0(i),'double');
            end
            plot(tl0{i},zeros(sl0(i),1)+7,'*r','MarkerSize',3);
        end
    else 
        for i = 1:nTrial
            subplot(nTrial,2,(i-1)*2+1)
            hold on
            ss(i) = fread(RasterFid,1,sizeSize);
            if ss(i) >0
                ts{i} = fread(RasterFid,ss(i),'double');
            end
            plot(ts{i},zeros(ss(i),1)+1,'.k','MarkerSize',10);

            sb(i) = fread(RasterFid,1,sizeSize);
            if sb(i) >0
                tb{i} = fread(RasterFid,sb(i),'double');
            end
            plot(tb{i},zeros(sb(i),1)+2,'.b','MarkerSize',10);

            sl(i) = fread(RasterFid,1,sizeSize);
            if sl(i) >0
                tl{i} = fread(RasterFid,sl(i),'double');
            end
            plot(tl{i},zeros(sl(i),1)+3,'.r','MarkerSize',10);

            sjb(i) = fread(RasterFid,1,sizeSize);
            if sjb(i) >0
                tjb{i} = fread(RasterFid,sjb(i),'double');
            end
            plot(tjb{i},zeros(sjb(i),1)+4,'ob','MarkerSize',3);

            sjl(i) = fread(RasterFid,1,sizeSize);
            if sjl(i) >0
                tjl{i} = fread(RasterFid,sjl(i),'double');
            end
            plot(tjl{i},zeros(sjl(i),1)+5,'or','MarkerSize',3);

            sb0(i) = fread(RasterFid,1,sizeSize);
            if sb0(i) >0
                tb0{i} = fread(RasterFid,sb0(i),'double');
            end
            plot(tb0{i},zeros(sb0(i),1)+6,'*b','MarkerSize',10);

            sl0(i) = fread(RasterFid,1,sizeSize);
            if sl0(i) >0
                tl0{i} = fread(RasterFid,sl0(i),'double');
            end
            plot(tl0{i},zeros(sl0(i),1)+7,'*r','MarkerSize',10);
        end
    end
    save(['../', RasterFn(1:length(RasterFn)-4),'.mat'], 'ts', 'tb', 'tl', 'tjb', 'tjl', 'tb0','tl0');
    fclose(RasterFid);
    subplot(2,2,2)
    hold on
    plot(inputLevel,ss./runTime*1000,'-*k');
    plot(inputLevel,sb./runTime*1000,'-*b');
    plot(inputLevel,sl./runTime*1000,'-*r');
    plot(inputLevel,sjb./runTime*1000,'-oc');
    plot(inputLevel,sjl./runTime*1000,'-om');
    plot(inputLevel,sb0./runTime*1000,'-*g');
    plot(inputLevel,sl0./runTime*1000,'-*y');
    xlim([0,inputLevel(nTrial)*1.1]);
    xlabel('input rate Hz');
    ylabel('firing rate Hz');

    subplot(2,2,4)
    hold on
    plot(inputLevel,cpuTime(1,:),'-*k');
    plot(inputLevel,cpuTime(2,:),'-*b');
    plot(inputLevel,cpuTime(3,:),'-*r');
    plot(inputLevel,cpuTime(4,:),'-oc');
    plot(inputLevel,cpuTime(5,:),'-om');
    plot(inputLevel,cpuTime(6,:),'-*g');
    plot(inputLevel,cpuTime(7,:),'-*y');
    xlim([0,inputLevel(nTrial)*1.1]);
    legend({'sim','bilinear','linear','jb','jl','bilinear0','linear0'});
    xlabel('input rate Hz');
    ylabel('cpuTime s');

    if ~isempty(ext)
        saveas(gcf,[p.theme,'-gainCurve.',ext],format);
    end
    if plotSubthreshold
        dataFn = datafile{1};
        jNDFn = datafile{4};
        tInFn = datafile{5};
        jbt = cell(nTrial,1);
        jbv = cell(nTrial,1);
        jlt = cell(nTrial,1);
        jlv = cell(nTrial,1);
        dendV = cell(nTrial,n);
        simV = cell(nTrial,1);
        biV = cell(nTrial,1);
        biV0 = cell(nTrial,1);
        liV0 = cell(nTrial,1);
        datafid = fopen(dataFn,'r');
        jNDfid = fopen(jNDFn,'r');
        if plotInput
            tInfid = fopen(tInFn,'r');
        end
        disp(tstep);
        if p.one
            if datafid
                for i=1:nTrial
                    simV{i} = fread(datafid, nDimSim(i), 'double');
                    if p.getDendV
                        for j=1:n
                            dendV{i,j} = fread(datafid, nDimSim(i), 'double');
                        end
                    end
                end
                for i=1:nTrial
                    biV{i} = fread(datafid, nDim(i), 'double');
                end
                for i=1:nTrial
                    liV{i} = fread(datafid, nDim(i), 'double');
                end
                for i=1:nTrial
                    biV0{i} = fread(datafid, nDim(i), 'double');
                end
                for i=1:nTrial
                    liV0{i} = fread(datafid, nDim(i), 'double');
                end
            end
            if jNDfid
                jbCrossT = cell(nTrial,1);
                jbCrossV = cell(nTrial,1);
                jbnCross = zeros(nTrial,1);
                jlCrossT = cell(nTrial,1);
                jlCrossV = cell(nTrial,1);
                jlnCross = zeros(nTrial,1);
                for i=1:nTrial
                    jbSize = fread(jNDfid,1,sizeSize);
                    jbt{i} = fread(jNDfid,jbSize,'double')*tstep;
                    jbv{i} = fread(jNDfid,jbSize,'double');
                    jbnCross(i) = fread(jNDfid,1,'int')
                    jbCrossT{i} = cell(jbnCross(i),1);
                    jbCrossV{i} = cell(jbnCross(i),1);
                    for j=1:jbnCross(i)
                        tmpSize = fread(jNDfid, 1, sizeSize)
                        jbCrossT{i}{j} = fread(jNDfid, tmpSize, 'double')*tstep;
                        jbCrossV{i}{j} = fread(jNDfid, tmpSize, 'double');
                    end
                end
                for i=1:nTrial
                    jlSize = fread(jNDfid,1,sizeSize);
                    jlt{i} = fread(jNDfid,jlSize,'double')*tstep;
                    jlv{i} = fread(jNDfid,jlSize,'double');
                    jlnCross(i) = fread(jNDfid,1,'int')
                    jlCrossT{i} = cell(jlnCross(i));
                    jlCrossV{i} = cell(jlnCross(i));
                    for j=1:jlnCross(i)
                        tmpSize = fread(jNDfid, 1, sizeSize)
                        jlCrossT{i}{j} = fread(jNDfid, tmpSize, 'double')*tstep;
                        jlCrossV{i}{j} = fread(jNDfid, tmpSize, 'double');
                    end
                end
            end
        else
            jbCrossT = cell(nTrial,1);
            jbCrossV = cell(nTrial,1);
            jbnCross = zeros(nTrial,1);
            jlCrossT = cell(nTrial,1);
            jlCrossV = cell(nTrial,1);
            jlnCross = zeros(nTrial,1);
            for i=1:nTrial
                if datafid
                    simV{i} = fread(datafid, nDimSim(i), 'double');
                    biV{i} = fread(datafid, nDim(i), 'double');
                    liV{i} = fread(datafid, nDim(i), 'double');
                    biV0{i} = fread(datafid, nDim(i), 'double');
                    liV0{i} = fread(datafid, nDim(i), 'double');
                    if p.getDendV
                        for j=1:n
                            dendV{i,j} = fread(datafid, nDimSim(i), 'double');
                        end
                    end
                end
                if jNDfid
                    jbSize = fread(jNDfid,1,sizeSize);
                    jbt{i} = fread(jNDfid,jbSize,'double')*tstep;
                    jbv{i} = fread(jNDfid,jbSize,'double');
                    jbnCross(i) = fread(jNDfid,1,'int')
                    jbCrossT{i} = cell(jbnCross(i),1);
                    jbCrossV{i} = cell(jbnCross(i),1);
                    for j=1:jbnCross(i)
                        tmpSize = fread(jNDfid, 1, sizeSize)
                        jbCrossT{i}{j} = fread(jNDfid, tmpSize, 'double')*tstep;
                        jbCrossV{i}{j} = fread(jNDfid, tmpSize, 'double');
                    end

                    jlSize = fread(jNDfid,1,sizeSize);
                    jlt{i} = fread(jNDfid,jlSize,'double')*tstep;
                    jlv{i} = fread(jNDfid,jlSize,'double');
                    jlnCross(i) = fread(jNDfid,1,'int')
                    jlCrossT{i} = cell(jlnCross(i),1);
                    jlCrossV{i} = cell(jlnCross(i),1);
                    for j=1:jlnCross(i)
                        tmpSize = fread(jNDfid, 1, sizeSize)
                        jlCrossT{i}{j} = fread(jNDfid, tmpSize, 'double')*tstep;
                        jlCrossV{i}{j} = fread(jNDfid, tmpSize, 'double');
                    end
                end
            end
        end
        for i=1:nTrial
            run_nt = round(runTime(i)/tstep) + 1;
            disp(['this is ', num2str(i), 'th trial']);
            disp(dt(i));
            t0 = linspace(0,nDimSim(i)-1,nDimSim(i)) * dt(i);
            t = linspace(0,nDim(i)-1,nDim(i)) * tstep;
            assert(t(end)==t0(end));
            xl = [0, t(end)];

            figure;
            subplot(2,1,1)
            hold on
            hs = plot(t0,simV{i},'k');
            handles = [hs];
            labels = {'sim'};
            if plotDendV && p.getDendV
                dline = [];
                sat = 0.8;
                val = 0.8;
                for j=1:n
                    hue = (double(j)-1.0)/(double(n)-1.0)*(2.0/3.0);
                    size(t0)
                    size(dendV{i,j})
                    dline(j) = plot(t0,dendV{i,j},'--','Color',hsv2rgb([hue,sat,val]));
                end
            end
            hb = plot(t,biV{i},'b');
            hl = plot(t,liV{i},'r');
            hb0 = plot(t,biV0{i},'g');
            hl0 = plot(t,liV0{i},'y');
            handles = [handles, hb, hl, hb0, hl0];
            labels = [labels, 'bi','li','bi0','hl0'];
            minV0 = min([min(simV{i}),min(biV{i}),min(liV{i}),min(biV0{i}),min(liV0{i})]);
            maxV0 = min([-50,max([max(simV{i}),max(biV{i}),max(liV{i}),max(biV0{i}),max(liV0{i})])]);
            vStretch = maxV0-minV0;
            minV = minV0-vStretch*0.1;
            maxV = maxV0+vStretch*0.1;
            yl = [minV,maxV];
            ylim(yl);
            if jNDfid
                hjb = plot(jbt{i},jbv{i},'.b','MarkerSize',3);
                for j=1:jbnCross(i)
                    plot(jbCrossT{i}{j}, jbCrossV{i}{j},':b');
                end
                hjl = plot(jlt{i},jlv{i},'.r','MarkerSize',3);
                for j=1:jlnCross(i)
                    plot(jlCrossT{i}{j}, jlCrossV{i}{j},':r');
                end
                handles = [handles, hjb, hjl];
                labels = [labels, 'jb','jl'];
            end
            if plotDendV && p.getDendV
                handles = [handles, dline];
                labels = [labels, strcat('dend[',cellstr(num2str(loc))',']')];
            end
            legend(handles,labels);
            if plotInput
                ntmp = fread(tInfid,[1,1],sizeSize);
                tin = fread(tInfid,[ntmp,1],'double');
                tID = fread(tInfid,[ntmp,1],sizeSize);
                tID = tID + 1;
                pE = (tID <= nE);
                pI = (tID > nE);
                tE = tin(pE);
                tI = tin(pI);
                Eid = tID(pE);
                Iid = tID(pI)-nE;
                Ein = length(tE);
                Iin = length(tI);
                vEtar = zeros(Ein,1);
                plot(tE,minV*ones(1,Ein),'.r');
                vItar = zeros(Iin,1);
                plot(tI,minV*ones(1,Iin),'.b');
                for j=1:Ein
                    iv = ceil(tE(j)/tstep);
                    if iv==0, iv=iv+1,end
                    jv = iv + 1;
                    vEtar = simV{i}(iv) + mod(tE(j),tstep)/tstep * (simV{i}(jv)-simV{i}(iv));
                    plot([tE(j),tE(j)],[minV,vEtar],':r');
                    text(tE(j),minV+(vEtar-minV)*0.3,num2str(Eid(j)),'Color','r','FontSize',textFontSize);
                    text(tE(j),minV+(vEtar-minV)*0.1,num2str(j),'Color','r','FontSize',textFontSize);
                 
                    iv = ceil((tE(j)+ldur)/tstep);
                    if iv==0, iv=iv+1,end
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vEtar = simV{i}(iv) + mod(tE(j),tstep)/tstep * (simV{i}(jv)-simV{i}(iv));
                        plot([tE(j),tE(j)]+ldur,[vEtar,maxV],':r');
                        text(tE(j)+ldur,maxV-(maxV-vEtar)*0.3,num2str(Eid(j)),'Color','r','FontSize',textFontSize);
                        text(tE(j)+ldur,maxV-(maxV-vEtar)*0.1,num2str(j),'Color','r','FontSize',textFontSize);
                    end
                end
                for j=1:Iin
                    iv = ceil(tI(j)/tstep);
                    if iv==0, iv=iv+1,end
                    jv = iv + 1;
                    vItar(j) = simV{i}(iv) + mod(tI(j),tstep)/tstep * (simV{i}(jv)-simV{i}(iv));
                    plot([tI(j),tI(j)],[minV,vItar(j)],':b');
                    text(tI(j),minV+(vItar(j)-minV)*0.3,num2str(Iid(j)),'Color','b','FontSize',textFontSize);
                    text(tI(j),minV+(vItar(j)-minV)*0.1,num2str(j),'Color','b','FontSize',textFontSize);
                        
                    iv = ceil((tI(j)+ldur)/tstep);
                    if iv==0, iv=iv+1,end
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vtar = simV{i}(iv) + mod(tI(j),tstep)/tstep * (simV{i}(jv)-simV{i}(iv));
                        plot([tI(j),tI(j)]+ldur,[vtar,maxV],':b');
                        text(tI(j)+ldur,maxV-(maxV-vtar)*0.3,num2str(Iid(j)),'Color','b','FontSize',textFontSize);
                        text(tI(j)+ldur,maxV-(maxV-vtar)*0.1,num2str(j),'Color','b','FontSize',textFontSize);
                    end
                end
            end
            xlim(xl);

            if dt(i) == tstep
                subplot(2,2,3)
            else 
                subplot(2,1,2)
            end
            hold on
            lh = [];
            lt = {};
            if ss(i) > 0
                ltmp = plot(ts{i},ones(1,ss(i)),'.k','MarkerSize',6,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'sims'];
            end
            if sb(i) > 0
                ltmp = plot(tb{i},ones(1,sb(i))+1,'.b','MarkerSize',6,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'bi'];
            end
            if sl(i) > 0
                ltmp = plot(tl{i},ones(1,sl(i))+2,'.r','MarkerSize',6,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'li'];
            end
            if sjb(i) > 0
                ltmp = plot(tjb{i},ones(1,sjb(i))+3,'ob','MarkerSize',3,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'jb'];
            end
            if sjl(i) > 0
                ltmp = plot(tjl{i},ones(1,sjl(i))+4,'or','MarkerSize',3,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'jl'];
            end
            if sb0(i) > 0
                ltmp = plot(tb0{i},ones(1,sb0(i))+5,'*b','MarkerSize',3,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'b0'];
            end
            if sl0(i) > 0
                ltmp = plot(tl0{i},ones(1,sl0(i))+6,'*r','MarkerSize',3,'LineStyle','none');
                lh = [lh, ltmp];
                lt = [lt,'l0'];
            end
            legend(lh,lt);
            ylim([0,8]);
            xlim(xl);

            if dt(i) == tstep
                subplot(2,2,4)
                hold on
                berr = biV{i} - simV{i};
                lerr = liV{i} - simV{i};
                b0err = biV0{i} - simV{i};
                l0err = liV0{i} - simV{i};
                hberr = plot(t0, berr, 'b');
                hlerr = plot(t0, lerr,'r');
                hb0err = plot(t0, b0err, 'g');
                hl0err = plot(t0, l0err, 'y');
                errorbar(t0(end)+1, mean(abs(berr)), std(berr),'*b');
                errorbar(t0(end)+2, mean(abs(lerr)), std(lerr),'*r');
                errorbar(t0(end)+3, mean(abs(b0err)), std(b0err),'*g');
                errorbar(t0(end)+4, mean(abs(l0err)), std(l0err),'*y');
                xlim(xl);
            end

            if ~isempty(ext)
                saveas(gcf,[p.theme,'-trial',num2str(i),'.',ext],format);
            end
        end
    end
    figure;
    for i=1:nTrial
        autoSim = autoCorr(simV{i});
        autoBi = autoCorr(biV{i});
        autoLi = autoCorr(liV{i});
        autoBi0 = autoCorr(biV0{i});
        autoLi0 = autoCorr(liV0{i});

        subplot(ceil(nTrial/2),2,i)
        hold on
        l0 = length(autoSim);
        l = length(autoBi);
        tau0 = linspace(-l0/2*tstep,l0/2*tstep,l0);
        tau = linspace(-l/2*tstep,l/2*tstep,l);
        autoSim = [autoSim(l0/2+1:l0);autoSim(1:l0/2)];
        autoBi = [autoBi(l/2+1:l);autoBi(1:l/2)];
        autoLi = [autoLi(l/2+1:l);autoLi(1:l/2)];
        autoBi0 = [autoBi0(l/2+1:l);autoBi0(1:l/2)];
        autoLi0 = [autoLi0(l/2+1:l);autoLi0(1:l/2)];
        plot(tau0, autoSim,'k');
        plot(tau, autoBi,'b');
        plot(tau, autoLi,'r');
        plot(tau, autoBi0,'c');
        plot(tau, autoLi0,'m');
        xlabel('\tau');
        ylabel('autocorr');
        title(['Trial ',num2str(i)]);
    end
    if ~isempty(ext)
        saveas(gcf,[p.theme,'-autoCorr.',ext],format);
    end
    fclose(datafid);
    if jNDfid
        fclose(jNDfid);
    end
    if plotInput
        fclose(tInfid);
    end
end
function auto = autoCorr(data)
    normedData = (data-mean(data))./std(data);
    L = length(normedData);
    l = 2^nextpow2(L);
    f = fft(data,l);
    psd = (abs(f)/l).^2;
    auto = ifft(psd);
end
