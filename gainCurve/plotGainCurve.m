function plotGainCurve(inputFn, ext, plotSubthreshold, plotInput, sizeSize)
    if nargin < 5
        sizeSize = 'int64';
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
    if isequal(ext, 'psc')
        format = 'epsc';
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
    load(p.libFile,'tstep','n','gList');

    nE = 0;
    nI = 0;
    for i=1:n
        if gList(i) > 0
            nE = nE + 1;
        else
            nI = nI + 1;
        end
    end
    datafilePattern = [p.theme,'\-s\d*\-(Data|jND|Raster|tIn)\.bin'];
    datafile = regexp(fnstr, datafilePattern,'match');
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
    for i = 1:nTrial
        ss(i) = fread(RasterFid,1,sizeSize);
        if ss(i) >0
            ts{i} = fread(RasterFid,ss(i),'double');
        end
        sb(i) = fread(RasterFid,1,sizeSize);
        if sb(i) >0
            tb{i} = fread(RasterFid,sb(i),'double');
        end
        sl(i) = fread(RasterFid,1,sizeSize);
        if sl(i) >0
            tl{i} = fread(RasterFid,sl(i),'double');
        end
        sjb(i) = fread(RasterFid,1,sizeSize);
        if sjb(i) >0
            tjb{i} = fread(RasterFid,sjb(i),'double');
        end
        sjl(i) = fread(RasterFid,1,sizeSize);
        if sjl(i) >0
            tjl{i} = fread(RasterFid,sjl(i),'double');
        end
        sb0(i) = fread(RasterFid,1,sizeSize);
        if sb0(i) >0
            tb0{i} = fread(RasterFid,sb0(i),'double');
        end
    end
    fclose(RasterFid);
    figure;
    hold on
    plot(inputLevel,ss./runTime,'-*k');
    plot(inputLevel,sb./runTime,'-*r');
    plot(inputLevel,sl./runTime,'-*b');
    plot(inputLevel,sjb./runTime,'-or');
    plot(inputLevel,sjl./runTime,'-ob');
    plot(inputLevel,sb0./runTime,'-og');
    xlim([0,inputLevel(nTrial)*1.1]);

    if ~isempty(ext)
        saveas(gcf,[p.theme,'-gainCurve.',ext],format);
    end
    if plotSubthreshold
        dataFn = datafile{1};
        jNDFn = datafile{3};
        tInFn = datafile{4};
        jbt = cell(nTrial,1);
        jbv = cell(nTrial,1);
        jbCrossT = cell(nTrial,1);
        jbCrossV = cell(nTrial,1);
        jlt = cell(nTrial,1);
        jlv = cell(nTrial,1);
        jlCrossT = cell(nTrial,1);
        jlCrossV = cell(nTrial,1);
        dendV = cell(nTrial,n);
        simV = cell(nTrial,1);
        biV = cell(nTrial,1);
        biV0 = cell(nTrial,1);
        datafid = fopen(dataFn,'r');
        jNDfid = fopen(jNDFn,'r');
        if plotInput
            tInfid = fopen(tInFn,'r');
        end
        disp(tstep);
        for i=1:nTrial
            disp(['this is ', num2str(i), 'th trial']);
            disp(dt(i));
            t0 = linspace(0,nDimSim(i)-1,nDimSim(i)) * dt(i);
            t = linspace(0,nDim(i)-1,nDim(i)) * tstep;
            assert(t(end)==t0(end));
            figure;
            subplot(2,1,1)
            hold on
            if ss(i) > 0
                plot(ts{i},ones(1,ss(i)),'.k','MarkerSize',6,'LineStyle','none');
            end
            if sb(i) > 0
                plot(tb{i},ones(1,sb(i))+1,'.b','MarkerSize',6,'LineStyle','none');
            end
            if sl(i) > 0
                plot(tl{i},ones(1,sl(i))+2,'.r','MarkerSize',6,'LineStyle','none');
            end
            if sjb(i) > 0
                plot(tjb{i},ones(1,sjb(i))+3,'ob','MarkerSize',3,'LineStyle','none');
            end
            if sjl(i) > 0
                plot(tjl{i},ones(1,sjl(i))+4,'or','MarkerSize',3,'LineStyle','none');
            end
            if sb0(i) > 0
                plot(tb0{i},ones(1,sb0(i))+5,'og','MarkerSize',3,'LineStyle','none');
            end
            ylim([0,6]);
            xl = [0, t(end)];
            xlim(xl);

            subplot(2,1,2)
            hold on
            if datafid
                simV{i} = fread(datafid, nDimSim(i), 'double');
                biV{i} = fread(datafid, nDim(i), 'double');
                liV{i} = fread(datafid, nDim(i), 'double');
                biV0{i} = fread(datafid, nDim(i), 'double');
                for j=1:n
                    dendV{i,j} = fread(datafid, nDimSim(i), 'double');
                end
            end
            plot(t0,simV{i},'k');
            plot(t,biV{i},'b');
            plot(t,liV{i},'r');
            plot(t,biV0{i},'g');
            minV = min([min(simV{i}),min(biV{i}),min(liV{i}),min(biV0{i})]);
            maxV = min([-50,max([max(simV{i}),max(biV{i}),max(liV{i}),max(biV0{i})]]);
            vStretch = maxV-minV;
            yl = [minV-vStretch*0.1,maxV+vStretch*0.1];
            ylim(yl);
            xlim(xl);
            if jNDfid
                jbSize = fread(jNDfid,1,sizeSize);
                jbt{i} = fread(jNDfid,jbSize,'double')*tstep;
                jbv{i} = fread(jNDfid,jbSize,'double');
                plot(jbt{i},jbv{i},'.b','MarkerSize',3);
                jbnCross = fread(jNDfid,1,'int');
                jbCrossT{i} = [];
                for j=1:jbnCross
                    tmpSize = fread(jNDfid, 1, sizeSize);
                    tmpCrossT = fread(jNDfid, tmpSize, 'double')*tstep;
                    tmpCrossV = fread(jNDfid, tmpSize, 'double');
                    jbCrossT{i} = [jbCrossT{i}, tmpCrossT];
                    jbCrossV{i} = [jbCrossV{i}, tmpCrossV];
                    plot(tmpCrossT, tmpCrossV,':b');
                end
                jlSize = fread(jNDfid,1,sizeSize);
                jlt{i} = fread(jNDfid,jlSize,'double')*tstep;
                jlv{i} = fread(jNDfid,jlSize,'double');
                plot(jlt{i},jlv{i},'.r','MarkerSize',3);
                jlnCross = fread(jNDfid,1,'int');
                jlCrossT{i} = [];
                for j=1:jbnCross
                    tmpSize = fread(jNDfid, 1, sizeSize);
                    tmpCrossT = fread(jNDfid, tmpSize, 'double')*tstep;
                    tmpCrossV = fread(jNDfid, tmpSize, 'double');
                    jlCrossT{i} = [jlCrossT{i}, tmpCrossT];
                    jlCrossV{i} = [jlCrossV{i}, tmpCrossV];
                    plot(tmpCrossT, tmpCrossV,':r');
                end
            end
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
                for i=1:Ein
                    iv = ceil(tE(i)/tstep);
                    if iv==0, iv=iv+1,end
                    jv = iv + 1;
                    vEtar = simV(iv) + mod(tE(i),tstep)/tstep * (simV(jv)-simV(iv));
                    plot([tE(i),tE(i)],[minV,vEtar],':r');
                    text(tE(i),minV+(vEtar-minV)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                    text(tE(i),minV+(vEtar-minV)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
                 
                    iv = ceil((tE(i)+ldur)/tstep);
                    if iv==0, iv=iv+1,end
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vEtar = simV(iv) + mod(tE(i),tstep)/tstep * (simV(jv)-simV(iv));
                        plot([tE(i),tE(i)]+ldur,[vEtar,maxV],':r');
                        text(tE(i)+ldur,maxV-(maxV-vEtar)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                        text(tE(i)+ldur,maxV-(maxV-vEtar)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
                    end
                end
                for i=1:Iin
                    iv = ceil(tI(i)/tstep);
                    if iv==0, iv=iv+1,end
                    jv = iv + 1;
                    vItar(i) = simV(iv) + mod(tI(i),tstep)/tstep * (simV(jv)-simV(iv));
                    plot([tI(i),tI(i)],[minV,vItar(i)],':b');
                    text(tI(i),minV+(vItar(i)-minV)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                    text(tI(i),minV+(vItar(i)-minV)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
                        
                    iv = ceil((tI(i)+ldur)/tstep);
                    if iv==0, iv=iv+1,end
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vtar = simV(iv) + mod(tI(i),tstep)/tstep * (simV(jv)-simV(iv));
                        plot([tI(i),tI(i)]+ldur,[vtar,maxV],':b');
                        text(tI(i)+ldur,maxV-(maxV-vtar)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                        text(tI(i)+ldur,maxV-(maxV-vtar)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
                    end
                end
            end
            if ~isempty(ext)
                saveas(gcf,[p.theme,'-trial',num2str(i),'.',ext],format);
            end
        end
    end
    fclose(datafid);
    fclose(jNDfid);
    fclose(tInfid);
end
