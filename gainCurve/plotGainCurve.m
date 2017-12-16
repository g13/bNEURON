function plotGainCurve(inputFn, ext, plotSubthreshold, sizeSize)
    if nargin < 4
        sizeSize = 'int64';
        if nargin < 3
            plotSubthreshold = false;
        end
    end
    if ext == psc
        format = 'epsc'
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
    fnstr = strjoin(fdr(1,:));
    load(p.libFile,'tstep','n');
    datafilePattern = [p.theme,'\-s\d*\-(Data|jND|Raster|tIn)\.bin'];
    datafile = regexp(fnstr, datafilePattern,'match','forceCellOutput');
    dimsFid = fopen(dimsFn,'r');
    if dimsFid
        nTrial = fread(dimsFid, 1, 'int');
        if plotSubthreshold
            nDimSim = fread(dimsFid, nTrial, 'unsigned long');
            nDim = fread(dimsFid, nTrial, 'unsigned long');
            dt = fread(dimsFid, nTrial, 'double');
            inputLevel = fread(dimsFid, nTrial, 'double');
            runTime = fread(dimsFid, nTrial, 'double');
        end
        fclose(dimsFid);
    end
    RasterFn = datafile{1}{2};
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
        ts{i} = fread(RasterFid,ss,'double');
        sb(i) = fread(RasterFid,1,sizeSize);
        tb{i} = fread(RasterFid,sb,'double');
        sl(i) = fread(RasterFid,1,sizeSize);
        tl{i} = fread(RasterFid,sl,'double');
        sjb(i) = fread(RasterFid,1,sizeSize);
        tjb{i} = fread(RasterFid,sjb,'double');
        sjl(i) = fread(RasterFid,1,sizeSize);
        tjl{i} = fread(RasterFid,sjl,'double');
        sb0(i) = fread(RasterFid,1,sizeSize);
        tb0{i} = fread(RasterFid,sb0,'double');
    end
    fclose(RasterFid);

    figure;
    hold on
    plot(inputLevel,ss./runTime,'*k');
    plot(inputLevel,sb./runTime,'*r');
    plot(inputLevel,sl./runTime,'*b');
    plot(inputLevel,sjb./runTime,'or');
    plot(inputLevel,sjl./runTime,'ob');
    plot(inputLevel,sb0./runTime,'og');

    saveas(gcf,[p.theme,'-gainCurve.',ext],format);
    if plotSubthreshold
        dataFn = datafile{1}{1};
        jNDFn = datafile{1}{3};
        tInFn = datafile{1}{4};
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
        for i=1:nTrial
            t0 = linspace(1,nDimSim(i)) * tstep;
            t = linspace(1,nDim(i)) * dt(i);
            figure;
            subplot(2,1,1)
            hold on
            plot(ts,ones(1,ss),'.k','MarkerSize',6,'LineStyle','none');
            plot(tb,ones(1,sb)+1,'.b','MarkerSize',6,'LineStyle','none');
            plot(tl,ones(1,sl)+2,'.r','MarkerSize',6,'LineStyle','none');
            plot(tjb,ones(1,sjb)+3,'ob','MarkerSize',3,'LineStyle','none');
            plot(tjl,ones(1,sjl)+4,'or','MarkerSize',3,'LineStyle','none');
            plot(tb0,ones(1,sb0)+5,'og','MarkerSize',3,'LineStyle','none');
            ylim([0,6]);
            xlim([0,t[nDim(i)]);

            subplot(2,1,2)
            hold on
            if datafid
                for i=1:nTrial
                    simV{i} = fread(datafid, nDimSim(i), 'double');
                    biV{i} = fread(datafid, nDim(i), 'double');
                    liV{i} = fread(datafid, nDim(i), 'double');
                    biV0{i} = fread(datafid, nDim(i), 'double');
                    for j=1:n
                        dendV{i,j} = fread(datafid, nDimSim(i), 'double');
                    end
                end
            end
            plot(t0,simV{i},'k');
            plot(t,biV{i},'r');
            plot(t,liV{i},'b');
            plot(t,biV0{i},'g');
            if jNDfid
                jbSize = fread(jNDfid,1,sizeSize);
                jbt{i} = fread(jNDfid,jbSize,'double');
                jbv{i} = fread(jNDfid,jbSize,'double');
                plot(jbt{i},jbv{i},'*r','LineStyle','none');
                jbnCross = fread(jNDfid,1,'int');
                for j=1:jbnCross
                    tmpSize = fread(jNDfid, 1, sizeSize);
                    tmpCrossT = fread(jNDfid, tmpSize, 'double');
                    tmpCrossV = fread(jNDfid, tmpSize, 'double');
                    jbCrossT{i} = [jbCrossT{i}, tmpCrossT];
                    jbCrossV{i} = [jbCrossV{i}, tmpCrossV];
                    plot(tmpCrossT, tmpCrossV,'b');
                end
                jlSize = fread(jNDfid,1,sizeSize);
                jlt{i} = fread(jNDfid,jlSize,'double');
                jlv{i} = fread(jNDfid,jlSize,'double');
                jlnCross = fread(jNDfid,1,'int');
                plot(jlt{i},jlv{i},'*b','LineStyle','none');
                for j=1:jbnCross
                    tmpSize = fread(jNDfid, 1, sizeSize);
                    tmpCrossT = fread(jNDfid, tmpSize, 'double');
                    tmpCrossV = fread(jNDfid, tmpSize, 'double');
                    jlCrossT{i} = [jlCrossT{i}, tmpCrossT];
                    jlCrossV{i} = [jlCrossV{i}, tmpCrossV];
                    plot(tmpCrossT, tmpCrossV,'r');
                end
            end
            saveas(gcf,[p.theme,'-trial',num2str(i),'.',ext],format);
        end
    end
    fclose(datafid);
    fclose(jNDfid);
end
