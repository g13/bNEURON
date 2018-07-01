function plotMultiGainCurve(inputFn, ext, plotSubthreshold, plotInput, plotDendV, sizeSize, plotAuto)
    textFontSize = 6;
    if nargin < 7
        plotAuto = true;
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
    end
    format = ext;
    if isequal(ext, 'psc')
        format = 'epsc';
    end
    if isequal(ext, 'jpg')
        format = 'jpeg';
    end
    p = read_cfg(inputFn);
    p.theme
    p.seed
    dimsFn = 'readoutDimension.bin';
    if isfield(p,'reformatInputFn')
        dimsFn = p.reformatInputFn;
    end
    if ~isstruct(p)
        return
    end
    fdr = struct2cell(dir());
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
    method = false(7,1);
    lcolor = {'k','b','r','.b','r','g','m'};
    dcolor = {'*k','*b','*r','ob','or','sb','sr'};
    mcolor = {'-*k','-*b','-*r','-ob','-or','-sb','-sr'};
    label = {'sim','bilinear','linear','jb','jl','bilinear0','linear0'};
    % name, fdr, date, bytes, isdir, datenum
    filter = [1,5];
    fdr = fdr(filter,:);
    fnstr = strjoin(fdr(1,:))
    dataFn = extractfn('Data',fnstr,p.theme,7);
    jNDFn = extractfn('jND',fnstr,p.theme,7);
    RasterFn = extractfn('Raster',fnstr,p.theme,7);
    cpuTimeFn = extractfn('cpuTime',fnstr,p.theme,7);
    filePattern = ['tIn-',p.theme,'\-s\d*\.bin'];
    match = regexp(fnstr, filePattern,'match');
    assert(length(match) == 1);
    tInFn = match{1}

    for i = 1:7
        m = regexp(strjoin(cpuTimeFn), ['-',num2str(i),'-',p.theme],'match')
        nm = length(m);
        if nm > 0
            assert(nm == 1);
            method(i) = true;
        end
    end
    nMethod = sum(method)
    in = find(method,1,'last')

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

    rasterSize = zeros(nTrial,7);
    rasterTime = cell(nTrial,7);
    figure;
    cm = 1;
    for im=1:7
        if method(im)
            RasterFid = fopen(RasterFn{im},'r');
            for i = 1:nTrial
                subplot(nTrial,2,(i-1)*2+1)
                hold on
                i
                im
                rasterSize(i,im) = fread(RasterFid,1,sizeSize);

                if rasterSize(i,im) > 0
                    rasterTime{i,im} = fread(RasterFid,rasterSize(i,im),'double');
                end
                plot(rasterTime{i,im},zeros(rasterSize(i,im),1)+cm,dcolor{im},'MarkerSize',2);
            end
            cm = cm + 1;
            fclose(RasterFid);
        end
    end
    save(['../Raster-',p.theme,'-s',num2str(p.seed),'-multi.mat'], 'rasterSize', 'rasterTime', 'method');
    subplot(2,2,2)
    hold on
    for im = 1:7
        if method(im)
            plot(inputLevel,rasterSize(:,im)./runTime*1000,mcolor{im});
        end
    end
    xlim([0,inputLevel(nTrial)*1.1]);
    xlabel('input rate Hz');
    ylabel('firing rate Hz');

    cpuTime = zeros(nTrial,7)-1;
    for im=1:7
        if method(im)
            cpuFid = fopen(cpuTimeFn{im},'r');
            cpuTime(:,im) = fread(cpuFid,[nTrial,1],'double');
        end
    end
    subplot(2,2,4)
    hold on
    line = [];
    for im = 1:7
        if method(im)
            plot(inputLevel,cpuTime(:,im),mcolor{im});
        end
    end
    xlim([0,inputLevel(nTrial)*1.1]);
    legend(label(method));
    xlabel('input rate Hz');
    ylabel('cpuTime s');

    if ~isempty(ext)
        saveas(gcf,[p.theme,'-gainCurve.',ext],format);
    end
    if plotSubthreshold || plotAuto
        dendV = cell(nTrial,n);
        v = cell(nTrial,7);
        if plotInput
            tInfid = fopen(tInFn,'r');
        end
        disp(tstep);
        for im = [1,2,3,6,7]
            datafid = fopen(dataFn{im},'r');
            for i=1:nTrial
                if im==1
                    v{i,1} = fread(datafid, nDimSim(i), 'double');
                    if p.getDendV
                        for j=1:n
                            dendV{i,j} = fread(datafid, nDimSim(i), 'double');
                        end
                    end
                else
                    v{i,im} = fread(datafid, nDim(i), 'double');
                end
            end
            fclose(datafid);
        end
    end
    if plotSubthreshold
        if method(4) 
            jNDfid = fopen(jNDFn{4},'r');
            jbt = cell(nTrial,1);
            jbCrossT = cell(nTrial,1);
            jbCrossV = cell(nTrial,1);
            jbnCross = zeros(nTrial,1);
            for i=1:nTrial
                jbSize = fread(jNDfid,1,sizeSize);
                jbt{i} = fread(jNDfid,jbSize,'double')*tstep;
                v{i,4} = fread(jNDfid,jbSize,'double');
                jbnCross(i) = fread(jNDfid,1,'int')
                jbCrossT{i} = cell(jbnCross(i),1);
                jbCrossV{i} = cell(jbnCross(i),1);
                for j=1:jbnCross(i)
                    tmpSize = fread(jNDfid, 1, sizeSize)
                    jbCrossT{i}{j} = fread(jNDfid, tmpSize, 'double')*tstep
                    jbCrossV{i}{j} = fread(jNDfid, tmpSize, 'double')
                end
            end
            fclose(jNDfid);
        end
        if method(5) 
            jNDfid = fopen(jNDFn{5},'r');
            jlt = cell(nTrial,1);
            jlCrossT = cell(nTrial,1);
            jlCrossV = cell(nTrial,1);
            jlnCross = zeros(nTrial,1);
            for i=1:nTrial
                jlSize = fread(jNDfid,1,sizeSize);
                jlt{i} = fread(jNDfid,jlSize,'double')*tstep;
                v{i,5} = fread(jNDfid,jlSize,'double');
                jlnCross(i) = fread(jNDfid,1,'int')
                jlCrossT{i} = cell(jlnCross(i),1);
                jlCrossV{i} = cell(jlnCross(i),1);
                for j=1:jlnCross(i)
                    tmpSize = fread(jNDfid, 1, sizeSize)
                    jlCrossT{i}{j} = fread(jNDfid, tmpSize, 'double')*tstep;
                    jlCrossV{i}{j} = fread(jNDfid, tmpSize, 'double');
                end
            end
            fclose(jNDfid);
        end
        for i=1:nTrial
            run_nt = round(runTime(i)/tstep) + 1;
            disp(['this is ', num2str(i), 'th trial']);
            disp(dt(i));
            t0 = linspace(0,nDimSim(i)-1,nDimSim(i)) * dt(i);
            t = linspace(0,nDim(i)-1,nDim(i)) * tstep;
            assert(t(end)==t0(end));
            xl = [0, t(end)];

            handles = [];
            figure;
            subplot(2,1,1)
            hold on
            minV0 = 1000;
            maxV0 = -1000;
            if method(1)
                handle = plot(t0,v{i,1},lcolor{1});
                handles = [handles,handle];
                [minV0, maxV0] = minmax(v{i,1},minV0,maxV0);
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
            end
            for im = [2,3,6,7]
                if method(im)
                    handle = plot(t,v{i,im},lcolor{im});
                    handles = [handles, handle];
                    [minV0, maxV0] = minmax(v{i,im},minV0,maxV0);
                end
            end
            vStretch = maxV0-minV0;
            minV = minV0-vStretch*0.1;
            maxV = maxV0+vStretch*0.1;
            yl = [minV,maxV];
            ylim(yl);
            if method(4) 
                handle = plot(jbt{i},v{i,4},'.b','MarkerSize',3);
                for j=1:jbnCross(i)
                    plot(jbCrossT{i}{j}, jbCrossV{i}{j},':b');
                end
                handles = [handles, handle];
            end
            if method(5) 
                handle = plot(jlt{i},v{i,5},'.r','MarkerSize',3);
                for j=1:jlnCross(i)
                    plot(jlCrossT{i}{j}, jlCrossV{i}{j},':r');
                end
                handles = [handles, handle];
            end
            if plotDendV && p.getDendV && method(1)
                handles = [handles, dline];
                labels = [label(method), strcat('dend[',cellstr(num2str(loc))',']')];
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
                    vEtar = v{i,in}(iv) + mod(tE(j),tstep)/tstep * (v{i,in}(jv)-v{i,in}(iv));
                    plot([tE(j),tE(j)],[minV,vEtar],':r');
                    text(tE(j),minV+(vEtar-minV)*0.3,num2str(Eid(j)),'Color','r','FontSize',textFontSize);
                    text(tE(j),minV+(vEtar-minV)*0.1,num2str(j),'Color','r','FontSize',textFontSize);
                 
                    iv = ceil((tE(j)+ldur)/tstep);
                    if iv==0, iv=iv+1,end
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vEtar = v{i,in}(iv) + mod(tE(j),tstep)/tstep * (v{i,in}(jv)-v{i,in}(iv));
                        plot([tE(j),tE(j)]+ldur,[vEtar,maxV],':r');
                        text(tE(j)+ldur,maxV-(maxV-vEtar)*0.3,num2str(Eid(j)),'Color','r','FontSize',textFontSize);
                        text(tE(j)+ldur,maxV-(maxV-vEtar)*0.1,num2str(j),'Color','r','FontSize',textFontSize);
                    end
                end
                for j=1:Iin
                    iv = ceil(tI(j)/tstep);
                    if iv==0, iv=iv+1,end
                    jv = iv + 1;
                    vItar(j) = v{i,in}(iv) + mod(tI(j),tstep)/tstep * (v{i,in}(jv)-v{i,in}(iv));
                    plot([tI(j),tI(j)],[minV,vItar(j)],':b');
                    text(tI(j),minV+(vItar(j)-minV)*0.3,num2str(Iid(j)),'Color','b','FontSize',textFontSize);
                    text(tI(j),minV+(vItar(j)-minV)*0.1,num2str(j),'Color','b','FontSize',textFontSize);
                        
                    iv = ceil((tI(j)+ldur)/tstep);
                    if iv==0, iv=iv+1,end
                    if iv+1 < run_nt
                        jv = iv + 1;
                        vtar = v{i,in}(iv) + mod(tI(j),tstep)/tstep * (v{i,in}(jv)-v{i,in}(iv));
                        plot([tI(j),tI(j)]+ldur,[vtar,maxV],':b');
                        text(tI(j)+ldur,maxV-(maxV-vtar)*0.3,num2str(Iid(j)),'Color','b','FontSize',textFontSize);
                        text(tI(j)+ldur,maxV-(maxV-vtar)*0.1,num2str(j),'Color','b','FontSize',textFontSize);
                    end
                end
            end
            xlim(xl);

            if dt(i) == tstep && method(1)
                subplot(2,2,3)
            else 
                subplot(2,1,2)
            end
            hold on
            cm = 1;
            for im = 1:7
                if method(im)
                    plot(rasterTime{i,im},zeros(1,rasterSize(i,im))+cm,dcolor{im});
                    cm = cm + 1;
                end
            end
            legend(label(method));
            ylim([0,cm]);
            xlim(xl);

            if dt(i) == tstep && method(1)
                subplot(2,2,4)
                hold on
                cm = 1;
                for im = [2,3,6,7]
                    err = v{i,im} - v{i,1};
                    plot(t0, err, lcolor{im});
                    errorbar(xl(2)+cm, mean(abs(err)), std(err),dcolor{im});
                    cm = cm + 1;
                end
                xlim([xl(1),xl(2)+cm]);
            end

            if ~isempty(ext)
                saveas(gcf,[p.theme,'-trial',num2str(i),'.',ext],format);
            end
        end
    end
    if plotAuto
        figure;
        for i=1:nTrial
            %subplot(ceil(nTrial/2),2,i)
            %hold on
            for im=[1,2,3,6,7]
                if method(im)
                    [auto,psd,normed] = autoCorr(v{i,im});
                    l = length(auto);
                    ns = 2/tstep+1;
                    auto = [auto(l/2+1:l);auto(1:l/2)];
                    tau = linspace(-l/2*tstep,l/2*tstep,l);
                    subplot(nTrial,3,(i-1)*3+1)
                    plot(linspace(0,(length(normed)-1)*tstep,length(normed)), normed, lcolor{im});
                    hold on
                    title(['Trial ',num2str(i)]);
                    subplot(nTrial,3,(i-1)*3+2)
                    semilogy(linspace(0,(l/2-1),l/2).*1000./tstep./l,smooth(psd(1:(l/2)),ns),lcolor{im});
                    hold on
                    xlabel('Hz')
                    subplot(nTrial,3,(i-1)*3+3)
                    plot(tau, auto,lcolor{im});
                    hold on
                end
            end
            %xlabel('\tau');
            %ylabel('autocorr');
            %title(['Trial ',num2str(i)]);
        end
        if ~isempty(ext)
            saveas(gcf,[p.theme,'-autoCorr.',ext],format);
        end
        if plotInput
            fclose(tInfid);
        end
    end
end
function [auto,psd,normedData] = autoCorr(data)
    normedData = (data-mean(data))./std(data);
    L = length(normedData);
    l = 2^nextpow2(L);
    f = fft(data,l);
    psd = (abs(f)/l).^2;
    auto = ifft(psd);
end
function files = extractfn(type,fnstr,theme,n)
    filePattern = [type,'-\d-',theme,'\-s\d*\.bin'];
    files0 = regexp(fnstr, filePattern,'match')
    files = cell(n,1);
    l = length(type);
    for i = 1:length(files0) 
        j = str2num(files0{i}(l+2));
        files{j} = files0{i};
    end
end
function [minV0, maxV0] = minmax(v,minV0,maxV0);
    nv = min(v);
    if nv < minV0
        minV0 = nv;
    end
    mv = max(v);
    if mv > maxV0
        maxV0 = mv;
    end
end
