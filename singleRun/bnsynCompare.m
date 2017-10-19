function bsynCompare(theme,rE,rI,libfile,seed,run_t)
    fclose('all');
    picformat = '';
    type = 'RS_exc_Rat';
    ith = 1;
    sizeofsize = 'int64';
    cutoff = true;
    ignore_t = 0;
    plotDelta = true;
    plotInput = true;
    load(libfile,'dtRange','vRange','sPSP','n','tstep','gList');
    nt0 = size(sPSP,1);
    ndt = length(dtRange);
    nv = length(vRange);
    edur = (nt0-1)*tstep;
    ldur = dtRange(ndt) - ignore_t;
    tstep = edur/(nt0-1);
    run_nt = round(run_t/tstep) + 1;
    prefix = ['E',num2str(floor(rE)),'-I',num2str(floor(rI)),'-corrL',num2str(floor(ldur)),'-',num2str(run_t),'-',num2str(seed)];
    filename = [prefix,'jND-',theme,'.bin'];
    fid = fopen(filename);
    
    jndbSize = fread(fid,[1,1],sizeofsize);
    if jndbSize > 0
        jndbt = fread(fid,[jndbSize,1], 'double');
        jndbt = jndbt*tstep;
        jndbv = fread(fid,[jndbSize,1], 'double');
    end
    ncrossb = fread(fid,[1,1], sizeofsize);
    crossb = struct('t',cell(ncrossb,1), ...
                   'v',cell(ncrossb,1));
    for i=1:ncrossb
        crossbSize = fread(fid,[1,1],sizeofsize);
        crossb(i).t = fread(fid,[crossbSize,1], 'double')*tstep;
        crossb(i).v = fread(fid,[crossbSize,1], 'double');
    end
    
    jndlSize = fread(fid,[1,1],sizeofsize);
    if jndlSize > 0
        jndlt = fread(fid,[jndlSize,1], 'double');
        jndlt = jndlt*tstep;
        jndlv = fread(fid,[jndlSize,1], 'double');
    end
    ncrossl = fread(fid,[1,1], sizeofsize);
    crossl = struct('t',cell(ncrossl,1), ...
                   'v',cell(ncrossl,1));
    for i=1:ncrossl
        crosslSize = fread(fid,[1,1],sizeofsize);
        crossl(i).t = fread(fid,[crosslSize,1], 'double')*tstep;
        crossl(i).v = fread(fid,[crosslSize,1], 'double');
    end
    
    fclose(fid);
    FontSize = 16;
    set(0,'DefaultAxesFontSize',FontSize);
    set(0,'DefaultTextFontSize',FontSize-2);
    nE = 0;
    nI = 0;
    for i=1:n
        if gList(i) > 0
            nE = nE + 1;
        else
            nI = nI + 1;
        end
    end
    filename = [prefix,'Data-',theme,'.bin'];
    fid = fopen(filename);
    outputMat = fread(fid,[run_nt,4],'double');
    dendOut = fread(fid,[run_nt,nE+nI],'double');
    fclose(fid);
    t = 0:tstep:run_t;
    size(outputMat)
    simV = outputMat(:,1);
    biV  = outputMat(:,2);
    liV  = outputMat(:,3);
    biV0  = outputMat(:,4);
    filename = [prefix,'tIn-',theme,'.bin'];
    fid = fopen(filename);
    ntmp = fread(fid,[1,1],sizeofsize);
    tin = fread(fid,[ntmp,1],'double');
    tID = fread(fid,[ntmp,1],sizeofsize);
    fclose(fid);
    tID = tID + 1;
    pE = (tID <= nE);
    pI = (tID > nE);
    tE = tin(pE);
    tI = tin(pI);
    Eid = tID(pE);
    Iid = tID(pI)-nE;
    Ein = length(tE);
    Iin = length(tI);
    figure;
    textFontSize = 8;
    minV = min([min(simV),min(biV),min(liV),min(biV0)]);
    maxV = max([max(simV),max(biV),max(liV),max(biV0)]);
    subplot(2,2,1);
    %plot(t,[simV,simVq,biV,liV]);
    plot(t,[simV,biV,liV,biV0]);
    hold on
    if jndbSize
        plot(jndbt,jndbv,'*k');
    end
    if jndlSize
        plot(jndlt,jndlv,'sk');
    end
    if ncrossb
        for i=1:ncrossb
            plot(crossb(i).t,crossb(i).v,'-k');
        end
    end
    if ncrossl
        for i=1:ncrossl
            plot(crossl(i).t,crossl(i).v,':k');
        end
    end
    % plot(t,[simV,liV]);
    % plot(t,simV);
    vEtar = zeros(Ein,1);
    plot(tE,minV*ones(1,Ein),'.r');
    vItar = zeros(Iin,1);
    plot(tI,minV*ones(1,Iin),'.b');
    if plotInput
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
    legend({'sim','bi','li','bi0','jb','jl'});
    %legend({'sim','simq','bi','li'});
    ylim([minV,maxV]);
    subplot(2,2,2);
    plot(t,simV);
    hold on
    plot(t,dendOut);
    legend('show');
    ylim([minV,maxV]);
    subplot(2,2,3);
    hold on;
    signBiV = sign(sum(biV-simV>0)-sum(biV-simV<0));
    mB = mean(abs(biV-simV))*signBiV;
    errorbar(2,mB,std(abs(biV-simV)),'b');
    signLiV = sign(sum(liV-simV>0)-sum(liV-simV<0));
    mL = mean(abs(liV-simV))*signLiV;
    errorbar(1,mL,std(abs(liV-simV)),'r');
    signBiV0 = sign(sum(biV0-simV>0)-sum(biV0-simV<0));
    mB0 = mean(abs(biV0-simV))*signBiV0;
    errorbar(3,mB0,std(abs(biV0-simV)),'c');
    plot(2,mB,'*b');
    plot(1,mL,'*r');
    plot(3,mB0,'*c');
    xlim([0,4]);
    set(gca,'xtick',[1,2,3]);
    legend({'biV','liV','biV0'});
    plot([0,4],[0,0],':k');
    if plotDelta
        subplot(2,2,4);
        %plot(t,simV-simVq);
        hold on
        dliV = liV - simV;
        plot(t,dliV);
        dbiV = biV - simV;
        plot(t,dbiV);
        dbiV0 = biV0 - simV;
        plot(t,dbiV0);
        %legend({'simV-simVq','\Delta(sim-li)','\Delta(sim-bi)'});
        legend({'\Delta(sim-li)','\Delta(sim-bi)','\Delta(sim-bi0)'});
        minV = min(min(dbiV),min(dliV));
        maxV = max(max(dbiV),max(dliV));
        vEtar = zeros(Ein,1);
        plot(tE,minV*ones(1,Ein),'.r');
        targetV = dliV;
        for i=1:Ein
            iv = ceil(tE(i)/tstep);
            if iv==0, iv=iv+1,end
            jv = iv + 1;
            vEtar = targetV(iv) + mod(tE(i),tstep)/tstep * (targetV(jv)-targetV(iv));
            plot([tE(i),tE(i)],[minV,vEtar],':r');
            text(tE(i),minV+(vEtar-minV)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
            text(tE(i),minV+(vEtar-minV)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
         
            iv = ceil((tE(i)+edur)/tstep);
            if iv==0, iv=iv+1,end
            if iv+1 < run_nt
                jv = iv + 1;
                vEtar = targetV(iv) + mod(tE(i),tstep)/tstep * (targetV(jv)-targetV(iv));
                plot([tE(i),tE(i)]+edur,[vEtar,maxV],':r');
                text(tE(i)+edur,maxV-(maxV-vEtar)*0.3,num2str(Eid(i)),'Color','r','FontSize',textFontSize);
                text(tE(i)+edur,maxV-(maxV-vEtar)*0.1,num2str(i),'Color','r','FontSize',textFontSize);
            end
        end
        vItar = zeros(Iin,1);
        plot(tI,minV*ones(1,Iin),'.b');
        for i=1:Iin
            iv = ceil(tI(i)/tstep);
            if iv==0, iv=iv+1,end
            jv = iv + 1;
            vItar(i) = targetV(iv) + mod(tI(i),tstep)/tstep * (targetV(jv)-targetV(iv));
            plot([tI(i),tI(i)],[minV,vItar(i)],':b');
            text(tI(i),minV+(vItar(i)-minV)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
            text(tI(i),minV+(vItar(i)-minV)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
                
            iv = ceil((tI(i)+edur)/tstep);
            if iv==0, iv=iv+1,end
            if iv+1 < run_nt
                jv = iv + 1;
                vtar = targetV(iv) + mod(tI(i),tstep)/tstep * (targetV(jv)-targetV(iv));
                plot([tI(i),tI(i)]+edur,[vtar,maxV],':b');
                text(tI(i)+edur,maxV-(maxV-vtar)*0.3,num2str(Iid(i)),'Color','b','FontSize',textFontSize);
                text(tI(i)+edur,maxV-(maxV-vtar)*0.1,num2str(i),'Color','b','FontSize',textFontSize);
            end
        end
    end
    
    saveas(gcf,[theme,'-bnsyn.fig']);
end
