function readNib4(dir)
    theme = dir
    file = [dir,'/',theme,'-vRange.bin'];
    fid = fopen(file);
    vRange = fread(fid,[1,inf],'double');
    nv0 = length(vRange);
    fclose(fid);
    file = [dir,'/',theme,'-dtRange.bin'];
    fid = fopen(file);
    dtRange = fread(fid,[1,inf],'double');
    ndt = length(dtRange);
    fclose(fid);
    file = [dir,'/',theme,'-p.bin'];
    fid = fopen(file);
    data = fread(fid,[1,inf],'double');
    tstep = data(1);
    dur = data(2);
    trans = data(3); 
    n = int32(data(4));
    disp([num2str(n),' Synapses']);
    fclose(fid);

    file = [dir,'/',theme,'-R.bin'];
    fid = fopen(file);
    loc = fread(fid,[n,1],'int64')
    pos = fread(fid,[n,1],'double')
    gList = fread(fid,[n,1],'double')
    fclose(fid);

    nt = round(dur/tstep)+1;
    
    kv = zeros(nt,ndt,n,n,ndt,nv0);
    tmax = zeros(n,ndt,nv0);
    sfire = zeros(n,ndt,nv0);
    sPSP = zeros(nt,n,ndt,nv0);
    dendv = zeros(nt,n,ndt,nv0);
    vleakage = zeros(nt,nv0);
    dendvleak = zeros(nt,n,nv0);
    bfire = zeros(ndt,n,n,ndt,nv0);
    for i=1:nv0
        disp(i);
        file = [dir,'/',theme,'-V',num2str(i-1),'.bin']
        fid = fopen(file);
        tmp  = fread(fid,[nt,1],'double'); 
        vleakage(:,i) = tmp;
        tmp  = fread(fid,[nt,n],'double'); 
        dendvleak(:,:,i) = tmp;
        for idt = 1:ndt
            tmp = fread(fid,[nt,n],'double'); 
            sPSP(:,:,idt,i) = tmp;
        end
        for idt = 1:ndt
            tmp = fread(fid,[n,1],'int64');
            tmax(:,idt,i) = tmp;
        end
        for idt = 1:ndt
            tmp = fread(fid,[n,1],'int64');
            sfire(:,idt,i) = tmp;
        end
        for idt = 1:ndt
            tmp = fread(fid,[nt,n],'double');
            dendv(:,:,idt,i) = tmp;
        end
        for idt = 1:ndt
            tmp = fread(fid,[nt,n*n*ndt],'double'); 
            kv(:,:,:,:,idt,i) = reshape(tmp,[nt,ndt,n,n]);
            max(max(abs(tmp)))
        end
        for idt = 1:ndt
            tmp = fread(fid,[n,n*ndt],'double'); 
            bfire(:,:,:,idt,i) = reshape(tmp,[ndt,n,n]);
        end
        fclose(fid);
    end
    kv0 = zeros(nt,n,n,ndt);
    file = [dir,'/',theme,'-kv0.bin']
    fid = fopen(file);
    for idt = 1:ndt
        tmp = fread(fid,[nt,n*n],'double');
        kv0(:,:,:,idt) = reshape(tmp,[nt,n,n]);
    end
    fclose(fid);
    sf = sum(sfire==0,3) 
    bf = sum(bfire==0,5);
    disp('sf:')
    size(sf)
    disp('bf:')
    size(bf)
    fireCap = zeros(size(sf'));
    for i=1:n
        size(reshape(bf(:,i,:,:),[ndt,n*ndt]));
        size(min([sf(i,:)',reshape(bf(:,i,:,:),[ndt,n*ndt])],2));
        fireCap(:,i) = min([sf(i,:)',reshape(bf(:,i,:,:),[ndt,n*ndt])],[],2);
    end
    disp('fireCap:')
    size(fireCap)
    %t = 0:tstep:dur;
    %figure
    %hold on
    %subplot(1,2,1);
    %plot(t,vleakage(:,1));
    %subplot(1,2,2);
    %hold on
    %plot(t,sPSP(:,1,1,1));
    %plot(t,kv(:,2,1,2,2,1));
    save(['../data/',theme,'-NEURON'],'kv0','kv','sPSP','tmax','vleakage','dendvleak','sf','bf','fireCap','dendv','loc','pos','n','gList','dtRange','vRange','tstep','-v7.3');
end
