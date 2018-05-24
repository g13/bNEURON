function exploreK(k,v,tRange,dt,dti,nv,ndt,nloc1,nloc2,run_nt,tstep,label,distanceMat)
    n0 = nloc1*nloc2;
    load(distanceMat);
    dist = zeros(n0,1);
    cc = 0;
    for i=1:nloc2
        for j=1:nloc1
            cc = cc + 1;
            dist(cc) = abs(distance(i)-distance(j));
            if distance(i) > distance(j)
                dist(cc) = dist(cc) + distance(i);
            else
                dist(cc) = dist(cc) + distance(j);
            end
        end
    end
    [~,did] = sort(dist);
    nt = length(tRange);
    sat = 0.8;
    val = 1.0;
    %%
    t0=(0:run_nt)*tstep;
    idt = 1;
    iv = 2;
    sel = dti(idt) + (100:300);
    t = t0(sel);
    tar = reshape(k(sel,:,:,idt,iv),length(sel),nloc1*nloc2);
    figure;
    subplot(2,2,1)
    hold on
    for i=1:nloc1*nloc2
        hue = (did(i)-1)/(n0-1) * (2/3);
        plot(t,tar(:,i),'Color',hsv2rgb([hue,sat,val]));
    end
    xlim([t(1)-5,t(end)+5]);
    xlabel('t ms');
    ylabel(['k',label]);
    title('dist:color');
    
    subplot(2,2,2)
    plot(mean(tar),std(tar),'*r');
    xlabel(['k',label,'<t>']);
    ylabel(['\sigma(k',label,'<t>)']);
    subplot(2,2,3)
    plot(std(tar)./abs(mean(tar)));
    xlabel('loc ID');
    ylabel('CV of t');
    subplot(2,2,4)
    hold on
    sat = 1.0;
    val = 1.0;
    for iv = nv:-1:1
        mks = iv*2;
        for idt = 1:ndt
            sel = dti(idt) + (100:300);
            hue = (idt-1)/(ndt-1) * (2/3);
            tar = reshape(k(sel,:,:,idt,iv),length(sel),nloc1*nloc2);
            plot(std(tar)./abs(mean(tar)),'o','MarkerSize',mks,'Color',hsv2rgb([hue,sat,val]));
        end
    end
    ylabel('CV of t');
    xlabel('pair #');
    title('v:mks, dt:color(r2b)')
    %%
    idt = 1;
    sel = dti(idt) + 200;
    tar = reshape(k(sel,:,:,idt,:),nloc1*nloc2,nv)';
    figure;
    subplot(2,2,1)
    hold on
    for i=1:nloc1*nloc2
        hue = (did(i)-1)/(n0-1) * (2/3);
        plot(v,tar(:,i),'Color',hsv2rgb([hue,sat,val]));
    end
    xlabel('v mV');
    ylabel(['k',label]);
    title('dist:color');
    subplot(2,2,2)
    plot(mean(tar),std(tar),'*r');
    xlabel(['k',label,'<v>']);
    ylabel(['\sigma(k',label,'<v>)']);
    subplot(2,2,3)
    plot(std(tar)./abs(mean(tar)));
    xlabel('loc ID');
    ylabel('CV of v');
    subplot(2,2,4)
    hold on
    sat = 1.0;
    val = 1.0;
    for it = 1:nt
        mks = 2*it;
        for idt = 1:ndt
            sel = dti(idt) + tRange(it);
            tar = reshape(k(sel,:,:,idt,:),nloc1*nloc2,nv)';
            hue = (idt-1)/(ndt-1) * (2/3);
            plot(std(tar)./abs(mean(tar)),'o','MarkerSize',mks,'Color',hsv2rgb([hue,sat,val]));
        end
    end
    ylabel('CV of v');
    xlabel('pair #');
    title('t:mks, dt:color(r2b)')
    %%
    sel = 200;
    iv = 2;
    tar0 = zeros(nloc1,nloc2,ndt,nv);
    for idt=1:ndt
        tar0(:,:,idt,iv) = k(dti(idt)+sel,:,:,idt,iv);
    end
    tar = reshape(tar0(:,:,:,iv),nloc1*nloc2,ndt)';
    figure;
    subplot(2,2,1)
    hold on
    for i=1:nloc1*nloc2
        hue = (did(i)-1)/(n0-1) * (2/3);
        plot(dt,tar(:,i),'Color',hsv2rgb([hue,sat,val]));
    end
    xlabel('dt ms');
    ylabel(['k',label]);
    title('dist:color');
    subplot(2,2,2)
    plot(mean(tar),std(tar),'*r');
    xlabel(['k',label,'<dt>']);
    ylabel(['\sigma(k',label,'<dt>)']);
    subplot(2,2,3)
    plot(std(tar)./abs(mean(tar)));
    xlabel('loc ID');
    ylabel('CV of dt');
    subplot(2,2,4)
    hold on
    sat = 1.0;
    val = 1.0;
    for iv = nv:-1:1
        mks = iv*2;
        for it = 1:nt
            for idt=1:ndt
                sel = dti(idt) + tRange(it);
                tar0(:,:,idt,iv) = k(sel,:,:,idt,iv);
            end
            hue = (it-1)/(nt-1) * (2/3);
            tar = reshape(tar0(:,:,:,iv),nloc1*nloc2,ndt)';
            plot(std(tar)./abs(mean(tar)),'o','MarkerSize',mks,'Color',hsv2rgb([hue,sat,val]));
        end
    end
    ylabel('CV of dt');
    xlabel('pair #');
    title('v:mks, t:color(r2b)') 
end