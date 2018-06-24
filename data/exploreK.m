function exploreK(k,Rs,v,tRange,dt,dti,nv,ndt,nloc1,nloc2,run_nt,tstep,label,distanceMat,ext,directory)
    format = ext;
    if isequal(ext, 'psc')
        format = 'epsc';
    end
    if isequal(ext, 'jpg')
        format = 'jpeg';
    end
    w = 1;
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
    %% t
    t0=(0:run_nt)*tstep;
    idt = 1;
    iv = 2;
    sel = dti(idt) + (100:300);
    t = t0(sel);
    tar = reshape(k(sel,:,:,idt,iv),length(sel),nloc1*nloc2);
    cvTar = squeeze(std(k(sel,:,:,idt,iv),w)./abs(mean(k(sel,:,:,idt,iv))));
    RsTar = reshape(Rs(sel,:,:,idt,iv),length(sel),nloc1*nloc2);
    cvRsTar = squeeze(std(Rs(sel,:,:,idt,iv),w)./abs(mean(Rs(sel,:,:,idt,iv))));
    figure;
    subplot(3,2,1)
        hold on
        for i=1:nloc1*nloc2
            hue = (did(i)-1)/(n0-1) * (2/3);
            plot(t,tar(:,i),'Color',hsv2rgb([hue,sat,val]));
        end
        xlim([t(1)-5,t(end)+5]);
        xlabel('t ms');
        ylabel(['k',label]);
        title('dist:color');
    subplot(3,2,3)
        hold on
        for i=1:nloc1*nloc2
            hue = (did(i)-1)/(n0-1) * (2/3);
            plot(t, RsTar(:,i),'Color',hsv2rgb([hue,sat,val]));
        end
        xlim([t(1)-5,t(end)+5]);
    %     xlabel('t ms');
        ylabel(['Rs',label]);
    subplot(2,2,2)
        plot(mean(tar),std(tar),'*r');
        xlabel(['k',label,'<t>']);
        ylabel(['\sigma(k',label,'<t>)']);
    subplot(3,4,9)
        imagesc(cvTar);
        disp(['maxCV<t>:', num2str(max(cvTar(:)))]);
        disp(['minCV<t>:', num2str(min(cvTar(:)))]);
        axis equal tight
        set(gca,'xaxisLocation','top');
        colorbar;
        ylabel(['loc',label(1)]);
        xlabel(['loc',label(2)]);
    subplot(3,4,10)
        imagesc(cvRsTar);
        disp(['maxRsCV<t>:', num2str(max(cvRsTar(:)))]);
        disp(['minRsCV<t>:', num2str(min(cvRsTar(:)))]);
        axis equal tight
        set(gca,'xaxisLocation','top');
        colorbar;
        ylabel(['loc',label(1)]);
        xlabel(['loc',label(2)]);    
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
    if ~isempty(ext)
        saveas(gcf,[directory,'-k-CVt',label,'.',ext],format);
    end
    %% v
    idt = 1;
    sel = dti(idt) + 100;
    tar = reshape(k(sel,:,:,idt,:),nloc1*nloc2,nv)';
    cvTar = std(squeeze(k(sel,:,:,idt,:)),w,3)'./abs(mean(squeeze(k(sel,:,:,idt,:)),3)');
    RsTar = reshape(Rs(sel,:,:,idt,:),nloc1*nloc2,nv)';
    cvRsTar = std(squeeze(Rs(sel,:,:,idt,:)),w,3)'./abs(mean(squeeze(Rs(sel,:,:,idt,:)),3)');
    figure;
    subplot(3,2,1)
        hold on
        for i=1:nloc1*nloc2
            hue = (did(i)-1)/(n0-1) * (2/3);
            plot(v,tar(:,i),'Color',hsv2rgb([hue,sat,val]));
        end
        xlabel('v mV');
        ylabel(['k',label]);
        title('dist:color');
    subplot(3,2,3)
        hold on
        for i=1:nloc1*nloc2
            hue = (did(i)-1)/(n0-1) * (2/3);
            plot(v, RsTar(:,i),'Color',hsv2rgb([hue,sat,val]));
        end
%         xlabel('v mV');
        ylabel(['Rs',label]);
    subplot(2,2,2)
        plot(mean(tar),std(tar),'*r');
        xlabel(['k',label,'<v>']);
        ylabel(['\sigma(k',label,'<v>)']);
    subplot(3,4,9)
        imagesc(cvTar);
        disp(['maxCV<v>:', num2str(max(cvTar(:)))]);
        disp(['minCV<v>:', num2str(min(cvTar(:)))]);
        axis equal tight
        set(gca,'xaxisLocation','top');
        colorbar;
        ylabel(['loc',label(1)]);
        xlabel(['loc',label(2)]);
    subplot(3,4,10)
        imagesc(cvRsTar);
        disp(['maxRsCV<v>:', num2str(max(cvRsTar(:)))]);
        disp(['minRsCV<v>:', num2str(min(cvRsTar(:)))]);
        axis equal tight
        set(gca,'xaxisLocation','top');
        colorbar;
        ylabel(['loc',label(1)]);
        xlabel(['loc',label(2)]);    
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
    if ~isempty(ext)
        saveas(gcf,[directory,'-k-CVv',label,'.',ext],format);
    end
    %% dt
    sel = 100;
    iv = 2;
    tar0 = zeros(nloc1,nloc2,ndt,nv);
    RsTar0 = zeros(nloc1,nloc2,ndt,nv);
    for idt=1:ndt
        tar0(:,:,idt,iv) = k(dti(idt)+sel,:,:,idt,iv);
        RsTar0(:,:,idt,iv) = Rs(dti(idt)+sel,:,:,idt,iv);
    end
    tar = reshape(tar0(:,:,:,iv),nloc1*nloc2,ndt)';
    cvTar = std(squeeze(tar0(:,:,:,iv)),w,3)'./abs(mean(squeeze(tar0(:,:,:,iv)),3)');
    RsTar = reshape(RsTar0(:,:,:,iv),nloc1*nloc2,ndt)';
    cvRsTar = std(squeeze(RsTar0(:,:,:,iv)),w,3)'./abs(mean(squeeze(RsTar0(:,:,:,iv)),3)');
    figure;
    subplot(3,2,1)
        hold on
        for i=1:nloc1*nloc2
            hue = (did(i)-1)/(n0-1) * (2/3);
            plot(dt,tar(:,i),'Color',hsv2rgb([hue,sat,val]));
        end
        xlabel('dt ms');
        ylabel(['k',label]);
        title('dist:color');
    subplot(3,2,3)
        hold on
        for i=1:nloc1*nloc2
            hue = (did(i)-1)/(n0-1) * (2/3);
            plot(dt, RsTar(:,i),'Color',hsv2rgb([hue,sat,val]));
        end
%         xlabel('dt ms');
        ylabel(['Rs',label]);
        subplot(2,2,2)
        plot(mean(tar),std(tar),'*r');
        xlabel(['k',label,'<dt>']);
        ylabel(['\sigma(k',label,'<dt>)']);
    subplot(3,4,9)
        imagesc(cvTar);
        disp(['maxCV<dt>:', num2str(max(cvTar(:)))]);
        disp(['minCV<dt>:', num2str(min(cvTar(:)))]);
        axis equal tight
        set(gca,'xaxisLocation','top');
        colorbar;
        ylabel(['loc',label(1)]);
        xlabel(['loc',label(2)]);
    subplot(3,4,10)
        imagesc(cvRsTar);
        disp(['maxRsCV<dt>:', num2str(max(cvRsTar(:)))]);
        disp(['minRsCV<dt>:', num2str(min(cvRsTar(:)))]);
        axis equal tight
        set(gca,'xaxisLocation','top');
        colorbar;
        ylabel(['loc',label(1)]);
        xlabel(['loc',label(2)]);  
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
    if ~isempty(ext)
        saveas(gcf,[directory,'-k-CVdt',label,'.',ext],format);
    end
end