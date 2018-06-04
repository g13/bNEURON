function meshK(k,tRange,tstep,dt,dti,ndt,v,nv,nloc1,nloc2,normTar,label,save,directory)
    if nargin < 11
        normTar = true;
    end
    figure;
    %% v averaged
    subplot(2,2,1)
    hold on
    tar0 = zeros(length(tRange),nloc1*nloc2,ndt,nv);
    [X,Y] = meshgrid(dt,tRange*tstep);
    for idt = 1:ndt
        tar0(:,:,idt,:) = reshape(k(dti(idt)+tRange,:,:,idt,:),length(tRange),nloc1*nloc2,nv);
    end
    for i = 1:nloc1*nloc2
        tar = mean(squeeze(tar0(:,i,:,:)),3);
        if normTar
            tar = (tar - mean(tar(:)))./std(tar(:));
            %tar = (tar - repmat(mean(tar),[length(tRange),1]))./repmat(std(tar),[length(tRange),1]);
        end
        Z = tar;
        mesh(X,Y,Z);
    end
    xlabel('dt');
    ylabel('t');
    if normTar
        title(['norm. k',label,'<v>']);
    else
        title(['k',label,'<v>']);
    end
    %% dt averged
    subplot(2,2,2)
    hold on
    [X,Y] = meshgrid(v,tRange*tstep);
    for idt = 1:ndt
        tar0(:,:,idt,:) = reshape(k(dti(idt)+tRange,:,:,idt,:),length(tRange),nloc1*nloc2,nv);
    end
    for i = 1:nloc1*nloc2
        tar = mean(squeeze(tar0(:,i,:,:)),2);
        if normTar
            tar = (tar - mean(tar(:)))./std(tar(:));
            %tar = (tar - repmat(mean(tar),[length(tRange),1]))./repmat(std(tar),[length(tRange),1]);
        end
        Z = squeeze(tar);
        mesh(X,Y,Z);
    end
    xlabel('v');
    ylabel('t');
    if normTar
        title(['norm. k',label,'<dt>']);
    else
        title(['k',label,'<dt>']);
    end
    %% t averaged
    subplot(2,2,3)
    hold on
    [X,Y] = meshgrid(v,dt);
    for idt = 1:ndt
        tar0(:,:,idt,:) = reshape(k(dti(idt)+tRange,:,:,idt,:),length(tRange),nloc1*nloc2,nv);
    end
    for i = 1:nloc1*nloc2
        tar = mean(squeeze(tar0(:,i,:,:)));
        if normTar
            tar = (tar-mean(tar(:)))./std(tar(:));
        end
        Z = squeeze(tar);
        mesh(X,Y,Z);
    end
    xlabel('v');
    ylabel('dt');
    if normTar
        title(['norm. k',label,'<t>']);
    else
        title(['k',label,'<t>']);
    end
    if save
        saveas(gcf,[directory,'-meshK',label,'.fig'],'fig');
    end
end