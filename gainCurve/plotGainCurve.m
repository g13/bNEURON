function plotGainCurve(inputFn, readSubthreshold = false)
    dimsFn = 'readoutDimension.bin';
    p = read_cfg(inputFn);
    if ~isstruct(p)
        return
    end
    fdr = struct2cell(dir());
    % name, fdr, date, bytes, isdir, datenum
    filter = [1,5];
    fdr = fdr(filter,:);
    fnstr = strjoin(fdr(1,:));
    load(libFile,'tstep','n');
    dataPrefixPattern = ['\s',p.theme,'-s\d*\-\d*\-(Data|jND|Raster|tIn)\.bin'];
    dataPrefix = regexp(fnstr, dataPrefixPattern,'once')
    dimsFid = fopen(dimsFn,'r');
    if dimsFid
        nTrial = fread(dimsFid, 1, 'int');
        fclose(dimsFid);
    end
    RasterFn = [dataPrefix,'Raster.bin'];
    if readSubthreshold
        dataFn = [dataPrefix,'Data.bin'];
        jNDFn = [dataPrefix,'jND.bin'];
        tInFn = [dataPrefix,'tIn.bin'];
        if isfield(p,'reformatInputFn')
            dimsFn = p.reformatInputFn;
        end
        dimsFid = fopen(dimsFn,'r');
        if dimsFid
            nTrial = fread(dimsFid, 1, 'int');
            nDimSim = fread(dimsFid, nTrial, 'unsigned long');
            nDim = fread(dimsFid, nTrial, 'unsigned long');
            tstep = fread(dimsFid, nTrial, 'double');
            fclose(dimsFid);
        end
        datafid = fopen(dataFn,'r');
        dendV = cell(nTrial,1);
        simV = cell(nTrial,1);
        biV = cell(nTrial,1);
        if datafid
            for i=1:nTrial
                simV = fread(datafid, nDimSim(i), 'double');
                biV = fread(datafid, nDim(i), 'double');
                liV = fread(datafid, nDim(i), 'double');
                biV0 = fread(datafid, nDim(i), 'double');
                for j=1:n
                    dendV = fread(datafid, nDimSim(i), 'double');
                end
            end
            data = fread(datafid, nDim(i), 'double');
            fclose(dimsFid);
        end
        jNDfid = fopen(jNDFn,'r');
        if jNDfid
            jbSize = fread(jNDfid,1,'uint64');
            jbt = fread(jNDfid,jbSize,'double');
            jbv = fread(jNDfid,jbSize,'double');
            jbnCross = fread(jNDfid,1,'uint64');
            for i=1:jbnCross
                tmpSize = fread(jNDfid, 1, 'uint64');
                jbCrossT = fread(jNDfid, tmpSize, 'double');
                jbCrossV = fread(jNDfid, tmpSize, 'double');
            end
            jlSize = fread(jNDfid,1,'uint64');
            jlt = fread(jNDfid,jlSize,'double');
            jlv = fread(jNDfid,jlSize,'double');
            jlnCross = fread(jNDfid,1,'uint64');
            for i=1:jbnCross
                tmpSize = fread(jNDfid, 1, 'uint64');
                jlCrossT = fread(jNDfid, tmpSize, 'double');
                jlCrossV = fread(jNDfid, tmpSize, 'double');
            end
        end
    end
end
