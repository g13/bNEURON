function setInput(write_inputTable, write_levels,cfgfile)
    nLevel = 1;
    inputLevels = [10,20,30,40];
    runTime = [100,200,300];
    dt = [0.1,0.1,0.01];
    inputs = [ones(6,1);zeros(2,1)+0.5];
    dt = zeros(nLevel,1);
    runTime = zeros(nLevel,1);
    inputLevels = zeros(nLevel,1);
    %inputs = ones(p.nInput,nLevel);
    if nargin < 3
        cfgfile = 'gainCurve.cfg';
    end
    p = read_cfg(cfgfile)
    assert(p.nInput == size(inputs,1));
    if ~write_inputTable 
        disp('nothing to write for inputTable')
        if p.pVar || p.tVar
            disp('although cfg file demand inputTable file');
        end
        return
    end
    if ~write_levels
        disp('nothing to write for levels')
        if p.dtVarLevels || p.tVarLevels || p.irregInputLevels
            disp('although cfg file demand levels data file');
        end
    end
    inputFn = 'inputTable.bin';
    levelsFn = 'levels.bin';
    if write_inputTable
        if p.pVar || p.tVar
            if sum(inputs(:)) == 0 
                disp('inputs data demanded, but not supplied');
                return
            end
            if (p.tVar && p.extendVar) || (p.pVar && p.extendVar) || (p.pVar && p.tVar)
                if length(size(inputs))<2 
                    disp('lack another dimension of data');
                    return
                end
            end
        end
        fid = fopen(inputFn,'w');
        fwrite(fid,inputs,'double');
        fclose(fid);
    end
    if write_levels
        fid = fopen(levelsFn,'w');
        if p.irregInputLevels
            if sum(inputLevels)==0
                disp('irregInputLevels demanded, but not supplied');
                return
            end
            assert(length(inputLevels) == nLevel);
            fwrite(fid,inputLevels,'double');
        end
        if p.tVarLevels
            if sum(runTime)==0
                disp('tVarLevels demanded, but not supplied');
                return
            end
            assert(length(runTime) == nLevel);
            fwrite(fid,runTime,'double');
        end
        if p.dtVarLevels
            if sum(dt)==0
                disp('dtVarLevels demanded, but not supplied');
                return
            end
            assert(length(dt) == nLevel);
            fwrite(fid,dt,'double');
        end
        fclose(fid);
    end
end
