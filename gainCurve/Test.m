theme = 'gainCurve-test';
nInput = 10;
%inputFn = '';
inputFn = 'inputTable.bin';
%levelsFn = '';
levelsFn = 'levels.bin';
pVar = true;
%pVar = false;
%tVar = true;
tVar = false;
extendVar = true;
%extendVar = false;
dtVarLevels = true;
tVarLevels = true;
irregInputLevels = true;
runTime0 = 1;
dt0 = 0.001;
if isempty(levelsFn)
    assert(~irregInputLevels && ~tVarLevels && ~dtVarLevels);
    nLevel = 10;
    inputLinspace = [1,nLevel,10]
    runTime = runTime0;
    dt = dt0;
else
    assert(irregInputLevels || tVarLevels || dtVarLevels);
    fid = fopen(levelsFn,'w');
    if irregInputLevels
        inputLevels = [1,2,4];
        fwrite(fid,inputLevels,'double');
        nLevel = length(inputLevels);
    end
    if tVarLevels
        runTime = [0.1,0.2,0.4];
        fwrite(fid,runTime,'double');
        nLevel = length(runTime);
        assert(nLevel==length(inputLevels)||length(inputLevels)==1);
    else
        runTime = runTime0;
    end
    if dtVarLevels
        dt = [0.001, 0.001, 0.001];
        fwrite(fid,dt,'double');
        nLevel = length(dt);
        assert((nLevel == length(runTime) || length(runTime) == 1) && (nLevel == length(inputLevels) || length(inputLevels) == 1));
    else
        dt = dt0;
    end
    fclose(fid);
    levelsFopt = [' --inputLevelsFn ', levelsFn];
end
if ~isempty(inputFn)
    fid = fopen(inputFn,'w');
    fclose(fid);
    fid = fopen(inputFn,'a');
    if ~pVar
        p.w = 1;
        p.b = 0.5;
        for i=1:nLevel
            if (~tVarLevels && ~dtVarLevels)
                t = 0:dt:runTime;
                % p %parameters can still be changed with extendVar
            else 
                if tVarLevels && dtVarLevels
                    t = 0:dt(i):runTime(i);
                else
                    if ~tVarLevels
                        t = 0:dt(i):runTime;
                    else
                        t = 0:dt:runTime(i);
                    end
                end
            end
            input = inputF(p,t)
            fwrite(fid,input,'double');
        end
    else
        p.w = 1;
        p.b = 0.5;
        if ~tVar
            if extendVar
                for i=1:nLevel
                    posR = ((1:nInput)+nLevel)/(nInput+nLevel);
                    fwrite(fid,posR,'double');
                end
            else
                posR = 1:nInput;
                fwrite(fid,posR,'double');
            end
        else
            for i=1:nInput     
                input = inputF(p,t);
                fwrite(fid,input,'double');
            end
        end
    end
    fclose(fid);
    inputFopt = [' --inputFn ', inputFn];
end
if pVar, spVar = ' --pVar'; else, spVar = ''; end
if tVar, stVar = ' --tVar'; else, stVar = ''; end
if extendVar, sextendVar = ' --extendVar'; else sextendVar = ''; end
if tVarLevels
    stVarLevels = ' --tVarLevels';
    runTimeOpt = '';
else 
    stVarLevels = ''; 
    runTimeOpt = [' --runTime ', num2str(runTime)];
end
if dtVarLevels
    sdtVarLevels = ' --dtVarLevels';
    dtOpt = '';
else 
    sdtVarLevels = '';
    dtOpt = [' --dt ', num2str(dt)];
end
if irregInputLevels, sirregInputLevels = ' --irregInputLevels'; else sirregInputLevels = ''; end
themeOpt = [' -m ', theme];
snInput = [' --nInput ', num2str(nInput)];
cmd = ['./gainCurve', themeOpt, snInput, spVar, stVar, sextendVar, sirregInputLevels, sdtVarLevels, stVarLevels, inputFopt, levelsFopt]
ret = system(cmd)
