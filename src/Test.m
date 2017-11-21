nInput = 10;
inputFn = 'inputTable.bin';
fid = fopen(inputFn,'w');
fclose(fid);
fid = fopen(inputFn,'a');
nLevel = 10;
runTime = ones(nLevel,1);
dt = ones(nLevel,1)*0.001;
inputLevel = ones(nLevel,1);
irregInputLevels = true;
pVar = true;
tVar = true;
extendVar = true;
dtVarLevels = true;
tVarLevels = true;
for i=0:nLevel
    t = 0:dt(i):runTime(i);
    input = inputF(p,t);
    fwrite(fid,input);
end
if pVar, spVar = '--pVar'; end
if tVar, stVar = '--tVar'; end
if extendvar, sextendVar = ' --extendVar'; end
if tVarLevels, stVarLevels = ' --tVarLevels'; end
if dtVarLevels, sdtVarLevels = ' --dtVarLevels'; end
if irregLevels, sirregInputLevels = ' --irregInputLevels'; end
fclose(fid);
cmd = 'test' + spVar + stVar + extendVar + sirregInputLevels + sdtVarLevels + stVarLevels;
system(cmd);
