nInput = 3;
inputFn = 'inputTable.bin';
levelsFn = 'levels.bin';
nLevel = 3;
fid = fopen(levelsFn,'w');
inputLevels = [10,20,40];
runTime = [10,20,30];
dt = [0.1,0.02,0.01];

assert(length(inputLevels) == nLevel);
fwrite(fid,inputLevels,'double');
assert(length(runTime) == nLevel);
fwrite(fid,runTime,'double');
assert(length(dt) == nLevel);
fwrite(fid,dt,'double');
fclose(fid);

fid = fopen(inputFn,'w');
inputs = rand(nInput,nLevel);
fwrite(fid,inputs,'double');