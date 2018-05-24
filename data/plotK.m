directory = 'active';
fid = fopen([directory,'/cKRsK.bin'],'r');
distanceMat = 'activeBranchDis.mat';
run_nt=2401;
nlocE=6;
nlocI=3;
tstep = 0.1;
v = [-74,-70,-66,-62];
nv = length(v);
dt = [0,10,20,30];
dti = round(dt/tstep);
ndt = length(dt);
range = 1:nlocE;
diag = (range-1)*nlocE + range;
nt = 4;
tRange = round(linspace(100,300,nt));
mesTrange = round((10:tstep:30)/tstep);

data = fread(fid,[run_nt,nv*nlocE*nlocE*ndt],'double');
kEE = reshape(data,run_nt,nlocE,nlocE,ndt,nv);
data = fread(fid,[run_nt,nv*nlocE*nlocE*ndt],'double');
RsEE = reshape(data,run_nt,nlocE,nlocE,ndt,nv);

data = fread(fid,[run_nt,nv*nlocE*nlocI*ndt],'double');
kEI = reshape(data,run_nt,nlocE,nlocI,ndt,nv);
data = fread(fid,[run_nt,nv*nlocE*nlocI*ndt],'double');
RsEI = reshape(data,run_nt,nlocE,nlocI,ndt,nv);

data = fread(fid,[run_nt,nv*nlocI*nlocE*ndt],'double');
kIE = reshape(data,run_nt,nlocI,nlocE,ndt,nv);
data = fread(fid,[run_nt,nv*nlocI*nlocE*ndt],'double');
RsIE = reshape(data,run_nt,nlocI,nlocE,ndt,nv);

data = fread(fid,[run_nt,nv*nlocI*nlocI*ndt],'double');
kII = reshape(data,run_nt,nlocI,nlocI,ndt,nv);
data = fread(fid,[run_nt,nv*nlocI*nlocI*ndt],'double');
RsII = reshape(data,run_nt,nlocI,nlocI,ndt,nv);

fclose(fid);

fid = fopen('cKRsK0.bin','r');

data = fread(fid,[run_nt,nlocE*nlocE*ndt],'double');
kEE0 = reshape(data,run_nt,nlocE,nlocE,ndt);
data = fread(fid,[run_nt,nlocE*nlocE*ndt],'double');
RsEE0 = reshape(data,run_nt,nlocE,nlocE,ndt);

data = fread(fid,[run_nt,nlocE*nlocI*ndt],'double');
kEI0 = reshape(data,run_nt,nlocE,nlocI,ndt);
data = fread(fid,[run_nt,nlocE*nlocI*ndt],'double');
RsEI0 = reshape(data,run_nt,nlocE,nlocI,ndt);

data = fread(fid,[run_nt,nlocI*nlocE*ndt],'double');
kIE0 = reshape(data,run_nt,nlocI,nlocE,ndt);
data = fread(fid,[run_nt,nlocI*nlocE*ndt],'double');
RsIE0 = reshape(data,run_nt,nlocI,nlocE,ndt);

data = fread(fid,[run_nt,nlocI*nlocI*ndt],'double');
kII0 = reshape(data,run_nt,nlocI,nlocI,ndt);
data = fread(fid,[run_nt,nlocI*nlocI*ndt],'double');
RsII0 = reshape(data,run_nt,nlocI,nlocI,ndt);
fclose(fid);
%%
% 2 outliers:
%(nlocE,nlocE):(3,3) (5,5); nlocE*nlocE: 15, 29
exploreK(kEE,v,tRange,dt,dti,nv,ndt,nlocE,nlocE,run_nt,tstep,'EE',distanceMat);
meshK(kEE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocE,false,'EE');
% meshK(kEE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocE,true,'EE');

exploreK(kEI,v,tRange,dt,dti,nv,ndt,nlocE,nlocI,run_nt,tstep,'EI',distanceMat);
meshK(kEI,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocI,false,'EI');
% meshK(kEI,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocI,true,'EI');
% 
exploreK(kIE,v,tRange,dt,dti,nv,ndt,nlocI,nlocE,run_nt,tstep,'IE',distanceMat);
meshK(kIE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocE,false,'IE');
% meshK(kIE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocE,true,'IE');
% 
exploreK(kII,v,tRange,dt,dti,nv,ndt,nlocI,nlocI,run_nt,tstep,'II',distanceMat);
meshK(kII,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocI,false,'II');
% meshK(kII,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocI,true,'II');

