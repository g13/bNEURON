function plotK(directory,distanceMat,EE,EI,IE,II,ext)
if nargin < 7
    ext = '';
    if nargin < 6
        II = true;
        if nargin < 5
            IE = true;
            if nargin < 4
                EI = true;
                if nargin < 3
                    EE = true;
                end
            end
        end
    end
end
fid = fopen([directory,'/cKRsK.bin'],'r');
run_nt=2401;
nlocE=6;
nlocI=3;
tstep = 0.1;
v = [-74,-70,-66,-62];
nv = length(v);
dt = [0,10,20,30];
dti = round(dt/tstep);
ndt = length(dt);
nt = 4;
tRange = round(linspace(100,300,nt));
mesTrange = round((10:tstep:30)/tstep);
% 
% timeline:  t1-----------t2--------------->
% location: loc2         loc1
% k(nt,loc1,loc2,dt,v);
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

% fid = fopen([directory,'/cKRsK0.bin'],'r');
% 
% data = fread(fid,[run_nt,nlocE*nlocE*ndt],'double');
% kEE0 = reshape(data,run_nt,nlocE,nlocE,ndt);
% data = fread(fid,[run_nt,nlocE*nlocE*ndt],'double');
% RsEE0 = reshape(data,run_nt,nlocE,nlocE,ndt);
% 
% data = fread(fid,[run_nt,nlocE*nlocI*ndt],'double');
% kEI0 = reshape(data,run_nt,nlocE,nlocI,ndt);
% data = fread(fid,[run_nt,nlocE*nlocI*ndt],'double');
% RsEI0 = reshape(data,run_nt,nlocE,nlocI,ndt);
% 
% data = fread(fid,[run_nt,nlocI*nlocE*ndt],'double');
% kIE0 = reshape(data,run_nt,nlocI,nlocE,ndt);
% data = fread(fid,[run_nt,nlocI*nlocE*ndt],'double');
% RsIE0 = reshape(data,run_nt,nlocI,nlocE,ndt);
% 
% data = fread(fid,[run_nt,nlocI*nlocI*ndt],'double');
% kII0 = reshape(data,run_nt,nlocI,nlocI,ndt);
% data = fread(fid,[run_nt,nlocI*nlocI*ndt],'double');
% RsII0 = reshape(data,run_nt,nlocI,nlocI,ndt);
% fclose(fid);
%%
% 2 outliers:
%(nlocE,nlocE):(3,3) (5,5); nlocE*nlocE: 15, 29
if ext
    save = true;
end
if EE
    exploreK(kEE,RsEE,v,tRange,dt,dti,nv,ndt,nlocE,nlocE,run_nt,tstep,'EE',distanceMat,ext,directory);
    meshK(kEE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocE,false,'EE',save,directory);
    meshK(kEE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocE,true,'EE',save,directory);
end
if EI
    exploreK(kEI,RsEI,v,tRange,dt,dti,nv,ndt,nlocE,nlocI,run_nt,tstep,'EI',distanceMat,ext,directory);
    meshK(kEI,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocI,false,'EI',save,directory);
    meshK(kEI,mesTrange,tstep,dt,dti,ndt,v,nv,nlocE,nlocI,true,'EI',save,directory);
end
if IE
    exploreK(kIE,RsIE,v,tRange,dt,dti,nv,ndt,nlocI,nlocE,run_nt,tstep,'IE',distanceMat,ext,directory);
    meshK(kIE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocE,false,'IE',save,directory);
    meshK(kIE,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocE,true,'IE',save,directory);
end
if II
    exploreK(kII,RsII,v,tRange,dt,dti,nv,ndt,nlocI,nlocI,run_nt,tstep,'II',distanceMat,ext,directory);
    meshK(kII,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocI,false,'II',save,directory);
    meshK(kII,mesTrange,tstep,dt,dti,ndt,v,nv,nlocI,nlocI,true,'II',save,directory);
end
