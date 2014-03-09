
%% Initialize
clear; load _recon; load _inits res fov Teps N n_spgr n_mt
cd ANALYZE

vsz = [fov.x/res.x fov.y/res.y fov.z];

%% Segment background noise

afi0 = imseg_hist_click(afi0,Teps,N);
afi1 = imseg_hist_click(afi1,Teps,N);
for i=1:n_spgr
    spgr(i).img = imseg_hist_click(spgr(i).img,Teps,N);
end
mtoff = imseg_hist_click(mtoff,Teps,N);
for i=1:n_mt
    mt(i).img = imseg_hist_click(mt(i).img,Teps,N);
end

% clean up
clear Teps N

%% Write Matlab variables to ANALYZE format

write_analyze(afi0, vsz, 'afi0', 0,1);
write_analyze(afi1, vsz, 'afi1', 0,1);
for i=1:n_spgr
    write_analyze(spgr(i).img, vsz, spgr(i).name,0,1);
end
write_analyze(mtoff, vsz, 'mtoff', 0,1);
for i=1:n_mt
    write_analyze(mt(i).img, vsz, mt(i).name,0,1);
end

%% Run Flirt coregistration

refName = 'mtoff';
cmdprfix   = '!flirt -searchcost mutualinfo -cost mutualinfo -bins 256 -searchrx -15 15 -searchry -15 15 -searchrz -15 15 -dof 6 -interp sinc ';

inName = 'afi0';
cmnd   = [cmdprfix '-in ' inName ' -ref ' refName ' -out '  inName '_out -omat ' inName '_rprm.txt'];
eval(cmnd);

inName = 'afi1';
cmnd   = ['!flirt -in ' inName ' -ref ' refName ' -out '  inName '_out -applyxfm -init afi0_rprm.txt'];
eval(cmnd)

for i = 1:n_mt
    inName = mt(i).name;
    cmnd   = [cmdprfix '-in ' inName ' -ref ' refName ' -out '  inName '_out -omat ' inName '_rprm.txt'];
    eval(cmnd);
end

for i = 1:n_spgr
    inName = spgr(i).name;
    cmnd   = [cmdprfix '-in ' inName ' -ref ' refName ' -out '  inName '_out -omat ' inName '_rprm.txt'];
    eval(cmnd);
end

%% Read Coregistered ANALYZE files to Matlab variables

afi0       = read_analyze('afi0_out' ,1,1);
afi1       = read_analyze('afi1_out' ,1,1);
for i=1:n_spgr
    spgr(i).img = read_analyze([spgr(i).name '_out'] ,1,1);
end

for i=1:n_mt
    mt(i).img = read_analyze([mt(i).name '_out'],1,1);
end

%% Clean and close up

clear inName cmnd cmdprfix refName i res fov vsz
cd .., save _recon_coreg;