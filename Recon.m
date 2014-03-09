clear;
load _inits afi_pfile spgr mt mtoff_pfile n_spgr n_mt pdir
cd(pdir)

%% Reconstruction

afi0 = recon_vipr(afi_pfile, 'afi0', 0, 2.2);
afi1 = recon_vipr(afi_pfile, 'afi1', 1, 2.2);

for i=1:n_spgr
    spgr(i).name = ['0' int2str(spgr(i).alpha)];
    spgr(i).name = ['spgr' spgr(i).name(end-1:end)];
    spgr(i).img  = recon_vipr(spgr(i).pfile, spgr(i).name, -1, 2.2);
end

mtoff = recon_vipr(mtoff_pfile, 'mtoff', -1, 2.2);

for i=1:n_mt
    mt(i).name = ['mt' int2str(mt(i).Off_Freq) int2str(mt(i).alpha_MT)];
    mt(i).img = recon_vipr(mt(i).pfile,  mt(i).name,  -1, 2.2);
end

clear i n_spgr n_mt pdir mtoff_pfile afi_pfile
cd ..
save _recon;