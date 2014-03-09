%% MT measurements from one series, sorted by freq then pwr

% Clean and Start up
clear;
load _inits n_mt mtoff_Off_Freq mtoff_alpha_MT
load _recon_coreg mt mtoff

% Pick a slice
z=9;

[mt ind] = nestedSortStruct(mt,{'Off_Freq','alpha_MT'});
mt(end+1).img = mtoff;
mt(end).Off_Freq = mtoff_Off_Freq;
mt(end).alpha_MT = mtoff_alpha_MT;
disp(['Will compare the following sets for slice: ' int2str(z)])
disp('   Offset Freq || alpha_MT')
mt = rmfield(mt, {'pfile','name'}); 
mt = struct2cell(mt);

qMTparams = zeros(size(mt,3),2);
for i=1:size(mt,3)
    qMTparams(i,1)  = mt{1,:,i};
    qMTparams(i,2)  = mt{2,:,i};
    tmp       = mt{3,:,i};
    tt(:,:,i) = tmp(:,:,z);
end

disp(qMTparams)

% Clean and Close up
clear tmp i mt ind mtoff n_mt mt_Off_Freq mtoff_alpha_MT
save _comparisons

%% MT measurements from one series, sorted by pwr then freq

% Clean and Start up
clear;
load _inits n_mt mtoff_Off_Freq mtoff_alpha_MT
load _recon_coreg mt mtoff

% Pick a slice
z=9;

[mt ind] = nestedSortStruct(mt,{'alpha_MT','Off_Freq'});
mt(end+1).img = mtoff;
mt(end).Off_Freq = mtoff_Off_Freq;
mt(end).alpha_MT = mtoff_alpha_MT;
disp(['Will compare the following sets for slice: ' int2str(z)])
disp('   Offset Freq || alpha_MT')
mt = rmfield(mt, {'pfile','name'}); 
mt = struct2cell(mt);

qMTparams = zeros(size(mt,3),2);
for i=1:size(mt,3)
    qMTparams(i,1)  = mt{1,:,i};
    qMTparams(i,2)  = mt{2,:,i};
    tmp       = mt{3,:,i};
    tt(:,:,i) = tmp(:,:,z);
end

disp(qMTparams)

% Clean and Close up
clear tmp i mt ind mtoff n_mt mt_Off_Freq mtoff_alpha_MT
save _comparisons

%% MTR from one series, sorted by freq then pwr

% Clean and Start up
clear;
load _inits n_mt
load _recon_coreg mt mtoff

% Pick a slice
z=9;

[mt ind] = nestedSortStruct(mt,{'Off_Freq','alpha_MT'});
disp(['Will compare the following sets for slice: ' int2str(z)])
disp('   Offset Freq || alpha_MT')
mt  = rmfield(mt, {'pfile','name'}); 
mt  = struct2cell(mt);
pts = find(mtoff);

qMTparams = zeros(size(mt,3),2);
for i=1:size(mt,3)
    qMTparams(i,1)  = mt{1,:,i};
    qMTparams(i,2)  = mt{2,:,i};
    tmp       = mt{3,:,i};
    tmp(pts)  = tmp(pts)./mtoff(pts); % find MTR
    tt(:,:,i) = tmp(:,:,z);
end

disp(qMTparams)

% Clean and Close up
clear tmp i mt ind mtoff n_mt pts
save _comparisons

%% MTR from one series, sorted by pwr then freq

% Clean and Start up
clear;
load _inits n_mt
load _recon_coreg mt mtoff

% Pick a slice
z=9;

[mt ind] = nestedSortStruct(mt,{'alpha_MT','Off_Freq'});
disp(['Will compare the following sets for slice: ' int2str(z)])
disp('   Offset Freq || alpha_MT')
mt  = rmfield(mt, {'pfile','name'}); 
mt  = struct2cell(mt);
pts = find(mtoff);

qMTparams = zeros(size(mt,3),2);
for i=1:size(mt,3)
    qMTparams(i,1)  = mt{1,:,i};
    qMTparams(i,2)  = mt{2,:,i};
    tmp       = mt{3,:,i};
    tmp(pts)  = tmp(pts)./mtoff(pts); % find MTR
    tt(:,:,i) = tmp(:,:,z);
end

disp(qMTparams)

% Clean and Close up
clear tmp i mt ind mtoff n_mt pts
save _comparisons

%% qMT params, AS

% Clean and Start up
clear, load _data fv2

imsc(131, fv2(:,:,1), [0 .20]);    clb, title('Bound pool fraction (f)')
imsc(132, fv2(:,:,2), [0  10]);    clb, title('Exchange rate constant (k)')
imsc(133, fv2(:,:,3), [0  20e-6]); clb, title('T2 of bound pool (t2b)')
% Grayscale colormap for Rick and Wally
colormap('gray')

%% qMT params, PM

