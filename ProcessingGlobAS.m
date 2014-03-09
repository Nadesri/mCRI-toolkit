% PROCESSINGGLOB Global Processing using Alexey's MT Fit
% Assuming standardized naming, will process AFI, SPGR, and qMT for the
%  current directory.
%
% Nade Sritanyaratana
% University of Wisconsin, Madison
% 11/3/10

clear
load _inits
clear pfile pdir
load _recon_coreg

%% FAM Processing
fam = afi_standard(afi0,afi1,n);
for i=1:size(fam,3)
    fam_s(:,:,i) = conv2(fam(:,:,i), kernel, 'same');
end

%% R1 Processing
data_spgr = zeros([res.x*res.y*res.z n_spgr]);
for i=1:n_spgr
    data_spgr(:,i)  = spgr(i).img(:);
end
[pd r1] = t1_fit_spgr_IRLS_afi(data_spgr, alpha_spgr, spgrTR, ones([res.x*res.y*res.z 1]), Ni);
r1 = reshape(r1, [res.x res.y res.z]);
pd = reshape(pd, [res.x res.y res.z]);
t1 = 1./r1;

%% qMT Processing Pt.1

mtoffz = mtoff(:,:,z);
mtz = mt; % Give mtz the same struct fields as mt
for i=1:n_mt
    mtz(i).img = mt(i).img(:,:,z);
end

r1z = r1(:,:,z);
fam_sz = fam_s(:,:,z);

%% Fat Suppression Based Segmentation
mask = (mtoffz-mtz(1).img); % we already happen to know that mtz(1) has the 
                             % most MT effect. We may have to fix this
                             % later..
mask = (mask>0).*mask; mask = mask./max(mask(:));
mask = mask>fatsatthresh;

%% Select Cartilage Segmentation
mask4 = zeros(size(mask));

figure
for i=1:length(z)
    imsc(mimsc(mask(:,:,i)))
    mask4(:,:,i) = roipoly;
end

mask = mask.*mask4;

%% MTR Based Segmentation
pts = find(mtoffz);

mtr = mt; % Give mtr the same struct fields as mt
for i=1:n_mt
    mtr(i).img = zeros(size(mtoffz));
    mtr(i).img(pts) = mtz(i).img(pts)./mtoffz(pts);
end

mask2 = ones(size(mask));
% MTR<1 segmentation
for i=1:n_mt
    mask2 = mask2.*(mtr(i).img<1);
end

% MTR monotonic segmentation
    % look for MT powers.
    alpha_MT_unique=[];
    for i=1:length(alpha_MT)
        if sum(alpha_MT(i)==alpha_MT_unique)==0
            alpha_MT_unique = [alpha_MT_unique alpha_MT(i)];
        end
    end
    
mask3 = ones(size(mask));
for i=1:length(alpha_MT_unique)
    % Gather only the MTRs that share the same alpha_MT
    k=1;
    for j=1:n_mt
        if mtr(j).alpha_MT==alpha_MT_unique(i)
            mtr1pwr(k) = mtr(j);
            k=k+1;
        end
    end
    
    % Sort the MTRs according to their frequency in ascending order
    mtr1pwr = sortStruct(mtr1pwr, 'Off_Freq');
 
    % Find the difference between 2 sets and ensure that it is above zero
    % (monotonically increasing)
    for j=1:length(mtr1pwr)-1
        mtrdiff = mtr1pwr(j+1).img-mtr1pwr(j).img;
        mask3 = mask3.*(mtrdiff>0);
    end
end

mask = mask.*mask2.*mask3;

% clean up
clear alpha_MT_unique pair1 pair2 pair1initflag pair2initflag mtrdiff mtr1pwr i j

%% Select pts for MT analysis

pts = find(mask);

%% Currently Unused Code (Monotonic Segmentation)

%     % look for MT Offsets.
%     MT_Off_Freq_unique=[];
%     for i=1:length(MT_Off_Freq)
%         if sum(MT_Off_Freq(i)==MT_Off_Freq_unique)==0
%             MT_Off_Freq_unique = [MT_Off_Freq_unique MT_Off_Freq(i)];
%         end
%     end
% for i=1:length(alpha_MT_unique)
%     pair1 = zeros(size(mask)); pair2 = zeros(size(mask)); 
%     pair1initflag = 0; pair2initflag = 0;% reset MTR pair
%     for j=1:n_mt
%         if mtr(j).alpha_MT==alpha_MT_unique(i)
%             if ~pair1initflag
%                 pair1 = mtr(j).img;
%                 pair1initflag = 1;
%             elseif ~pair2initflag
%                 pair2 = mtr(j).img;
%                 pair2initflag = 1;
%             end
%             if pair1initflag&&pair2initflag
%                 if pair1.Off_Freq>pair2.Off_Freq
%                     high = pair1;
%                     low  = pair2;
%                 else
%                     high = pair2;
%                     low  = pair1;
%                 end
%                 mtrdiff = high.img-low.img;
%                 mask3 = mask3.*((mtrdiff)>0);
%                 pair1 = pair2; %move the pair up to the next frequency
%                 pair2initflag = 0; pair2 = zeros(size(pair1));
%             end
%         end
%     end
% end

%% qMT Processing Pt. 2

data_MT = zeros([length(pts) n_mt]);
for i=1:n_mt
    data_MT(:,i) = mtr(i).img(pts);
end
B1 = fam_sz(pts)/nomflip-1;

%  Fitting MT-weighted data to a qMT model
%  Function signature:
%       [fv, rnrm]=mt_model_fit(data, mtoff, r1, flip0, mtflip, tm, TR, pulse_type, line_shape, b1err, T2b, ig);
%
%   INPUTS:
%       data [Np x Nms]   - MT weighted data array (Np and Nms - number of points and measurements, respectively)
%       mtoff [1 x Nms]   - values of MT pulse frequence offsets in Hz
%       r1 [Np x 1]       - vector of relaxation rates (1/T1)
%       flip0 [1 x 1]     - excitation flip angle of MT sequence (in degrees)
%       mtflip [1 x Nms]  - flip angle of MT saturation pulse
%       tm [1 x 1]        - the width of MT saturation pulse in
%       TR [1 x 1]        - repetition time of MT SPGR sequence
%       pulse_type        -  type of MT pulse ('sinc' or 'fermi')
%       line_shape        -  shape of saturation line of bound pool ('gaussian' or
%                                  'superLorentian)
%       mode -    0: calculation of SS magnetization; data array assumed to
%                           contain rows of corresponding f,k,T2b parameters
%                       1:  standard multipoint fit
%                       2: two point mapping of bound pool fraction
%
%       OPTIONAL Parameters:
%       b1err [Np x 1]   -  vector of relative B1 error (OPTIONAL)
%
%       b0err [Np x 1]   -  vector of off-resonance B0 (in Hz) OPTIONAL
%
%       T2b   -   value of T2 of bound pool (in sec); if empty, the parameter is
%               obtained from the data fit (less precise, but more accurate estimates
%               of f and k) (OPTIONAL)
%
%
%       ig [3(2)x1]   -   vector of initial guess [f0 k0 (T2b0)] (OPTIONAL)
%
%   OUTPUTS:
%       fv [Np x 3(2)]  -   vector of fitted parameter with column order [f k (T2b)]
%       rnrm [Np x 1]   -   residual values of the fits
%
%
%   Reference: Yarnykh VL, MRM (2002)
%
%   Alexey Samsonov
%   University of Wisconsin, Madison
%   August 2005 - September 2007
%
%   v1.1 21-Oct-2009
%   Modified by Samuel A. Hurley
%            v1.0 - include B0 correction
[fv, rnrm]=mt_model_fit(data_MT, MT_Off_Freq, r1z(pts), alpha_0, alpha_MT, t_m, TR, pulse_type, line_shape, B1);

fv2 = zeros([res.y res.x length(z) size(fv,2)]);
for i=1:size(fv,2)
    im  = zeros([res.y res.x length(z)]);
    tmp = squeeze(fv(:,i));
    im(pts) = tmp;
    fv2(:,:,:,i) = im;
end

rnrm2 = zeros([res.y res.x length(z)]);
rnrm2(pts) = rnrm;

%% RNRM and k-Based Segmentation

fv2(:,:,:,1) = fv2(:,:,:,1).*(rnrm2<rnrmthresh).*(fv2(:,:,:,2)<kthresh);
fv2(:,:,:,2) = fv2(:,:,:,2).*(rnrm2<rnrmthresh).*(fv2(:,:,:,2)<kthresh);
fv2(:,:,:,3) = fv2(:,:,:,3).*(rnrm2<rnrmthresh).*(fv2(:,:,:,2)<kthresh);

%% Clean and Close up

clear B1 data_MT tmp im fv rnrm i mask2 mask3
save _data fam fam_s r1 t1 pd mask fv2 rnrm2 mtr