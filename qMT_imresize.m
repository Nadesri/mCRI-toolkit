% qMT Image Resizing
% Can only be run after creating all raw images (See Recon.m)

%% Initialize

clear
load _recon mt spgr mtoff afi0 afi1
load _inits res n_mt n_spgr

%% Resize Images

afi0r   = zeros([res.y res.x res.z]); afi1r   = afi0r;
for i=1:n_spgr
    spgr_r{i} = zeros([res.y res.x res.z]);
end
mtoffr = zeros([res.y res.x res.z]);
for i=1:n_mt
    mt_r{i} = zeros([res.y res.x res.z]);
end

for i=1:res.z
      afi0r(:,:,i) = imresize(  afi0(:,:,i), [res.y res.x]);
      afi1r(:,:,i) = imresize(  afi1(:,:,i), [res.y res.x]);
      for j=1:n_spgr
          spgr_r{j}(:,:,i) = imresize(spgr(j).img(:,:,i), [res.y res.x]);
      end
      mtoffr(:,:,i)  = imresize(mtoff(:,:,i), [res.y res.x]);
      for j=1:n_mt
          mt_r{j}(:,:,i) = imresize(mt(j).img(:,:,i), [res.y res.x]);
      end
end

afi0=afi0r; afi1=afi1r;
for i=1:n_spgr
    spgr(i).img = spgr_r{i};
end
mtoff = mtoffr;
for i=1:n_mt
    mt(i).img = mt_r{i};
end

%% Clean and Close up

clear i j mt_r spgr_r afi0r afi1r mtoffr
save _recon mt spgr mtoff afi0 afi1 -append