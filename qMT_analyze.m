% QMT ANALYZE
%   Script to follow after qMT Processing script. Make sure mask and fv2 
%   exist prior to running.

%% Cartilage Segmentation: Post-fit

maskroi = zeros(size(mask));
figure
for i=1:length(z)
    imsc(mimsc(fv2(:,:,i,1)),[0 0.2])
    maskroi(:,:,i) = roipoly;
end

mask = mask.*maskroi;

%% Produce Values

f_mean = mean(nonzeros(fv2(:,:,:,1).*mask))
f_std = std(nonzeros(fv2(:,:,:,1).*mask))
k_mean = mean(nonzeros(fv2(:,:,:,2).*mask))
k_std = std(nonzeros(fv2(:,:,:,2).*mask))
t2b_mean = mean(nonzeros(fv2(:,:,:,3).*mask))
t2b_std = std(nonzeros(fv2(:,:,:,3).*mask))

save _values mask f_mean f_std k_mean k_std t2b_mean t2b_std