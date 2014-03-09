% SCRIPT: qMT_Spectra
% Calculate z-spectrum
%
% Depending on the cell, will evaluate z-spectra through various sampling
% methods. Requires tt, which is a comparison variable obtained through
% qMT_comparisons.m
%
% Nade Sritanyaratana
% University of Wisconsin, Madison
% Created: June 6, 2011

%% Z-spectrum along a pt

% Clean and Start up
clear, close all, load _comparisons

figure, imsc(mimsc(tt(:,:,1)));
title(['MT measurement at Offset: ' int2str(qMTparams(1,1)) ' and alpha: ' int2str(qMTparams(1,2))])

% Choose a point
[X Y] = ginput(1); X = round(X); Y=round(Y);
disp(['Choosing the point at Y=' int2str(Y) ' and X=' int2str(X)]);

% look for MT powers.
alpha_MT_unique=[];
for i=1:size(qMTparams,1)
    if sum(qMTparams(i,2)==alpha_MT_unique)==0
        alpha_MT_unique = [alpha_MT_unique qMTparams(i,2)];
    end
end

zspectrum = zeros(1,length(alpha_MT_unique));
offset    = zeros(1,length(alpha_MT_unique));
for i=1:length(alpha_MT_unique)
    k(i)=1;
    for j=1:size(tt,3)
        if qMTparams(j,2)==alpha_MT_unique(i)
            zspectrum(k(i),i) = tt(Y,X,j);
            offset(k(i),i)    = qMTparams(j,1);
            k(i)=k(i)+1;
        end
    end
end

k = k-1; % #significant rows in zspectrum/offset for each power

for i=1:length(alpha_MT_unique)
    % prep the plotting vectors
    xvect = offset(:,i);    xvect = xvect(1:k(i));
    yvect = zspectrum(:,i); yvect = yvect(1:k(i)); % k removes trailing zeros
    plot(xvect,yvect,'o'), hold all
    legstr{i} = int2str(alpha_MT_unique(i)); %string for legend
end
title(['Point based z-spectrum for X=' int2str(X) ' and Y=' int2str(Y)])
legend(legstr,'Location','Best')
xlabel('Offset Frequency')
ylabel('MT(on)/MT(off)')
clear xvect yvect legstr i j

%% Z-spectrum along an ROI

% Clean and Start up
clear, close all, load _comparisons

imsample=tt(:,:,1);
figure, imsc(mimsc(imsample));
title(['MT measurement at Offset: ' int2str(qMTparams(1,1)) ' and alpha: ' int2str(qMTparams(1,2))])

% Choose an ROI
mask = roipoly; close,
figure, imsc(mimsc(imsample+mask*max(imsample(:))))
title(['MT measurement at Offset: ' int2str(qMTparams(1,1)) ' and alpha: ' int2str(qMTparams(1,2))])
clear imsample

% look for MT powers.
alpha_MT_unique=[];
for i=1:size(qMTparams,1)
    if sum(qMTparams(i,2)==alpha_MT_unique)==0
        alpha_MT_unique = [alpha_MT_unique qMTparams(i,2)];
    end
end

zspectrum = zeros(1,length(alpha_MT_unique));
offset    = zeros(1,length(alpha_MT_unique));
for i=1:length(alpha_MT_unique)
    k(i)=1;
    for j=1:size(tt,3)
        if qMTparams(j,2)==alpha_MT_unique(i)
            tmp = tt(:,:,j).*mask;
            zspectrum(k(i),i) = mean(nonzeros(tmp(:)));
            offset(k(i),i)    = qMTparams(j,1);
            k(i)=k(i)+1;
        end
    end
end

k = k-1; % #significant rows in zspectrum/offset for each power
figure,
for i=1:length(alpha_MT_unique)
    % prep the plotting vectors
    xvect = offset(:,i);    xvect = xvect(1:k(i));
    yvect = zspectrum(:,i); yvect = yvect(1:k(i)); % k removes trailing zeros
    plot(xvect,yvect,'o'), hold all
    legstr{i} = int2str(alpha_MT_unique(i)); %string for legend
end
title('ROI based z-spectrum')
legend(legstr,'Location','Best')
xlabel('Offset Frequency')
ylabel('MT(on)/MT(off)')
clear xvect yvect legstr i j

%% Z-spectrum along a prepared mask

% Clean and Start up
clear, close all, load _comparisons
zthis = z;

% CHOOSE A SAVED MASK HERE
load _data mask
load _inits z
zthat = z;
mask = mask(:,:,zthat==zthis);

imsample = tt(:,:,1);
figure, imsc(mimsc(imsample+mask*max(imsample(:))));
title(['MT measurement at Offset: ' int2str(qMTparams(1,1)) ' and alpha: ' int2str(qMTparams(1,2))])
clear imsample

% look for MT powers.
alpha_MT_unique=[];
for i=1:size(qMTparams,1)
    if sum(qMTparams(i,2)==alpha_MT_unique)==0
        alpha_MT_unique = [alpha_MT_unique qMTparams(i,2)];
    end
end

zspectrum = zeros(1,length(alpha_MT_unique));
offset    = zeros(1,length(alpha_MT_unique));
for i=1:length(alpha_MT_unique)
    k(i)=1;
    for j=1:size(tt,3)
        if qMTparams(j,2)==alpha_MT_unique(i)
            tmp = tt(:,:,j).*mask;
            zspectrum(k(i),i) = mean(nonzeros(tmp(:)));
            offset(k(i),i)    = qMTparams(j,1);
            k(i)=k(i)+1;
        end
    end
end

k = k-1; % #significant rows in zspectrum/offset for each power
figure,
for i=1:length(alpha_MT_unique)
    % prep the plotting vectors
    xvect = offset(:,i);    xvect = xvect(1:k(i));
    yvect = zspectrum(:,i); yvect = yvect(1:k(i)); % k removes trailing zeros
    plot(xvect,yvect,'o'), hold all
    legstr{i} = int2str(alpha_MT_unique(i)); %string for legend
end
title('Mask based z-spectrum')
legend(legstr,'Location','Best')
xlabel('Offset Frequency')
ylabel('MT(on)/MT(off)')
clear xvect yvect legstr i j