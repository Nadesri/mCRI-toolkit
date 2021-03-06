% FUNCTION [fv rnrm] = qMT_model_fit_final( data_MT, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT, mode, in16 (data_SPGR or PD), in17 (alpha_SPGR or R1), TR_SPGR, t_m_SPGR)
%
% Fits MT & SPGR signal to a two-componant T1 & T2 model using the MT technique, or just fits MT if PD and R1 map are provided
%
%
%           data_MT          [Np x No]   - MT-Weighted Data ( Np : Number of Points, No : Number of Offest Frequencies for all  MT Data)
%           MT_Off_Freq      [No x 1 ]   - MT Data Offset Frequency Values in Hz
%           alpha_MT         [No x 1 ]   - MT Saturation Flip Angle, in Degrees
%           alpha_0          [No x 1 ]   - Excitation Flip Angle of MT Data, in Degrees
%           t_r              [No x 1 ]   - Delay Time After the Exitation Pulse in s (Vector of Time in case TR has been changed duo to SAR limitation)
%           t_m              [No x 1 ]   - Duration of an Off-Resonance RF Pulse in s (Vector)
%           t_s                          - Delay Time Before the Exitation Pulse in s
%           pulse_type                   - Type of MT pulse('sinc' or 'fermi')
%           line_shape                   - Shape of Saturation Line of Bound Pool ('Gaussian' or 'SuperLorentzian')
%           table_update                 - Update SuperLorentian/Lorentian Table ( 0 : No Update, 1 : Update )
%           B0               [Np x 1 ]   - Vector(or Matrix [Np x No] ) of B0 Error Deviation(-,+)
%           B1               [Np x 1 ]   - Vector(or Matrix [Np x (No+Nf)]) of B1 Error Deviation(-,+)
%           B1_MT            [Np x 1 ]   - Vector(or Matrix [Np x No]) of B1 Error Deviation(-,+) for MT Pulse
%           mode                          - Fitting Mode  0 : Fit Both MT and SPGR data for 5 Parameters [ PD R1 f k T2_B]
%                                                         1 : MT Fit Only,Using PD and R1 Map for 3 Parameters [ f k T2_B]
%                                                         2 : Fit Both MT and SPGR data for 6 Parameters [ PD R1 f k T2_B T2_F]
%                                                         3 : MT Fit Only,Using PD and R1 Map for 4 Parameters [ f k T2_B T2_F]
%           data_SPGR        [Np x Nf]    - SPGR Data ( Np : Number of Points, Nf : Number of Flip Angles)
%           PD               [Np x 1]     - Vector of Proton Density ( Np : Number of Points)
%           alpha_SPGR       [Nf x 1]     - SPGR Flip Angles, in Degrees ( Nf : Number of Flip Angles)
%           R1               [Np x 1]     - Vector of Relaxation Rate(1/T1)( Np : Number of Points)
%
%
% Outputs:
%           fv               [Np x N]    - [PD R1 f k T2_B] or [f k T2_B] or [PD R1 f k T2_B T2_F] or [f k T2_B T2_F] N = 5,3,6,4 for different mode
%           rnrm             [Np x 1]    - residuals
%
% Pouria Mossahebi
% University of Wisconsin
% v1.0 07-Sep-2010
% v1.1 10-Sep-2010 add mode 2 and 3
% v1.2 15-Sep-2010 getting t_r as a vector
% v2.1 28-Sep-2010 adding B0 correcetion
% v3.0 15-Jan-2011 supply different MT data with different flip angles
% v3.1 01-Feb-2011 supply different B1 map for MT pulse
% v4.0 27-Apr-2011 supply the option to update SuperLorentian/Lorentian table
% v5.0 28-Apr-2011 supply all MT data as one input
% v6.0 27-Oct-2011 getting t_m as a vector
% v7.0 15-Mar-2012 getting B0 and B1 as matrices for using correct values of them due to motion

function [fv rnrm] = qMT_model_fit_final_smartig( data_MT, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT,  mode, in15, in16, TR_SPGR, t_m_SPGR)

disp('WARNING: Using new initial guess matrix as per Pouria''s suggestions. TODO: Compare/Contrast fit with previous invivo data!!!!! (NS 10/06/13)')

if nargin ~= 18
    error('Number of inputs is incorrect.');
end

% Calculate w_1rms Using Alexey's MATLAB Function (getB1eff_vasily)
w_1rms0   = zeros([numel(MT_Off_Freq) 1]);
for ii = 1:numel(alpha_MT)
    w_1rms0(ii) = getB1eff_vasily(alpha_MT(ii), t_m(ii), pulse_type);
end

if size(B1_MT,2) ~= 1 && size(B1_MT,2) ~=length(alpha_0)
        error('Size of B1_MT is Incorrect.');
end

if size(B0,2) ~= 1 && size(B0,2) ~=length(alpha_0)
        error('Size of B0 is Incorrect.');
end

if mode == 1 || mode == 3
    PD = in15;
    R1 = in16;
%     max_val     = max(data_MT(:));
%     data_MT     = data_MT/max_val;
    
    if size(B1,2) ~= 1 && size(B1,2) ~=length(alpha_0)
        error('Size of B1 is Incorrect.');
    end
    
elseif mode == 0 || mode == 2
    data_SPGR   = in15;
    alpha_SPGR  = in16;
    max_val     = max(max(data_MT(:)),max(data_SPGR(:)));
    data_MT     = data_MT/max_val;
    data_SPGR   = data_SPGR/max_val;
    
    % Convert flip angles to radians
    alpha_SPGR  = alpha_SPGR  .* pi/180;
    
    if size(B1,2) ~= 1 && size(B1,2) ~=length(alpha_0)+length(alpha_SPGR)
        error('Size of B1 is Incorrect.');
    end
else
    error('Incorrect mode');
end

npts  = size(data_MT, 1);
    

% Apply B1 Correction to Flip Angles
alpha_0_tmp = alpha_0 .* pi/180;
alpha_0     = zeros([npts numel(alpha_0_tmp)]);

if size(B1,2) == 1
    alpha_0  = (1+B1) * alpha_0_tmp;
else
    for ii = 1:numel(alpha_0_tmp)
        alpha_0(:,ii) = (1+B1(:,ii))* alpha_0_tmp(ii);
    end
end

if mode == 0 || mode == 2
    alpha_tmp   = alpha_SPGR;
    alpha_SPGR  = zeros([npts numel(alpha_tmp)]);
    
    if size(B1,2) == 1
        alpha_SPGR  = (1+B1) * alpha_tmp;
    else
        for ii = 1 : numel(alpha_tmp)
            alpha_SPGR(:,ii) = (1+B1(:,numel(alpha_0_tmp)+ii))* alpha_tmp(ii);
        end
    end
end
    

OPT = optimset('lsqnonlin');
OPT = optimset(OPT, 'Display', 'off');
OPT.DiffMinChange = 1e-16;
OPT.TolX          = 1e-16;
OPT.TolFun        = 1e-16;


% Preallocate Vectors
switch mode
    case 0
        fv    = zeros([npts 5]);
    case 1
        fv    = zeros([npts 3]);
    case 2
        fv    = zeros([npts 6]);
    case 3
        fv    = zeros([npts 4]);
    otherwise
        error('Incorrect Mode');
end

% Preallocate Vectors
rnrm  = zeros([npts 1]);
W_B   = zeros([numel(MT_Off_Freq) 1]);
W_F   = zeros([numel(MT_Off_Freq) 1]);

% w = waitbar(0, 'MT Fitting...');
fprintf('MT Fit...');

if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
    g_B = superLorentian_fastPM([],[],1,table_update,t_m(1));
end

%% Find IG using the fitted mean of data

vox_data_MT = mean(data_MT,1)';
B0mean = mean(B0,1);
B1mean = mean(B1,1);
B1_MTmean = mean(B1_MT);

alpha_0mean = mean(alpha_0,1);

if (mode==0)||(mode==2)
    alpha_SPGRmean = mean(alpha_SPGR,1);
end

if sum(vox_data_MT) ~= 0

    vox_MT_Off_Freq = zeros(size(MT_Off_Freq'));

    % Grab MT Offset Frequency for Current Voxel

    if size(B0,2) == 1
        vox_MT_Off_Freq = MT_Off_Freq'- B0mean;
    else
        for kk = 1:numel(MT_Off_Freq)
            vox_MT_Off_Freq(kk) = MT_Off_Freq(kk)- B0mean(kk);
        end
    end

    if size(B1_MT,2) == 1
        w_1rms = w_1rms0*(1 + B1_MTmean);
    else
        for kk = 1:numel(MT_Off_Freq)
            w_1rms(kk,1) = w_1rms0(kk)*(1 + B1_MTmean(kk));
        end
    end

    if mode == 0 || mode == 1
        W_F0   = ((w_1rms./(2*pi* vox_MT_Off_Freq)).^2)/0.022;
    end

    W   = pi * (w_1rms).^2;

    % Reset guess to initial guess
    if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
        guess = InitialGuessHuman();

    elseif strcmp(line_shape, 'Gaussian')
        guess = InitialGuessPhantom();
    end

    % Grab MT Flip Angle for Current Voxel
    vox_alpha_0    = alpha_0mean(:)';

    if mode == 0 || mode == 2
        % Grab SPGR Data for Current Voxel
        vox_data_SPGR = mean(data_SPGR,1)';

        % Grab SPGR and MT Flip Angles for Current Voxel
        vox_alpha_SPGR = alpha_SPGRmean(:)';

        if mode == 0
            [xx] = lsqnonlin(@qMT_SPGR_Signal_Wrapper, guess(1:5,1), guess(1:5,2), guess(1:5,3), OPT);
        elseif mode == 2
            [xx] = lsqnonlin(@qMT_SPGR_Signal_Wrapper, guess(:,1), guess(:,2), guess(:,3), OPT);
        end

        xx(1)    = xx(1)*max_val;

    elseif mode == 1 || mode == 3

        % Grab PD and R1 Map for Current Voxel
        vox_PD = mean(PD)';
        vox_R1 = mean(R1)';

        if mode == 1
            [xx] = lsqnonlin(@qMT_Signal_Wrapper, guess(3:5,1),guess(3:5,2),guess(3:5,3),OPT);
        elseif mode == 3
            [xx] = lsqnonlin(@qMT_Signal_Wrapper, guess(3:6,1),guess(3:6,2),guess(3:6,3),OPT);
        end

    else
        error('Unknown Mode. Should be either 0-2(MT-SPGR Fit) or 1-3(MT Fit)');
    end

    smart_ig = xx;

end

if mode==0
    smart_ig = [smart_ig;guess(1,end)];
elseif mode==1
    smart_ig = [vox_PD;vox_R1;smart_ig;guess(1,end)];
elseif mode==2
    % do nothing - smart_ig is the right size!
elseif mode==3
    smart_ig = [vox_PD;vox_R1;smart_ig];
end 
smart_ig = [smart_ig, guess(:,2), guess(:,3)];


%% Main Body

% Loop Through Each Voxel in the Image
for ii = 1:npts               % Loop over voxels
    % Grab MT Data and Offset Freq for Current Voxel
    progressbar(ii/npts);
    
    vox_data_MT = data_MT(ii,:)';
    
    if sum(vox_data_MT) ~= 0
        
        vox_MT_Off_Freq = zeros(size(MT_Off_Freq'));
        
        % Grab MT Offset Frequency for Current Voxel
        
        if size(B0,2) == 1
            vox_MT_Off_Freq = MT_Off_Freq'- B0(ii);
        else
            for kk = 1:numel(MT_Off_Freq)
                vox_MT_Off_Freq(kk) = MT_Off_Freq(kk)- B0(ii,kk);
            end
        end
        
        if size(B1_MT,2) == 1
            w_1rms = w_1rms0*(1 + B1_MT(ii));
        else
            for kk = 1:numel(MT_Off_Freq)
                w_1rms(kk,1) = w_1rms0(kk)*(1 + B1_MT(ii,kk));
            end
        end
 
        if mode == 0 || mode == 1
            W_F0   = ((w_1rms./(2*pi* vox_MT_Off_Freq)).^2)/0.022;
        end
        
        W   = pi * (w_1rms).^2;
        
        % Reset guess to initial guess
        if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
            guess = smart_ig;
            
        elseif strcmp(line_shape, 'Gaussian')
            guess = smart_ig;
        end
        
        % Grab MT Flip Angle for Current Voxel
        vox_alpha_0    = alpha_0(ii,:);
        
        if mode == 0 || mode == 2
            % Grab SPGR Data for Current Voxel
            vox_data_SPGR = data_SPGR(ii,:)';
            
            % Grab SPGR and MT Flip Angles for Current Voxel
            vox_alpha_SPGR = alpha_SPGR(ii,:)';
            
            if mode == 0
                [xx residual] = lsqnonlin(@qMT_SPGR_Signal_Wrapper, guess(1:5,1), guess(1:5,2), guess(1:5,3), OPT);
            elseif mode == 2
                [xx residual] = lsqnonlin(@qMT_SPGR_Signal_Wrapper, guess(:,1), guess(:,2), guess(:,3), OPT);
            end
            
            rnrm(ii) = sqrt(residual)/(sum(vox_data_SPGR) + sum(vox_data_MT));
            
            xx(1)    = xx(1)*max_val;
            
        elseif mode == 1 || mode == 3
            
            % Grab PD and R1 Map for Current Voxel
            vox_PD = PD(ii)';
            vox_R1 = R1(ii)';
            
            if mode == 1
                [xx residual] = lsqnonlin(@qMT_Signal_Wrapper, guess(3:5,1),guess(3:5,2),guess(3:5,3),OPT);
            elseif mode == 3
                [xx residual] = lsqnonlin(@qMT_Signal_Wrapper, guess(3:6,1),guess(3:6,2),guess(3:6,3),OPT);
            end
            
            rnrm(ii) = sqrt(residual)/(sum(vox_data_MT));
            
        else
            error('Uknown Mode. Should be either 0-2(MT-SPGR Fit) or 1-3(MT Fit)');
        end
        
        fv(ii,:) = xx;
        
    end
    
end     % End voxels

if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
    g_B = superLorentian_fastPM([],[],-1);
end

disp('done!');

%% Helper functions and initial IG

% Wrapper for qMT_residuals
    function residual = qMT_SPGR_Signal_Wrapper(x)
        
        % Calculate Saturation Rate for Free and Bound Pool
        if mode == 0
            R1_B = 1;
            R1_F = x(2,1) - x(4,1).*(R1_B - x(2,1))./(R1_B - x(2,1) + x(4,1).*(1-x(3,1))./x(3,1)) ;
            W_F  = R1_F * W_F0;
            
        elseif mode == 2
            g_F = x(6,1)./(pi*(1 + (2*pi* vox_MT_Off_Freq.*x(6,1)).^2));
            W_F = W .*g_F;            
        end
        
        % Calculate the Line Shape for Bound Pool
        if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
%             T2b_init_val = x(5,1)
            g_B = superLorentian_fastPM(vox_MT_Off_Freq,x(5,1));
            W_B = W .*g_B;
            
        elseif strcmp(line_shape, 'Gaussian')
            E = pi * (w_1rms).^2 .*exp(-(2*pi*x(5,1)*(vox_MT_Off_Freq)).^2/2);
            W_B = x(5,1)./(sqrt(2*pi)).*E;
            
        else
            error('Uknown line shape. Available choices: SuperLorentzian (Neural Tissues), Gaussian (Agar)');
        end
        
        S_MT = qMT_Signal_final(x(3,1), x(4,1), vox_alpha_0', t_m', t_s, t_r', x(1,1), x(2,1), W_B, W_F);
        
        S_SPGR = zeros(size(vox_alpha_SPGR'));
        
        for ll = 1:numel(vox_alpha_SPGR)
            S_SPGR(ll) = qMT_Signal_final(x(3,1), x(4,1), vox_alpha_SPGR(ll), t_m_SPGR, t_s, TR_SPGR(ll)-(t_m_SPGR+t_s), x(1,1), x(2,1), 0, 0);            
        end
        
        % Compute residual based
        residual=[S_SPGR-vox_data_SPGR' S_MT-vox_data_MT'];
        
        residual(isnan(residual)) = 1e99;
        residual(isinf(residual)) = 1e99;
        
    end

    function residual = qMT_Signal_Wrapper(x)
        
        % Calculate Saturation Rate for Free and Bound Pool
        if mode == 1
            R1_B = 1;
            R1_F = vox_R1 - x(2,1).*(R1_B - vox_R1)./(R1_B - vox_R1 + x(2,1).*(1-x(1,1))./x(1,1)) ;
            W_F  = R1_F * W_F0;
            
        elseif mode == 3
            g_F = x(4,1)./(pi*(1 + (2*pi* vox_MT_Off_Freq.*x(4,1)).^2));
            W_F = W .*g_F;
        end
        % Calculate the Line Shape for Bound Pool
        if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
            
            g_B = superLorentian_fastPM(vox_MT_Off_Freq,x(3,1));
            W_B = W .*g_B;
            
        elseif strcmp(line_shape, 'Gaussian')
            
            E = pi * (w_1rms).^2 .*exp(-(2*pi*x(3,1)*(vox_MT_Off_Freq)).^2/2);
            W_B = x(3,1)./(sqrt(2*pi)).*E;
            
        else
            error('Uknown line shape. Available choices: SuperLorentzian (Neural Tissues), Gaussian (Agar)');
        end
        
        S_MT = qMT_Signal_final(x(1,1), x(2,1), vox_alpha_0', t_m', t_s, t_r', vox_PD, vox_R1, W_B, W_F);
        
        % Compute residual
        residual = S_MT-vox_data_MT';
        
        residual(isnan(residual)) = 1e99;
        residual(isinf(residual)) = 1e99;
        
    end

% ------ Helper Functions Below ------
% Initial Guess for 3.0T Magnet
    function ig = InitialGuessHuman()
        %         [ mean    min       max     ]
%         ig(1,:) = [0.005    0         100000  ]; % PD
%         ig(2,:) = [1.000    0.2       10.00   ]; % R1
%         ig(3,:) = [0.100    0.015     0.500   ]; % f
%         ig(4,:) = [2.000    0.0100    100.00  ]; % k 
%         ig(5,:) = [11e-6    1.1e-6    49.9e-6 ]; % T2_B
%         ig(6,:) = [50e-3    10e-3     4000e-3 ]; % T2_F
        
        % Per Pouria's suggestions
        ig(1,:) = [0.005    0         100000  ]; % PD
        ig(2,:) = [0.300    0.2       10.00   ]; % R1
        ig(3,:) = [0.016    0.015     0.500   ]; % f
        ig(4,:) = [0.020    0.0100    50.00  ];  % k
        ig(5,:) = [1.2e-6   1.1e-6    49.9e-6 ]; % T2b
        ig(6,:) = [11e-3    10e-3     4000e-3 ]; % T2f
    end

    function ig = InitialGuessPhantom()
        %         [ mean    min        max     ]
        % PD
        ig(1,:) = [0.005    0          100000  ];
        % R1
        ig(2,:) = [1.000    0.1        20.00   ];
        % f
        ig(3,:) = [0.00001  0.000001   0.5     ];
        % k
        ig(4,:) = [2.000    0.000100   100.00  ];
        % T2_B
        ig(5,:) = [24e-6    1.1e-6     49.9e-6 ];
        % T2_F
        ig(6,:) = [50e-3    10e-3      4000e-3 ];
    end

% /END of MAIN
end
