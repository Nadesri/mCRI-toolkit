% FUNCTION [fv rnrm] = qMT_model_fit_final_noPD_fixedT2b( data_MT, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT, mode, in16 (data_SPGR or PD), in17 (alpha_SPGR or R1), TR_SPGR, t_m_SPGR, alpha_ref, TR_ref)
%
% Fits MT & SPGR signal to a two-componant T1 & T2 model using the MT technique, or just fits MT if PD and R1 map are provided
%
%
%           data_MT          [Np x No]   - Normalized MT-Weighted Data by Reference Data ( Np : Number of Points, No : Number of Offest Frequencies for all  MT Data)
%           MT_Off_Freq      [No x 1 ]   - MT Data Offset Frequency Values in Hz
%           alpha_MT         [No x 1 ]   - MT Saturation Flip Angle, in Degrees
%           alpha_0          [No x 1 ]   - Excitation Flip Angle of MT Data, in Degrees
%           t_r              [No x 1 ]   - Delay Time After the Exitation Pulse in s (Vector of Time in case TR has been changed duo to SAR limitation)
%           t_m              [No x 1 ]   - Duration of an Off-Resonance RF Pulse in s (Vector)
%           t_s                          - Delay Time Before the Exitation Pulse in s
%           pulse_type                   - Type of MT pulse('sinc' or 'fermi')
%           line_shape                   - Shape of Saturation Line of Bound Pool ('Gaussian' or 'SuperLorentzian')
%           table_update                 - Update SuperLorentian/Lorentian Table ( 0 : No Update, 1 : Update )
%           B0               [Np x 1 ]   - Vector of B0 Error Deviation(-,+)
%           B1               [Np x 1 ]   - Vector of B1 Error Deviation(-,+)
%           B1_MT            [Np x 1 ]   - Vector of B1 Error Deviation(-,+) for MT Pulse
%           mode                         - Fitting Mode  0 : Fit Both MT and SPGR data for 5 Parameters [ PD R1 f k T2_B]
%                                                        1 : MT Fit Only,Using PD and R1 Map for 3 Parameters [ f k T2_B]
%                                                        2 : Fit Both MT and SPGR data for 6 Parameters [ PD R1 f k T2_B T2_F]
%                                                        3 : MT Fit Only,Using PD and R1 Map for 4 Parameters [ f k T2_B T2_F]
%           data_SPGR        [Np x Nf]   - Normalized SPGR Data by Reference Data( Np : Number of Points, Nf : Number of Flip Angles)
%           PD               [Np x 1 ]   - Vector of Proton Density ( Np : Number of Points)
%           alpha_SPGR       [Nf x 1 ]   - SPGR Flip Angles, in Degrees ( Nf : Number of Flip Angles)
%           R1               [Np x 1 ]   - Vector of Relaxation Rate(1/T1)( Np : Number of Points)
%           TR_SPGR          [Nf x 1 ]   - Repitation Time of SPGR Data, in s
%           t_m_SPGR         [1  x 1 ]   - Duration of an Off-Resonance RF Pulse in s for SPGR Data
%           alpha_ref        [1  x 1 ]   - Excitation Flip Angle of Reference Data, in Degrees
%           TR_ref           [1  x 1 ]   - Repitation Time of Reference Data, in s
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
% v6.1 16-Feb-2012 Output the fitted signal
% v6.2 01-Mar-2012 getting separate information of reference data


function [fv rnrm fvals] = qMT_model_fit_final_noPD_fixedT2b( data_MT, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT,  mode, in15, in16, TR_SPGR, t_m_SPGR, alpha_ref, TR_ref,T2b)

if nargin ~= 21
    error('Number of inputs is incorrect.');
end

% Calculate w_1rms Using Alexey's MATLAB Function (getB1eff_vasily)
w_1rms0   = zeros([numel(MT_Off_Freq) 1]);
for ii = 1:numel(alpha_MT)
    w_1rms0(ii) = getB1eff_vasily(alpha_MT(ii), t_m(ii), pulse_type);
end

% Convert flip angles to radians
alpha_0    = alpha_0    .* pi/180;
alpha_ref  = alpha_ref  .* pi/180;

% Apply B1 Correction to Flip Angles
alpha_0    = (1+B1) * alpha_0;
alpha_ref  = (1+B1) * alpha_ref;

if mode == 1 || mode == 3
    PD = in15;
    R1 = in16;
    
elseif mode == 0 || mode == 2
    data_SPGR   = in15;
    alpha_SPGR  = in16;
    
    % Convert flip angles to radians
    alpha_SPGR  = alpha_SPGR  .* pi/180;
    % Apply B1 Correction to Flip Angles
    alpha_SPGR  = (1+B1) * alpha_SPGR;
    
else
    error('Incorrect mode');
end

OPT = optimset('lsqnonlin');
OPT = optimset(OPT, 'Display', 'off');
OPT.DiffMinChange = 1e-16;
OPT.TolX          = 1e-16;
OPT.TolFun        = 1e-16;

% Loop Through Each Voxel in the Image
npts  = size(data_MT, 1);

% Preallocate Vectors
switch mode
    case 0
        fv    = zeros([npts 4]);
    case 1
        fv    = zeros([npts 3]);
    case 2
        fv    = zeros([npts 5]);
    case 3
        fv    = zeros([npts 4]);
    otherwise
        error('Incorrecet Mode');
end

% Preallocate Vectors
if exist('fvals','var')
    fvals = zeros([npts numel(alpha_SPGR)+numel(MT_Off_Freq)]);
end

rnrm  = zeros([npts 1]);
W_B   = zeros([numel(MT_Off_Freq) 1]);
W_F   = zeros([numel(MT_Off_Freq) 1]);

% w = waitbar(0, 'MT Fitting...');
fprintf('MT Fit...');

if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
    g_B = superLorentian_fastPM([],[],1,table_update,t_m(1));
end

for ii = 1:npts               % Loop over voxels
    % Grab MT Data and Offset Freq for Current Voxel
    progressbar(ii/npts);
    
    vox_data_MT = data_MT(ii,:)';
    
    if sum(vox_data_MT) ~= 0
        
        % Grab SPGR Flip Angle and MT Offset Frequency for Current Voxel
        vox_MT_Off_Freq = MT_Off_Freq'- B0(ii);
        
        w_1rms = w_1rms0*(1 + B1_MT(ii));
        
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
        vox_alpha_0   = alpha_0(ii,:);
        vox_alpha_ref = alpha_ref(ii);
        
        if mode == 0 || mode == 2
            % Grab SPGR Data for Current Voxel
            vox_data_SPGR = data_SPGR(ii,:)';
            
            % Grab SPGR and MT Flip Angles for Current Voxel
            vox_alpha_SPGR = alpha_SPGR(ii,:)';
            
            if mode == 0
                [xx residual] = lsqnonlin(@qMT_SPGR_Signal_Wrapper, guess(1:4,1), guess(1:4,2), guess(1:4,3), OPT);
            elseif mode == 2
                [xx residual] = lsqnonlin(@qMT_SPGR_Signal_Wrapper, guess(:,1), guess(:,2), guess(:,3), OPT);
            end
            
            rnrm(ii) = sqrt(residual)/(sum(vox_data_SPGR) + sum(vox_data_MT));
            
            if exist('fvals','var')
                fvals(ii,:)=qMT_SPGR_Signal_Wrapper(xx);
            end
            

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
            
            if exist('fvals','var')
                fvals(ii,:)=qMT_Signal_Wrapper(xx);
            end
           

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

% Wrapper for qMT_residuals
    function residual = qMT_SPGR_Signal_Wrapper(x)
        
        % Fix input T2b value
        x(4,1) = T2b;
        
        % Calculate Saturation Rate for Free and Bound Pool
        if mode == 0
            R1_B = 1;
            R1_F = x(1,1) - x(3,1).*(R1_B - x(1,1))./(R1_B - x(1,1) + x(3,1).*(1-x(2,1))./x(2,1)) ;
            W_F  = R1_F * W_F0;
            
        elseif mode == 2
            g_F = x(5,1)./(pi*(1 + (2*pi* vox_MT_Off_Freq.*x(5,1)).^2));
            W_F = W .*g_F;
        end
        
        % Calculate the Line Shape for Bound Pool
        if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
            g_B = superLorentian_fastPM(vox_MT_Off_Freq,x(4,1));
            W_B = W .*g_B;
            
        elseif strcmp(line_shape, 'Gaussian')
            E = pi * (w_1rms).^2 .*exp(-(2*pi*x(4,1)*(vox_MT_Off_Freq)).^2/2);
            W_B = x(4,1)./(sqrt(2*pi)).*E;
            
        else
            error('Uknown line shape. Available choices: SuperLorentzian (Neural Tissues), Gaussian (Agar)');
        end
        
        S_MT     = qMT_Signal_final(x(2,1), x(3,1), vox_alpha_0'  , t_m'  , t_s, t_r'   , 1, x(1,1), W_B, W_F  );
        S_MT_Off = qMT_Signal_final(x(2,1), x(3,1), vox_alpha_ref , 0     , 0  , TR_ref , 1, x(1,1), 0  , 0  );
        
        S_MT_norm = S_MT/S_MT_Off ;
        
        S_SPGR = zeros(size(vox_alpha_SPGR'));
        
        for ll = 1:numel(vox_alpha_SPGR)
            S_SPGR(ll) = qMT_Signal_final(x(2,1), x(3,1), vox_alpha_SPGR(ll), t_m_SPGR, t_s, TR_SPGR(ll)-(t_m_SPGR+t_s), 1, x(1,1), 0, 0);
        end
        
        S_SPGR_norm = S_SPGR/S_MT_Off ;
        
        % Compute residual based
        residual=[S_SPGR_norm-vox_data_SPGR' S_MT_norm-vox_data_MT'];
        
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
%     function ig = InitialGuessHuman()
%         %         [ mean    min       max     ]
% %         % PD
% %         ig(1,:) = [0.005    0         100000  ];
%         % R1
%         ig(1,:) = [1.000    0.01       10.00   ];
%         % f 
%         ig(2,:) = [0.1000    0.001     0.500   ];
%         % k 
%         ig(3,:) = [1.000    0.0100    100.00  ];
%         % T2_B
%         ig(4,:) = [7.0e-6  1.0e-6    50.0e-6 ];
%         % T2_F
%         ig(5,:) = [100e-3    10e-3     200e-3 ];
%     end

    function ig = InitialGuessHuman()
        %         [ mean    min       max     ]
%         % PD
%         ig(1,:) = [0.005    0         100000  ];
        % R1
        ig(1,:) = [0.9063    0.01       10.00   ];
        % f 
        ig(2,:) = [0.1414    0.001     0.500   ];
        % k 
        ig(3,:) = [1.9828    0.0100    100.00  ];
        % T2_B
        ig(4,:) = [T2b  1.0e-6    50.0e-6 ];
        % T2_F
        ig(5,:) = [100e-3    10e-3     200e-3 ];
    end

    function ig = InitialGuessPhantom()
        %         [ mean    min        max     ]
%         % PD
%         ig(1,:) = [0.005    0          100000  ];
        % R1
        ig(1,:) = [1.000    0.1        20.00   ];
        % f
        ig(2,:) = [0.00001  0.000001   0.5     ];
        % k
        ig(3,:) = [2.000    0.000100   100.00  ];
        % T2_B
        ig(4,:) = [T2b   1.0e-6     50.0e-6 ];
        % T2_F
        ig(5,:) = [50e-3    10e-3      4000e-3 ];
    end

% /END of MAIN
end
