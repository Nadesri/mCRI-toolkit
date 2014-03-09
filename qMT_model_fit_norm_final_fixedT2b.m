% FUNCTION [fv rnrm] = qMT_model_fit_norm_final_fixedT2b( data_MT_norm, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT, mode, R1, alpha_ref, TR_ref)
%
% Fits MT signal to a two-componant T1 & T2 model using the MT technique 
%
%
%           data_MT_norm    [Np x No]   - MT-Weighted Data or MT Data Normalized to Reference Data ( Np : Number of Points, No : Number of Offest Frequencies for MT Data)
%           MT_Off_Freq     [No x 1 ]   - Data MT Offset Frequency Values in Hz
%           alpha_MT        [No x 1 ]   - MT Saturation Flip Angle, in Degrees
%           alpha_0         [No x 1 ]   - Excitation Flip Angle of First  MT Data, in Degrees
%           t_r             [No x 1 ]   - Delay Time After the Exitation Pulse in s (Vector of Time in case TR has been changed due to SAR limitation)
%           t_m             [No x 1 ]   - Duration of an Off-Resonance RF Pulse in s (Vector)
%           t_s                         - Delay Time Before the Exitation Pulse in s
%           pulse_type                  - Type of MT pulse('sinc' or 'fermi')
%           line_shape                  - Shape of Saturation Line of Bound Pool ('Gaussian' or 'SuperLorentzian')
%           table_update                - Update SuperLorentian/Lorentian Table ( 0 : No Update, 1 : Update )
%           B0              [Np x 1 ]   - Vector of B0 Error Deviation(-,+)
%           B1              [Np x 1 ]   - Vector of B1 Error Deviation(-,+)
%           B1_MT           [Np x 1 ]   - Vector of B1 Error Deviation(-,+) of MT Pulse
%           mode                        - 0 : Normalized MT Data Fit, for 3 Parameters [f k T2_B] 
%                                         1 : Normalized MT Data Fit, for 4 Parameters [ f k T2_B T2_F]
%           R1              [Np x 1]    - Vector of Relaxation Rate(1/T1)( Np : Number of Points)
%           alpha_ref       [1  x 1 ]   - Excitation Flip Angle of Reference Data, in Degrees
%           TR_ref          [1  x 1 ]   - Repitation Time of Reference Data, in s


% Outputs:
%           fv              [Np x N]    - [f k T2_B] or [f k T2_B  T2_F]
%           rnrm            [Np x 1]    - residuals
%
% Pouria Mossahebi
% University of Wisconsin
% v1.0 28-Oct-2011
% v2.0 01-Mar-2012 getting separate information of reference data



function [fv rnrm] = qMT_model_fit_norm_final_fixedT2b( data_MT_norm, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT, mode, R1, alpha_ref, TR_ref,T2b)

if nargin ~= 18
    error('Mode 0 and 1 require 18 inputs.');
end

    
% Calculate w_1rms Using Alexey's MATLAB Function (getB1eff_vasily)
w_1rms0   = zeros([numel(MT_Off_Freq) 1]);
for ii = 1:numel(alpha_MT)
    w_1rms0(ii) = getB1eff_vasily(alpha_MT(ii), t_m(ii), pulse_type);
end

% Convert flip angles to radians
alpha_0   = alpha_0   .* pi/180;
alpha_ref = alpha_ref .* pi/180;

% Apply B1 Correction to Flip Angles
alpha_0   = (1+B1) * alpha_0;
alpha_ref = (1+B1) * alpha_ref;

OPT = optimset('lsqnonlin');
OPT = optimset(OPT, 'Display', 'off');
OPT.DiffMinChange = 1e-16;
OPT.TolX          = 1e-16;
OPT.TolFun        = 1e-16;

% Loop Through Each Voxel in the Image
npts  = size(data_MT_norm, 1);

switch mode
    case 0 
        fv    = zeros([npts 3]);
    case 1 
        fv    = zeros([npts 4]);
    otherwise
        error('Incorrecet Mode');
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

for ii = 1:npts               % Loop over voxels
    % Grab MT Data and Offset Freq for Current Voxel
    progressbar(ii/npts);
    
    vox_data_MT_norm = data_MT_norm(ii,:)';

    if sum(vox_data_MT_norm) ~= 0
        
        % Grab SPGR Flip Angle and MT Offset Frequency for Current Voxel
        vox_MT_Off_Freq = MT_Off_Freq'- B0(ii);
        
        w_1rms = w_1rms0*(1 + B1_MT(ii));
      
        if mode == 0 
            W_F0   = 0*((w_1rms./(2*pi* vox_MT_Off_Freq)).^2)/0.022;
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
        
        switch mode
            case 0
                vox_R1 = R1(ii)';
                [xx residual] = lsqnonlin(@qMT_Signal_Normalized_Wrapper, guess(1:3,1),guess(1:3,2),guess(1:3,3),OPT);
            case 1
                vox_R1 = R1(ii)';
                [xx residual] = lsqnonlin(@qMT_Signal_Normalized_Wrapper, guess(:,1),guess(:,2),guess(:,3),OPT);
            
            otherwise
                error('Incorrecet Mode');
        end
        rnrm(ii) = sqrt(residual)/(sum(vox_data_MT_norm));
        fv(ii,:) = xx;
        
    end
    
end     % End voxels

if ~exist('line_shape','var') || isempty(line_shape) ||  strcmp(line_shape, 'SuperLorentzian')
    g_B = superLorentian_fastPM([],[],-1);
end

disp('done!');

% Wrapper for qMT_Normalized residuals
    function residual = qMT_Signal_Normalized_Wrapper(x)
        
        % Fix input T2b value
        x(3,1) = T2b;

        % Calculate Saturation Rate for Free and Bound Pool
        if mode == 0 
            R1_B = 1;
            R1_F = vox_R1 - x(2,1).*(R1_B - vox_R1)./(R1_B - vox_R1 + x(2,1).*(1-x(1,1))./x(1,1)) ;
            W_F  = R1_F * W_F0;
        elseif mode == 1
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
        
        S_MT     = qMT_Signal_final(x(1,1), x(2,1), vox_alpha_0'  , t_m'  , t_s, t_r'  , 1, vox_R1, W_B, W_F  );
        S_MT_Off = qMT_Signal_final(x(1,1), x(2,1), vox_alpha_ref , 0     , 0  , TR_ref, 1, vox_R1, 0  , 0  );
%         S_MT_Off  = sin(vox_alpha_ref).*(1-exp(-TR_ref*vox_R1))./(1-cos(vox_alpha_ref).*exp(-TR_ref*vox_R1));
        
        
        S_MT_norm = S_MT/S_MT_Off ;
        
        residual = S_MT_norm - vox_data_MT_norm';

        
        residual(isnan(residual)) = 1e99;
        residual(isinf(residual)) = 1e99;
            
    end

% ------ Helper Functions Below ------
% Initial Guess for 3.0T Magnet
%     function ig = InitialGuessHuman()
%         %         [ mean    min       max     ]
%         % f
%         ig(1,:) = [0.100    0.001     0.500   ];
%         % k
%         ig(2,:) = [1.000    0.0100    100.00  ];
%         % T2_B
%         ig(3,:) = [7.0e-6   1.0e-6    50.0e-6 ];
%         % T2_F
%         ig(4,:) = [100e-3    10e-3     200e-3 ];
%     end

    function ig = InitialGuessHuman()
        %         [ mean    min       max     ]
        % f
        ig(1,:) = [0.1414    0.001     0.500   ];
        % k
        ig(2,:) = [1.9828    0.0100    100.00  ];
        % T2_B
        ig(3,:) = [T2b   1.0e-6    50.0e-6 ];
        % T2_F
        ig(4,:) = [100e-3    10e-3     200e-3 ];
    end

    function ig = InitialGuessPhantom()
        %         [ mean    min        max     ]
        % f
        ig(1,:) = [0.00001  0.000001   0.5     ];
        % k
        ig(2,:) = [2.000    0.000100   100.00  ];
        % T2_B
        ig(3,:) = [T2b   1.0e-6    50.0e-6 ];
        % T2_F
        ig(4,:) = [50e-3    10e-3      4000e-3 ];
    end

% /END of MAIN
end
