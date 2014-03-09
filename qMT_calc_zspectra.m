% qMT_CALC_ZSPECTRA Calculate Theoretical Z-Spectra
% S1 = qMT_calc_zspectra(mtoffs,mtfs,initvars) calculates the theoretical
% qMT curves for a given set of qMT properties and pulse sequence
% parameters. qMT_calc_zspectra assumes the mCRI model proposed by
% Mossahebi P and Samsonov A in MRM 2013. 

function Z = qMT_calc_zspectra(mtoffs, mtfs, initvars)
% Inputs: mtfs, mtoff, initvars

% Initialize SuperLorentian Table
superLorentian_fastPM([],[],1,0,initvars.t_m);

PD=initvars.PD;     R1=initvars.R1;   R1_B = initvars.R1_B;
f=initvars.f;       k=initvars.k;               T2b=initvars.T2b;
T2f=initvars.T2f;   t_m=initvars.t_m;           t_s=initvars.t_s;
t_r=initvars.t_r;	alphaMTR=initvars.alphaMTR;

if isfield(initvars,'mode')
    mode = initvars.mode;
else
    mode = 0;
end

% alpha_SPGR = initvars.alpha_SPGR;
% TR_SPGR = initvars.TR_SPGR;

w_1rms = getB1eff_vasily(mtfs, t_m, 'fermi');
W = pi * (w_1rms).^2;


%     g_B = find_g_B_NS(mtoffs,T2b);
g_B = superLorentian_fastPM(mtoffs,T2b);
W_B = W .*g_B;

MToff = qMT_Signal_final(f, k, alphaMTR/180*pi, 0, 0, t_r+t_m+t_s, PD, R1, 0,0);

% -- Signal
for kk=1:numel(mtfs)
    if (mode==0)||(mode==1)
        S1(kk) = qMT_Signal_final(f, k, alphaMTR/180*pi, t_m, t_s, t_r, PD, R1, W_B(kk),0);
    elseif (mode==2)||(mode==3)
        W_F = calc_W_F(w_1rms,mtoffs,f,k,R1,R1_B,T2f);
        S1(kk) = qMT_Signal_final(f, k, alphaMTR/180*pi, t_m, t_s, t_r, PD, R1, W_B(kk),W_F(kk));
    else
        error('mode not specified. Attempt to default to mode=0 failed. Check code.')
    end
end

Z = S1/MToff;

% Close up SuperLorentian Table
superLorentian_fastPM([],[],-1);

end