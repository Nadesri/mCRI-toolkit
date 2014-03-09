% qMT_model_fit_final_par Parallel Processing qMT Wrapper Function
% 
% [fv rnrm] = qMT_model_fit_final_par(vargin) accepts all input arguments
% of qMT_model_fit_final normally, and produces the same ouputs fv and
% rnrm. qMT_model_fit_final_par is simply a wrapper function that divides
% data_MT, B0, B1, B1_MT, and data_SPGR(or R1/PD) into the number of
% available labs for a parfor function.
%
% As a fair warning, qMT_model_fit_final_par does NOT perform parity checks
% (correct sizes/dimensions of inputs). Therefore, it is suggested that 
% qMT_model_fit_final be used first to confirm parity, and then update the
% function call to qMT_model_fit_final_par.
%
% For more help, see qMT_model_fit_final. qMT_model_fit_final was written
% and developed by Pouria Mossahebi.
%
% Nade Sritanyaratana
% University of Wisconsin-Madison
% October 09, 2013
% v1.0

function [fv rnrm] = qMT_model_fit_final_par( data_MT, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0, B1, B1_MT,  mode, in15, in16, TR_SPGR, t_m_SPGR)

    labsize = matlabpool('size');
    if ~labsize
        matlabpool;
    else
        disp(['Parallel functionality using ' int2str(labsize) ' labs.']);
    end

    sectionsize = round(size(data_MT,1)/labsize);
    
    % Initialize setions
    data_MTp = cell([labsize 1]);
    B0p      = cell([labsize 1]);
    B1p      = cell([labsize 1]);
    B1_MTp   = cell([labsize 1]);
    in15p    = cell([labsize 1]);
    in16p    = cell([labsize 1]);
    % Separate data into sections
    for i=1:labsize
        if i~=labsize
            parindices = ((i-1)*sectionsize+1):(i*sectionsize);
        else
            parindices = ((i-1)*sectionsize+1):size(data_MT,1);
        end
        data_MTp{i} = data_MT(parindices,:);
        B0p{i}    = B0(parindices);
        B1p{i}    = B1(parindices);
        B1_MTp{i} = B1_MT(parindices);
        if (mode==0)||(mode==2)
            in15p{i}  = in15(parindices,:); %data_SPGR
            in16p{i}  = in16; %alpha_SPGR
        else
            in15p{i}  = in15(parindices); %PD
            in16p{i}  = in16(parindices); %R1
        end
    end
    
    parfor i=1:labsize
        [fvp{i} rnrmp{i}] = qMT_model_fit_final(data_MTp{i}, MT_Off_Freq, alpha_MT, alpha_0, t_r, t_m, t_s, pulse_type, line_shape, table_update, B0p{i}, B1p{i}, B1_MTp{i},  mode, in15p{i}, in16p{i}, TR_SPGR, t_m_SPGR);
    end
    
    % recombine fv and rnrm from fvp and rnrmp
    nparam = size(fvp{1},2);
    fv   = zeros([size(data_MT,1) nparam]);
    rnrm = zeros([size(data_MT,1) 1]);
    
    for i=1:labsize
        
        if i~=labsize
            parindices = ((i-1)*sectionsize+1):(i*sectionsize);
        else
            parindices = ((i-1)*sectionsize+1):size(data_MT,1);
        end
        
        for i_param = 1:nparam
            fv(parindices,i_param) = fvp{i}(:,i_param);
        end
        
        rnrm(parindices) = rnrmp{i};
    end
    
end
