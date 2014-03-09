% qMT_MODEL_FIT_FINAL_GETMODE Get qMT mode
% mode = qMT_model_fit_final_getmode(fv) simply returns the mode of the
% qMT_model_fit_final used, based on the 4th dimension of fv.
%
% If fv has 5 parameters, mode is 0
% If fv has 3 parameters, mode is 1
% If fv has 6 parameters, mode is 2
% If fv has 4 parameters, mode is 3

function mode = qMT_model_fit_final_getmode(fv)

switch size(fv,4)
    case 5
        mode = 0;
    case 3
        mode = 1;
    case 6
        mode = 2;
    case 4
        mode = 3;
end

end