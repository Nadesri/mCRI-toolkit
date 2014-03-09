disp('Initializing: Recon')
Recon;
disp('Recon: done!');

disp('Initializing: Resizing')
qMT_imresize;
disp('Resizing: Done!')

disp('Initializing: Coregistration')
CoregAnalyze;
disp('Coregistration: Done!')

disp('Initializing: ProcessingAS')
ProcessingGlobAS;
disp('ProcessingAS: Done!')

disp('Initlializing: qMT Analysis')
qMT_analyze
disp('qMT Analysis: Done!')