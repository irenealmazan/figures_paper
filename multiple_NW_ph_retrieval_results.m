% this script shows the results

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


jitterlevel_1 = [100];
noiselevel_array_1 = [0];%[1 2 3];

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

flipflag_array = [1];

NW_phasecorrection_128angles;

for mm = 1:numel(noiselevel_array_1)
    
   
    
    noiselevel_str = num2str(noiselevel_array_1(mm));
    
    for jjj = 1:numel(jitterlevel_1)
        
       percent = jitterlevel_1(jjj);
       
       flipflag = flipflag_array(jjj);
       
        savefolder = [parentfolder 'allresults_blueshift/jitter_' num2str(jitterlevel_1(jjj)) '_noiselevel_' noiselevel_str];

        load([savefolder '/results.mat']);

        %load([parentfolder '/data_ERHIO_ini/struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(percent)]);

        
        %newobj.chi = struct_best_ERHIO.chi;
        
        struct_best_ERHIO.chi = newobj.chi;
        NW_figures_single_jitter;
        

    end
    
end