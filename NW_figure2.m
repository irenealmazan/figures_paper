% In this script to create the  figure 2 for the paper:

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

jitterlevel_1 = 10;
noiselevelvector = [0 1 3 4];



for kkkk = 1:numel(noiselevelvector)
    
    noiselevel_str = num2str(noiselevelvector(kkkk));
    savefolder = [parentfolder 'allresults_blueshift/jitter_' num2str(jitterlevel_1) '_noiselevel_' noiselevel_str];
    
    load([savefolder '/results.mat']);
    
    if flipflag
        rho_ini_plot = ifftn(conj(fftn(rho_ini)));
        support_plot = abs(ifftn(conj(fftn(support_iter))));
    else
        rho_ini_plot = rho_ini;
        support_plot = support_iter;
    end
    
    [rho_ini_shift,shift_direct_space_ini] = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_ini_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    [support_shift,shift_support_vect] = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
    midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];
    
    
    phase_rho_ini = angle(rho_ini_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    
    phaseoffset_rho_ERHIO = 1.54;
    
    NW_ERHIO(kkkk).NW1 = NW;
    NW_ERHIO(kkkk).NW1_nophase = NW.*conj(NW);
    NW_ERHIO(kkkk).NW2 = rho_ini_shift.*exp(-1i*phase_rho_ini).*support_shift_fin/sqrt(mncntrate/mn);
    NW_ERHIO(kkkk).NW2_nophase = NW_ERHIO(kkkk).NW2.*conj(NW).*exp(-1i*phaseoffset_rho_ERHIO);
end


