% In this script to create the  figure 3 for the paper:

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

jitterlevel_1 = 10;
noiselevel_str = '0';

savefolder = [parentfolder 'allresults_blueshift/jitter_' num2str(jitterlevel_1) '_noiselevel_' noiselevel_str];

load([savefolder '/results.mat']);

%%%% Check for the flipflag:

DisplayResults.show_rho_theta_update(6,errlist,rho.*support_iter,midsl,angles_list,delta_thscanvals'+dth_disp,norm_grad_rho(1:nrho),beta_rho(1:nrho),norm_grad_theta(1:cnt_ntheta-1),beta_theta(1:cnt_ntheta-1),'theta');


flipflag = 0;

%%% Panel a)



if flipflag 
    rho_ini_plot = ifftn(conj(fftn(rho_ini)));
    rho_plot = ifftn(conj(fftn(rho)));
    support_plot = abs(ifftn(conj(fftn(support_iter))));
else
    rho_ini_plot = rho_ini;
    rho_plot = rho;
    support_plot = support_iter;
end


[rho_shift,shift_direct_space] = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
[rho_ini_shift,shift_direct_space_ini] = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_ini_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);


[support_shift,shift_support_vect] = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z); 
support_shift_abs = abs(support_shift);
support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));

midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];


phase_rho_shift = angle(rho_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
phase_rho_ini = angle(rho_ini_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));

phaseoffset_rho = 1.54;
phaseoffset_rho_ERHIO = 1.54;

% structure:
NW_struct.NW1 = NW*exp(-1i*phase_NW);
NW_struct.NW1_nophase = NW.*conj(NW);
NW_struct.NW2 = rho_ini_shift.*exp(-1i*phase_rho_ini).*support_shift_fin/sqrt(mncntrate/mn);
NW_struct.NW2_nophase = NW_struct.NW2.*conj(NW)*exp(-1i*phase_rho_ini).*exp(-1i*phaseoffset_rho_ERHIO);
NW_struct.NW3 = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin/sqrt(mncntrate/mn);
NW_struct.NW3_nophase = NW_struct.NW3.*conj(NW).*exp(-1i*phaseoffset_rho);


% test phase ofset:


intenscolor = [0 0.03];
phasecolor = [-2 2];

DisplayResults.compare_two_objects(NW_struct.NW1_nophase,NW_struct.NW2_nophase,'','',intenscolor,phasecolor,[1 128],[midpoint_1(2) midpoint_1(3)],'23',2);
DisplayResults.compare_two_objects(NW_struct.NW1_nophase,NW_struct.NW3_nophase,'object','retrieved object',intenscolor,phasecolor,[1 128],[midpoint_1(2) midpoint_1(3)],'23',3);





% figures:

intenscolor = [0 1.1];
phasecolor = [0 2];
phasecolor2 = [-0.1 0.1];
dimension = 3;
FiguresForPaper.figure4_panel_a(NW_struct,intenscolor,phasecolor,phasecolor2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),4);

savefig([parentfolder 'pictures_paper/fig4_phaseobject_v2.fig']);

%%% Panel b)

[theta_iter] = DisplayResults.read_angles_iterations(data_exp,delta_thscanvals,delta_thscanvals);

DisplayResults.display_all_angles_oneiterations_errorrel(theta_iter,data_exp,dth_disp,[1 cnt_ntheta],'absolute',1025);
