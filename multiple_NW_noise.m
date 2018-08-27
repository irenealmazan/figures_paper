% this script shows the results

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions/'));

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

addpath(genpath(parentfolder));

jitterlevel_2 = [0 0];
mncrate_array = [1e3 1e2 1e4];
noiseflag_array = [1 1 1];
noiselevel_array_2 = [1 3 4];
%%%% no noise


NW_flags;

for mmm = 1:numel(mncrate_array)
    
    noiseflag = noiseflag_array(mmm);
    
    mncntrate = mncrate_array(mmm);
    
    noiselevel_str = num2str(noiselevel_array_2(mmm));
    
    load('../results_files/Original_Sample.mat');
    NW = img;
    jjjj = 1;
    
    midpoint_1 = [round(size(NW,1)/2)+1 round(size(NW,2)/2)+1 round(size(NW,3)/2)+1];
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    
    
    struct_noise(mmm).struct_toplot(jjjj).rho_shift = NW*exp(-1i*phase_NW);
    struct_noise(mmm).struct_toplot(jjjj).support_shift_fin = abs(NW);
    struct_noise(mmm).struct_toplot(jjjj).rho_nophase = NW.*conj(NW);
    struct_noise(mmm).struct_toplot(jjjj).mn = 1e10;
    struct_noise(mmm).struct_toplot(jjjj).mncrate = 1e10;
    struct_noise(mmm).struct_toplot(jjjj).middpind = 65;
    
    
    for jjjj = 2:numel(jitterlevel_2)
        
        percent = jitterlevel_2(jjjj);
        
        load([parentfolder '/data_initial_alljitter/data_cnt_' noiselevel_str '_jitter_' num2str(percent) '.mat']);
        NW_add_dp_noise;
        
        midpoint_1 = [round(size(NW,1)/2)+1 round(size(NW,2)/2)+1 round(size(NW,3)/2)+1];
        phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
        
        phaseoffset_rho = angle(noise_NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
        
        struct_noise(mmm).struct_toplot(jjjj).rho_shift = noise_NW*exp(-1i*phaseoffset_rho);
        struct_noise(mmm).struct_toplot(jjjj).support_shift_fin = abs(NW);
        struct_noise(mmm).struct_toplot(jjjj).rho_nophase = noise_NW.*conj(NW).*struct_noise(mmm).struct_toplot(jjjj).support_shift_fin;
        struct_noise(mmm).struct_toplot(jjjj).mn = mn;
        struct_noise(mmm).struct_toplot(jjjj).mncrate = mncntrate;
        struct_noise(mmm).struct_toplot(jjjj).middpind = middpind;
        struct_noise(mmm).struct_toplot(jjjj).snr = snr;
        
        figure_num = jjjj*10;
        
        intesncolor = [0 0.15];
        phasecolor = [0 2];
        
        dimension = '3';
        %DisplayResults.compare_two_objects(NW*exp(-1i*phase_NW)/sqrt(mn/mncntrate),noise_NW*exp(-1i*phaseoffset_rho),'Original object',['Compatible object at ' num2str(percent) ' % jitter no noise'],intesncolor,phasecolor,[50 80 50 80],64,dimension,figure_num);
        DisplayResults.compare_two_objects(NW*exp(-1i*phase_NW),noise_NW*exp(-1i*phaseoffset_rho),'Original object',['Compatible object at ' num2str(percent) ' % jitter no noise'],intesncolor,phasecolor,[50 80 50 80],64,dimension,figure_num);
        
        dimension = '13';
        %DisplayResults.compare_two_objects(NW*exp(-1i*phase_NW)/sqrt(mn/mncntrate),noise_NW*exp(-1i*phaseoffset_rho),'Original object',['Compatible object at ' num2str(percent) ' % jitter no noise'],intesncolor,phasecolor,[50 80],[65 65],dimension,figure_num+1);
        DisplayResults.compare_two_objects(NW*exp(-1i*phase_NW),noise_NW*exp(-1i*phaseoffset_rho),'Original object',['Compatible object at ' num2str(percent) ' % jitter no noise'],intesncolor,phasecolor,[50 80],[65 65],dimension,figure_num+1);
        
        DisplayResults.compare_two_objects(NW.*conj(NW),noise_NW.*conj(NW),'Original object',['Compatible object at ' num2str(percent) ' % jitter no noise'],intesncolor,phasecolor,[50 80],[65 65],dimension,figure_num+2);
        
        DisplayResults.show_rocking_curve(delta_thscanvals+dth_disp',rock_curve_noise,'new',figure_num+3,'r','true using 2D FT');
        
        
    end
    fig_num = 30+mmm;
    phase_color = [0 2];
    phase_color_2 = [-0.2 0.2];
    intenscolor = [0 1.2];
    dimension = 3;
    FiguresForPaper.figure5_bottompanel(struct_noise(mmm).struct_toplot,intenscolor,phase_color,phase_color_2,[45 85 45 85],[midpoint_1(dimension)],num2str(dimension),fig_num);

end





%figure(30);
%savefig([parentfolder 'allresults_blueshift/rho_extremejitter_noiselevel' noiselevel_str '.fig']);


