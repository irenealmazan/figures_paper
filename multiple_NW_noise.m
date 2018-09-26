% this script shows the results

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions/'));

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

addpath(genpath(parentfolder));

jitterlevel_3 = [0 0 10 50 100 500];
mncrate_array_2 = [1e3 1e3 1e2 1e4];
noiseflag_array_2 = [0 1 1 1 1];
noiselevel_array_2 = [0 1 3 4];
%%%% no noise


NW_flags;

flagOnlyPlot = 0;

for mmm = 1:numel(mncrate_array_2)
    
  
    
    mncntrate = mncrate_array_2(mmm);
    
    noiselevel_str = num2str(noiselevel_array_2(mmm));
    
    load('../results_files/Original_Sample.mat');
    NW = img;
    kkkk = 1;
    
    midpoint_1 = [round(size(NW,1)/2)+1 round(size(NW,2)/2)+1 round(size(NW,3)/2)+1];
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    
    
    struct_noise(mmm).struct_toplot(kkkk).rho_shift = NW*exp(-1i*phase_NW);
    struct_noise(mmm).struct_toplot(kkkk).support_shift_fin = abs(NW);
    struct_noise(mmm).struct_toplot(kkkk).rho_nophase = NW.*conj(NW);
    struct_noise(mmm).struct_toplot(kkkk).mn = 1e10;
    struct_noise(mmm).struct_toplot(kkkk).mncrate = 1e10;
    struct_noise(mmm).struct_toplot(kkkk).middpind = 65;
    
    
    for kkkk = 2:numel(jitterlevel_3)
        
        percent = jitterlevel_3(kkkk);
        
        load([parentfolder '/data_initial_alljitter/data_cnt_0_jitter_' num2str(percent) '.mat']);
        display(num2str(kkkk))
        display(num2str(jitterlevel_3))
        noiseflag = noiseflag_array_2(mmm);
        
        NW_add_dp_noise;
        
        midpoint_1 = [round(size(NW,1)/2)+1 round(size(NW,2)/2)+1 round(size(NW,3)/2)+1];
        phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
        
        phaseoffset_rho = angle(noise_NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
        
        struct_noise(mmm).struct_toplot(kkkk).rho_shift = noise_NW*exp(-1i*phaseoffset_rho);
        struct_noise(mmm).struct_toplot(kkkk).support_shift_fin = abs(NW);
        struct_noise(mmm).struct_toplot(kkkk).rho_nophase = noise_NW.*conj(NW).*struct_noise(mmm).struct_toplot(kkkk).support_shift_fin;
        struct_noise(mmm).struct_toplot(kkkk).mn = mn;
        struct_noise(mmm).struct_toplot(kkkk).mncrate = mncntrate;
        struct_noise(mmm).struct_toplot(kkkk).middpind = middpind;
        struct_noise(mmm).struct_toplot(kkkk).snr = snr;
        
    end
    
end


for mmm = 1:numel(mncrate_array_2)
    
    
    fig_num = 300+mmm;
    phase_color = [0 2];
    phase_color_2 = [-0.02 0.02];
    intenscolor = [0 1.2];
    dimension = 3;
    FiguresForPaper.figure5_bottompanel(struct_noise(mmm).struct_toplot,intenscolor,phase_color,phase_color_2,[45 85 45 85],[midpoint_1(dimension)],num2str(dimension),fig_num);

end

for mmm = 1:numel(mncrate_array_2)
    noiselevel_str = num2str(noiselevel_array_2(mmm));
    
    fig_num = 300+mmm;
    %figure(30);
    savefig([parentfolder 'allresults_blueshift/alljitter_noiselevel' noiselevel_str '.fig']);
    
end
