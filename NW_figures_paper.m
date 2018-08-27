% In this script the figures for the paper are generated:

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


%%%%% Figure 1 %%%%%%%%%%%%%%


fig_num = 1;
slice_array = [1:floor(round(numel(angles_list)/4)):numel(angles_list)];
[Psi_mod_matrix] = FiguresForPaper.display_slice_dp(NW,probe,delta_thscanvals,slice_array,ki_o,kf_o,d2_bragg,th,X,Y,Z,fig_num);

figure(2);
clf;

for jj = 1:numel(slice_array)
    %subplot(1,numel(slice_array),jj);
    figure;
    imagesc(Psi_mod_matrix(:,:,slice_array(jj)));
    axis image;
    colorbar;
    title(num2str(slice_array(jj)));
end

fig_num = 3;
DisplayFunctions.display_diff_geom(NW,ki,kf,qbragg,fig_num,X,Y,Z);


%% Fig 2: detrimental effect of non regular delta theta spacing in 3DFT for noiselevel 4:

jitter_summary = [0 5 10 20 40];
noiselevel_str = '0';

load('../results_files/Original_Sample.mat');
NW = img;

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';


% figure(90);
% clf;
% hold on;

jj = 1;
NW_struct(jj).rho_shift = NW;
support_abs = abs(NW_struct(jj).rho_shift);
support_final = (support_abs > 0.1*max(support_abs(:)));
NW_struct(jj).support_shift_fin = support_final;
NW_struct(jj).rho_nophase = NW_struct(jj).rho_shift.*conj(NW);
NW_struct(jj).percent = 0;
NW_struct(jj).dth_disp = zeros(size(NW,3),1);

% plot(NW_struct(jj).dth_disp,'r','LineWidth',3.0);
% legend_str{jj} = ['Original object'];


for jj = 2:numel(jitter_summary)+1
   load([parentfolder 'data_initial_alljitter/data_cnt_' noiselevel_str '_jitter_' num2str(jitter_summary(jj-1)) '.mat']);
   display(['percent' num2str(percent)])
   
%    for kk = 1:numel(FT_Proj_vol)
%       data_dp(:,:,kk) = FT_Proj_vol(kk).FT_Psij; 
%    end
    
   NW_struct(jj).percent = percent;
   NW_struct(jj).rho_shift = noise_NW*sqrt(mn/mncntrate);%fftshift(ifftn(fftshift(data_dp)));
   support_abs = abs(NW);
   support_final = (support_abs > 0.1*max(support_abs(:)));
   NW_struct(jj).support_shift_fin = support_final;
   NW_struct(jj).rho_nophase = NW_struct(jj).rho_shift.*conj(NW);
   NW_struct(jj).dth_disp = dth_disp;
   
   midpoint_1 = [round(size(noise_NW,1)/2)+1 round(size(noise_NW,2)/2)+1 round(size(noise_NW,3)/2)+1];

%    plot(NW_struct(jj).dth_disp,'-o');
%    legend_str{jj} = ['jitter = ' num2str(jitter_summary(jj-1))];
   
end
% figure(90);
% legend(legend_str);

fig_num = 40;
phase_color = [-2 2];
phase_color_2 = [-0.1 0.1];
intenscolor = [0 1.1];
dimension = 3;
%FiguresForPaper.figure5_bottompanel(NW_struct,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);

FiguresForPaper.figure2_noiseeffect(NW_struct,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);

figure(90);
clf;
hold on;
for jj = 1:numel(NW_struct)
    plot(NW_struct(jj).dth_disp,'-o');
   legend_str{jj} = ['jitter = ' num2str(NW_struct(jj).percent)];
   pause();
   
end
legend(legend_str);




%% Figure 1: low panel %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%% Figure 5: comparison of the quality of the reconstruction for
%%%%%%%% different jittering and different levels of noise

jitterlevel_summary = [0 5 10 20 40 500];
noiselevel_str = '0';

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/allresults_blueshift/';




for jj = 1:numel(jitterlevel_summary)
    
    savefolder = [parentfolder 'jitter_' num2str(jitterlevel_summary(jj)) '_noiselevel_' noiselevel_str ];
    load([savefolder '/results.mat']);

    [theta_iter] = DisplayResults.read_angles_iterations(data_exp,delta_thscanvals,delta_thscanvals);

    
    struct_err(jj).chi = [newobj.chi' errlist];
    %struct_err(kk).chid_direct = [err_ERHIO errlist_direct];
    struct_err(jj).rho = rho;
    struct_err(jj).support_iter = support_iter;
    struct_err(jj).rho_3DFT = rho_3DFT;
    struct_err(jj).support_new = support_new;
    struct_err(jj).mn = mn;
    struct_err(jj).mncntrate = mncntrate;
    struct_err(jj).theta_iter = theta_iter;
end


jj = numel(jitterlevel_summary) + 1;
struct_err(jj).chi = [zeros(numel([newobj.chi' errlist]))];
%struct_err(jj).chid_direct = [err_ERHIO errlist_direct];
struct_err(jj).rho = NW;
struct_err(jj).support_iter = abs(NW);
struct_err(jj).rho_3DFT = NW;
struct_err(jj).support_new = abs(NW);
struct_err(jj).mn = mn;
struct_err(jj).mncntrate = mncntrate;
struct_err(jj).theta_iter = [];

%mkdir results_sim_blueshift;

save([parentfolder 'struct_err_level' noiselevel_str '.mat'],'struct_err');

%%%%%%%%%%% error metric at different angular jittering percentages
jitterlevel_summary = [0 5 10 20 40 100 500];
noiselevel_str = '0';

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/allresults_blueshift/';

figure(200);
clf;
hold on;

for kk = 1:numel(struct_err)-1
   
    chi_final(kk) = struct_err(kk).chi(end)*sqrt(struct_err(kk).mn/struct_err(kk).mncntrate);
    %chi_direct_final(kk) = struct_err(kk).chid_direct(end);
    
    plot(log10(struct_err(kk).chi*sqrt(struct_err(kk).mn/struct_err(kk).mncntrate)),'LineWidth',3.0);
    legend_str{kk} = [num2str(jitterlevel_summary(kk)) ];
end
legend(legend_str);

figure(200);
savefig([parentfolder 'error_chi_alljiiter_' noiselevel_str '_scaled.fig']);


%%%%%%%%%%% error metric at different angular jittering percentages at the
%%%%%%%%%%% 2000 th iteration
jitterlevel_summary = [0 5 10 20 40 500];
noiselevel_str = '0';

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/allresults_blueshift/';



figure(2);
%subplot(121);
plot(jitterlevel_summary,log10(chi_final(1:end)),'-o');
%title('recip');
xlabel('% of angular jitter');ylabel('log(\epsilon)');
ax = gca;
set(ax,'FontSize',20);

figure(2);
savefig([parentfolder 'chi_vs_jitter_levelnoise_' noiselevel_str '_scaled.fig']);


%%%%%%%%%%% error metric at different angular jittering percentages at the
%%%%%%%%%%% 2000 th iteration with respect to the true object
jitterlevel_summary = [0 5 10 20 40];
noiselevel_str = '4';

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/allresults_blueshift/';


for jj = 1:numel(jitterlevel_summary)
    
    savefolder = [parentfolder 'jitter_' num2str(jitterlevel_summary(jj)) '_noiselevel_' noiselevel_str ];
    load([savefolder '/results.mat']);
    
    data_fin = data_exp;
    for kk = 1:numel(data_exp)
        data_fin(kk).I = data_exp(kk).simI/mn * mncntrate;
    end
    [errfinal(jj)] = DiffractionPatterns.calc_error_multiangle(probe, rho, data_fin,angles_list,ki_o,kf_o,X,Y,Z);
   
end

save([parentfolder 'errfinal_trueobject_' noiselevel_str],'errfinal');

figure(4);
%subplot(121);
plot(jitterlevel_summary,log10(errfinal(1:end)),'-o');
%title('recip');
xlabel('% of angular jitter');ylabel('log(\epsilon)');
ax = gca;
set(ax,'FontSize',20);

figure(4);
savefig([parentfolder 'chi_vs_jitter_levelnoise_' noiselevel_str '_trueobject.fig']);


%%%%%%%%%%% Object reconstructed at differnt jittering percentages and one
%%%%%%%%%%% noiselevel



% for noiselevel = 0;
phaseoffset = [1.49 1.54 1.69 1.545 1.4 0];%[1.55 1.58 1.58 1.545 1.01];
flipflag_list = [0 1 1 1 1 0];%[1 1 1 0 1];

for kk = 1:numel(struct_err)
    
    flipflag = flipflag_list(kk);
    
    if flipflag
        rho_plot = ifftn(conj(fftn(struct_err(kk).rho)));
        support_plot = abs(ifftn(conj(fftn(struct_err(kk).support_iter))));
    else
        rho_plot = struct_err(kk).rho;
        support_plot = struct_err(kk).support_iter;
    end
    
    
    rho_shift = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    
    support_shift = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
    midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];


    phase_rho_shift = angle(rho_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    
    
    
    struct_toplot(kk).rho_shift = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin;
    struct_toplot(kk).support_shift_fin = support_shift_fin;
    struct_toplot(kk).rho_nophase = rho_shift.*conj(NW)*exp(-1i*phase_rho_shift).*support_shift_fin*exp(-1i*phaseoffset(kk));
    
    
    
end

%DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(1).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);


fig_num = 30;
phase_color = [0 2.0];
phase_color_2 = [-0.1 0.1];
intenscolor = [0 0.027];
dimension = 3;
FiguresForPaper.figure5_bottompanel(struct_toplot,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);
% check phase offset

DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(1).rho_nophase,'','',[1 128],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(2).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(3).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(4).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);
DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(5).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);


figure(fig_num);
savefig([parentfolder 'allresults_blueshift/rho_alljitter_noiselevel' noiselevel_str '.fig']);

save([parentfolder 'allresults_blueshift/rhostructtoplot_level' noiselevel_str '.mat'],'struct_toplot','flipflag_list','phaseoffset');



%%%%%%%%%%%%%%%% Figure 2: effect of jittering in ER_HIO results


% for noiselevel = 1;
flipflag_list = [0 1 1 1 1 0];
%phaseoffset = [1.55 1.58 1.43 1.57 1.57];

% for noiselevel = 0;



phaseoffset_ERHIO =  [1.52 1.64 1.68 1.53 1.57 0];%[1.59 1.65 1.67 1.57 1.40];
%flipflag_list = [1 1 1 0 1];

for kk = 1:numel(struct_err)-1
    
    flipflag = flipflag_list(kk);
    
    if flipflag
        rho_plot = ifftn(conj(fftn(struct_err(kk).rho_3DFT)));
        support_plot = abs(ifftn(conj(fftn(struct_err(kk).support_new))));
    else
        rho_plot = struct_err(kk).rho_3DFT;
        support_plot = struct_err(kk).support_new;
    end
    
    
    rho_shift = DiffractionPatterns.shift_object(NW*sqrt(mncntrate/mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    
    support_shift = DiffractionPatterns.shift_object(abs(NW*sqrt(mncntrate/mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
   midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];

    
     phase_rho_shift = angle(rho_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
           
    
    struct_toplot(kk).rho_shift = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin;
    struct_toplot(kk).support_shift_fin = support_shift_fin;
    struct_toplot(kk).rho_nophase = rho_shift.*conj(NW)*exp(-1i*phase_rho_shift).*support_shift_fin*exp(-1i*phaseoffset_ERHIO(kk));
    
    %DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(kk).rho_nophase,'','',[40 90],[65 65],'23',32+kk);

    
end

%DisplayResults.compare_two_objects(NW*sqrt(mncntrate/mn).*conj(NW),struct_toplot(5).rho_nophase,'','',[40 90],[midpoint_1(2) midpoint_1(3)],'23',42);


fig_num = 40;
phase_color = [0 2.0];
phase_color_2 = [-0.1 0.1];
intenscolor = [0 0.027];
dimension = 3;
FiguresForPaper.figure5_bottompanel(struct_toplot,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);

figure(40);
savefig([parentfolder 'allresults_blueshift/rhoERHIO_alljitter_noiselevel' noiselevel_str '.fig']);



%%%%%% Compare the angle retrieval for different degrees of noise levels:


jitterlevel_single = [10];
noiselevel_summary = [3 1 4 0];
legend_summary = {'noise, cnt\_rate = 1e2','noise, cnt\_rate = 1e3','noise, cnt\_rate = 1e4','no noise'};

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';


for jj = 1:numel(noiselevel_summary)
    
    noiselevel_str = num2str(noiselevel_summary(jj));
    
    savefolder = [parentfolder 'allresults_blueshift/jitter_' num2str(jitterlevel_single) '_noiselevel_' noiselevel_str ];
    load([savefolder '/results.mat']);

    [theta_iter] = DisplayResults.read_angles_iterations(data_exp,delta_thscanvals,delta_thscanvals);
    
    
    
    struct_err(jj).chi = [newobj.chi' errlist];
    %struct_err(kk).chid_direct = [err_ERHIO errlist_direct];
    struct_err(jj).rho = rho;
    struct_err(jj).support_iter = support_iter;
    struct_err(jj).rho_3DFT = rho_3DFT;
    struct_err(jj).support_new = support_new;
    struct_err(jj).mn = mn;
    struct_err(jj).mncntrate = mncntrate;
    struct_err(jj).theta_iter = theta_iter;
    struct_err(jj).NW = NW;
    
    theta_noise_struct(jj).theta_iter = theta_iter;
    theta_noise_struct(jj).legendstr = legend_summary{jj};
    theta_noise_struct(jj).dth_disp = dth_disp;
    
    DisplayResults.display_all_angles_oneiterations_errorrel(theta_noise_struct(jj).theta_iter,data_exp,dth_disp,[1 cnt_ntheta],'absolute',1000+jj);

    
end


DisplayResults.display_all_angles_oneiterations_errorrel_levelnoise(theta_noise_struct,data_exp,[2000 2000 2000 2000],'absolute',1027);

figure(1027);
savefig([parentfolder 'allresults_blueshift/error_correction_alllevels_jitter_' num2str(jitterlevel_single) '.fig' ]);

%%%%% Compare the image quality for the same level of jittering and
%%%%% different leves of noise

%jitterlevel_single = [40];
%noiselevel_summary = [3 1 4 0];

%flipflag_list = [0 1 1 1]; % 0 percent
%flipflag_list = [1 0 0 1]; % 5 percent
flipflag_list = [1 0 0 0]; % 10 percent
%flipflag_list = [1 0 1 0]; % 20 percent
%flipflag_list = [1 0 0 1]; % 40 percet

phaseoffset = [1.64 1.64 1.53 1.54];

for jj = 1:numel(noiselevel_summary)
    
   

    
     flipflag = flipflag_list(jj);
    
    if flipflag
        rho_plot = ifftn(conj(fftn(struct_err(jj).rho)));
        support_plot = abs(ifftn(conj(fftn(struct_err(jj).support_iter))));
    else
        rho_plot = struct_err(jj).rho;
        support_plot = struct_err(jj).support_iter;
    end
    
    
    rho_shift = DiffractionPatterns.shift_object(NW*sqrt(struct_err(jj).mncntrate/struct_err(jj).mn),rho_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    
    
    support_shift = DiffractionPatterns.shift_object(abs(NW*sqrt(struct_err(jj).mncntrate/struct_err(jj).mn)),support_plot,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
    support_shift_abs = abs(support_shift);
    support_shift_fin = (support_shift_abs>0.1*max(support_shift_abs(:)));
    
    midpoint_1 = [round(size(rho_shift,1)/2)+1 round(size(rho_shift,2)/2)+1 round(size(rho_shift,3)/2)+1];


    phase_rho_shift = angle(rho_shift(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    phase_NW = angle(NW(midpoint_1(1),midpoint_1(2),midpoint_1(3)));
    
    
    
    struct_toplot(jj).rho_shift = rho_shift*exp(-1i*phase_rho_shift).*support_shift_fin*sqrt(struct_err(jj).mn/struct_err(jj).mncntrate);
    struct_toplot(jj).support_shift_fin = support_shift_fin;
    struct_toplot(jj).rho_nophase = rho_shift.*conj(NW)*exp(-1i*phase_rho_shift).*support_shift_fin*exp(-1i*phaseoffset(jj));
    struct_toplot(jj).NW = NW;
end

fig_num = 60;
phase_color = [0 2.0];
phase_color_2 = [-0.1 0.1];
intenscolor = [0 1];
dimension = 3;
FiguresForPaper.object_samejitter_differentnoise(struct_toplot,intenscolor,phase_color,phase_color_2,[40 90 40 90],[midpoint_1(dimension)],num2str(dimension),fig_num);

fig_num = 70;
FiguresForPaper.object_samejitter_differentnoise_support(struct_toplot,X,Y,Z,fig_num);

fig_num = 80;
FiguresForPaper.object_samejitter_differentnoise_realobject(struct_toplot,X,Y,Z,fig_num);

figure(60);
savefig([parentfolder 'allresults_blueshift/object_support_alllevels_jitter_' num2str(jitterlevel_single) '.fig' ]);


figure(70);
savefig([parentfolder 'allresults_blueshift/object_support_3Dalllevels_jitter_' num2str(jitterlevel_single) '.fig' ]);

figure(80);
savefig([parentfolder 'allresults_blueshift/object_trueobject_3Dalllevels_jitter_' num2str(jitterlevel_single) '.fig' ]);


%%%%%% Plot all the jitter for all the noise levels:


jitterlevel_summary = [0 5 10 20 40];
noiselevel_summary = [3 1 4 0];
legend_summary = {'noise, cnt\_rate = 1e2','noise, cnt\_rate = 1e3','noise, cnt\_rate = 1e4','no noise'};

parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';


figure(200);
clf;
hold on;
for jj = 1:numel(noiselevel_summary)
    load([parentfolder 'allresults_blueshift/struct_err_level' num2str(noiselevel_summary(jj)) '.mat']);
    
    for kk = 1:numel(struct_err)-1
        
        struct_noise(jj).chi_final(kk) = struct_err(kk).chi(end)*sqrt(struct_err(kk).mn/struct_err(kk).mncntrate);
        %chi_direct_final(kk) = struct_err(kk).chid_direct(end);
        
       % plot(log10(struct_err(kk).chi*sqrt(struct_err(kk).mn/struct_err(kk).mncntrate)),'LineWidth',3.0);
    end
    
    figure(200);
    plot(jitterlevel_summary,log10( struct_noise(jj).chi_final(1:end)),'-o');
   %legend_str{jj} = [legend_summary(jj) ];

end
legend(legend_summary);
xlabel('angular jitter %');
ylabel('log(\epsilon)');
set(gca,'FontSize',20);

figure(1);
savefig([parentfolder 'allresults_blueshift/chi_vs_jitter_allnoise_scaled.fig' ]);

