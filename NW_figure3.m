% In this script to create the  figure 3 for the paper:

addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));


parentfolder = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

jitterlevel_1 = 10;
noiselevel_str = '1';

savefolder = [parentfolder 'allresults_blueshift/jitter_' num2str(jitterlevel_1) '_noiselevel_' noiselevel_str ];

load([savefolder '/results.mat']);





h2 = figure(7);
clf;
plot(log10([struct_best_ERHIO.chi' errlist]),'LineWidth',3.0);
hold on;
plot(log10(struct_best_ERHIO.chi),'LineWidth',3.0);
%title('error metric recip. space');
xlabel('Iterations');
ylabel('log(\epsilon)');

hold on;plot([326:426],log10(struct_best_ERHIO.chi(326:426)'),'LineWidth',3.0,'Color','k')
hold on;plot([610:710],log10(struct_best_ERHIO.chi(610:710)'),'LineWidth',3.0,'Color','k')
hold on;plot([882:932],log10(struct_best_ERHIO.chi(882:932)'),'LineWidth',3.0,'Color','k')
hold on;plot([331:431],log10(struct_best_ERHIO.chi(331:431)),'LineWidth',3.0,'Color','k')
hold on;plot([611:711],log10(struct_best_ERHIO.chi(611:711)),'LineWidth',3.0,'Color','k')
hold on;plot([891:941],log10(struct_best_ERHIO.chi(891:941)'),'LineWidth',3.0,'Color','k')


ax = gca;
set(ax,'FontSize',20);