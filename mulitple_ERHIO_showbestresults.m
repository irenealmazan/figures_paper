addpath(genpath('../m_scripts/'));
addpath(genpath('../calc_functions'));

parentfolderstr = '../Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

addpath(genpath(parentfolderstr));

folderstr = [parentfolderstr '/data_ERHIO_ini/'];



jitterlevel_1 = [100];
mncrate_array_1 = [1e4];%[2e2 2e2 1e3];
noiseflag_array_1 = [1];%[0 1 1];
noiselevel_array_1 = [4];%[1 3 4];%[1 2 3]
%%%% no noise
flipflag_array = [0];

counter = 1;

for mm = 1:numel(noiselevel_array_1)
    
    noiselevel_str = num2str(noiselevel_array_1(mm));
    mncntrate = mncrate_array_1(mm);
    addNWstrain = 1;
    figure(100);
    clf;
    hold on;
    for jjj = 1:numel(jitterlevel_1)
        
        percent = jitterlevel_1(jjj);
        
        plotResults = 0;
        NW_make_mesgrhid_and_kvectors;
        
        
        
        load([folderstr 'struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(percent)]);
        
        %figure(mm);
        figure(100);
        plot(log10(struct_best_ERHIO.chi),'LineWidth',3.0);
        legend_str{jjj} = ['noiselevel = ' noiselevel_str ' jitter = ' num2str(percent)];
        
        flipflag = flipflag_array(jjj);
        
        if flipflag
            rho_3DFT_toplot= (ifftn(conj(struct_best_ERHIO.dp)));
        else
            rho_3DFT_toplot= (ifftn(struct_best_ERHIO.dp));
        end
        
        rho_2DFT = DiffractionPatterns.From3DFT_to_2DFT(rho_3DFT_toplot,delta_thscanvals,probe,ki_o,kf_o,X,Y,Z);
        [rho_2DFT_shift,rho_2DFT_shift_vector] = DiffractionPatterns.shift_object(NW,rho_2DFT,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
        
        
        midpoint = [round(size(rho_2DFT_shift,1)/2)+1 round(size(rho_2DFT_shift,2)/2)+1 round(size(rho_2DFT_shift,3)/2)+1];

        dimension = 3;
        intenscolor = [0 0.1];
        phasecolor = [-2 2];
        DisplayResults.compare_two_objects(NW/sqrt(mncntrate/mm),rho_2DFT_shift,'','',intenscolor,phasecolor,[40 90 40 90],[midpoint(dimension)],num2str(dimension),72+counter);
        DisplayResults.compare_two_objects(NW/sqrt(mncntrate/mm),rho_2DFT_shift,'','',intenscolor,phasecolor,[40 90],[midpoint(dimension) midpoint(dimension)],'13',82+counter);
        %DisplayResults.compare_two_objects(ifftn(struct_best_ERHIO.dp),NW,'','',[1 128 1 128],[35],'3',52+counter);

        counter = counter + 1;
    end
    
    figure(mm);
    legend(legend_str);
    
end