addpath(genpath('./m_scripts/'));
addpath(genpath('./calc_functions'));


jitterlevel_1 = [5];%[0 5 10 20 40];
mncrate_array_1 = [1e3];%[2e2 2e2 1e3];
noiseflag_array_1 = [1];%[0 1 1];
noiselevel_array_1 = [1];%[1 2 3]
%flipflag_array = [0 0 1 1 1];
flipflag_array = [1];
%%%% no noise

folderpath =  '../Theta_annealing_blueshift_multiplesim_paperFig3/data_ERHIO_ini/';

counter = 1;

for mm = 1:numel(mncrate_array_1)
    
    noiselevel_str = num2str(noiselevel_array_1(mm));
    mncntrate = mncrate_array_1(mm);
    addNWstrain = noiseflag_array_1(mm);
    
    for jjj = 1:numel(jitterlevel_1)
        
        percent = jitterlevel_1(jjj);
        
        plotResults = 0;
        NW_make_mesgrhid_and_kvectors;
        
        
        
        load([folderpath 'struct_ERHIO_ini' noiselevel_str '_jitter_' num2str(percent)]);
        
        figure(mm);
        clf;
        hold on;
        for kk = 1:numel(struct_ER_HIO)
            figure(mm);
            plot(log10(struct_ER_HIO(kk).chi),'LineWidth',3.0);
            legend_str{kk} = ['noiselevel = ' noiselevel_str ' jitter = ' num2str(percent) 'loop kk = ' num2str(kk)];
            
            
            flipflag = flipflag_array(jjj);
            
            if flipflag
                rho_3DFT_toplot= (ifftn(conj(struct_ER_HIO(kk).dp)));
            else
                rho_3DFT_toplot= (ifftn(struct_ER_HIO(kk).dp));
            end
            
            rho_2DFT = DiffractionPatterns.From3DFT_to_2DFT(rho_3DFT_toplot,delta_thscanvals,probe,ki_o,kf_o,X,Y,Z);
            [rho_2DFT_shift,rho_2DFT_shift_vector] = DiffractionPatterns.shift_object(NW,rho_2DFT,delta_thscanvals,ki_o,kf_o,kf_o-ki_o,d2_bragg,X,Y,Z);
            
            
            midpoint = [round(size(rho_2DFT_shift,1)/2)+1 round(size(rho_2DFT_shift,2)/2)+1 round(size(rho_2DFT_shift,3)/2)+1];
            
            dimension = 3;
            DisplayResults.compare_two_objects(NW,rho_2DFT_shift,'','',[40 90 40 90],[midpoint(dimension)],num2str(dimension),42+counter);
            %DisplayResults.compare_two_objects(NW,rho_ini,'','',[40 90 40 90],[65],'3',100);

            
            counter = counter + 1;
        end
    end
    
    figure(mm);
    legend(legend_str);
    
end