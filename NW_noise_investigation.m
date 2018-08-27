

jitterval = [0 5 10 20 40];
noiselevelarray_4 = [0 1 4];
cntrate_array_4 = [1e3 1e3 1e4];

folder = '/Users/ialmazn/Box Sync/Nanowire_ptychography/NSLS II/NSLS II March 2017/Irene_Analysis/Theta_annealing_blueshift_multiplesim_paperFig3_128angles/';

for mmm = 1:numel(noiselevelarray_4)
    
    load([folder 'allresults_blueshift/jitter_' num2str(jitterval()) '_noiselevel_' num2str(noiselevelarray_4(mmm)) '/results.mat']);
    
    [Psi_mod_noiselevel(mmm).Psi_mod,rock_data_exp_noiselevel(mmm).rc] = DiffractionPatterns.calc_rock_curve_2DFT(rho,probe,angles_list,ki_o,kf_o,kf_o-ki_o,X,Y,Z);
     
    for jj = 1:numel(data_exp)        
        err_per_slice(jj) = sum((sqrt(Psi_mod_noiselevel(mmm).Psi_mod(jj).I(:))-sqrt(data_exp(jj).I(:))).^2./numel(data_exp(jj).I));             
    end
        
    errperslice_noiselevel(mmm).err_per_slice = err_per_slice;
    
    error_noiselevel(mmm).errlist = [newobj.chi' errlist];
    
    mn_noiselevel(mmm).mn = mn;
end

 figure(1);clf;hold on;
for mmm = 1:numel(noiselevelarray_4)
    plot(delta_thscanvals,errperslice_noiselevel(mmm).err_per_slice);
    legend_str{mmm} = ['noiselevel = ' num2str(noiselevelarray_4(mmm)) 'count rate = ' num2str(cntrate_array_4(mmm))]; 
end

legend(legend_str);


mmm = 3;

figure(400);clf
    for jj = numel(data_exp)/2+1
        
        load([folder 'allresults_blueshift/jitter_' num2str(jitterval) '_noiselevel_' num2str(noiselevelarray_4(mmm)) '/results.mat']);

        
        subplot(311);
        imagesc(Psi_mod_noiselevel(mmm).Psi_mod(jj).I);
        axis image;
        colorbar;
        %title(['noise jj' num2str(jj)]);
        
        subplot(312);
        imagesc(data_exp(jj).I);
        axis image;
        colorbar;
        %title('no noise');
        
        subplot(313);
        imagesc((sqrt(Psi_mod_noiselevel(mmm).Psi_mod(jj).I)-sqrt(data_exp(jj).I)).^2);
        axis image;
        colorbar;
        
    end
    
    figure(1000);
    clf;
    
    for mmm = 1:numel(noiselevelarray_4)
       subplot(121);
       hold on;
       plot(log10(error_noiselevel(mmm).errlist),'LineWidth',3.0); 
       
       subplot(122);
       hold on;
       plot(log10(error_noiselevel(mmm).errlist*sqrt(mn_noiselevel(mmm).mn/cntrate_array_4(mmm))),'LineWidth',3.0); 

       
      legend_str{mmm} = ['noiselevel = ' num2str(noiselevelarray_4(mmm)) 'count rate = ' num2str(cntrate_array_4(mmm))]; 

      
      error_jitter(mmm) = error_noiselevel(mmm).errlist(end)*sqrt(mn_noiselevel(mmm).mn/cntrate_array_4(mmm));
      
    end
legend(legend_str);

