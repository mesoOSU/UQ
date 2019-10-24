function nuis_par_dist = nuis_par(strain_inc,FFT_stress,k,store,S,yin,nuis_par_dist,figure_path)


        figure %plot the true virtural ground truth data (before sample variance introduced
        set(gcf,'Visible','off')
        plot(strain_inc,FFT_stress(:,:),'r--','Linewidth',1);
        hold on        
       
        nuis_par_dist(:,k/store*(S)-S + 1:k/store*(S)) = yin;
        
        for s = 1:k/store*S
            xflip = [strain_inc' fliplr(strain_inc')];
            yflip = [nuis_par_dist(:,s)' fliplr(nuis_par_dist(:,s)')];
            p = patch(xflip,yflip,'b','EdgeAlpha',.1,...
                'FaceColor','none','Linesmoothing','on','Linewidth',.5);
        end
        title('Nuisance Parameter Distribution')
        xlabel('strain')
        ylabel('stress (MPa)')
        
        print( fullfile(figure_path,'Nuisance_parameter_distribution'),'-dpng')
        
        %print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Nuisance_Parameter_Distribution'],'-dpng')
end