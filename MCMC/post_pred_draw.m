function new_traj = post_pred_draw(N,S,D,curr_theta,curr_Delta,curr_delta,block_index,strain_inc,FFT_stress,new_traj,k,store,figure_path)
    
    %The current parameters at iteration k are approximate samples from the
    %marginal posterior over theta*,and the precisions
    
    Kosher = 0;
           while Kosher ==0
                
               %sample a random effect from the likelihood given 
               %theta*, and the precisions
               rand_eff = mvnrnd(curr_theta(S*D + 1:end),curr_Delta\eye(D));
                          
               if rand_eff(4) >= rand_eff(3)
                   continue
               end

               if sum(rand_eff <= 0) ~= 0
                   continue
               end

               Kosher = 1;
           end
        
        %evaluate vpsc at sampled random effect
        rand_eff_vpsc = VPSC(rand_eff(1),rand_eff(2),rand_eff(3),rand_eff(4),strain_inc);
        
        %sample new trajectory from the likelihood given vpsc evaluated at
        %random effect and precisions
        new_traj(:,k/store) = mvnrnd(rand_eff_vpsc,eye(N)*(1/curr_delta));
        
        %'new' trajectory
        figure %plot the true virtural ground truth data (before sample variance introduced
        set(gcf,'Visible','off')
        h1 = plot(strain_inc,FFT_stress(:,:),'r--','Linewidth',1.2);
        hold on             
                
        for s = 1:k/store
            xflip = [strain_inc' fliplr(strain_inc')];
            yflip = [new_traj(:,s)' fliplr(new_traj(:,s)')];
            p = patch(xflip,yflip,'b','EdgeAlpha',.1,...
                'FaceColor','none','Linesmoothing','on','Linewidth',1);
        end
        title('Draws from Posterior Predictive Distribution')
        xlabel('strain')
        ylabel('stress (MPa)')
        
        print( fullfile(figure_path,'Draws from Post Pred'),'-dpng')
        %print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Draws from Posterior Predictive'],'-dpng')
		
		%New trajectory distribution
		figure
		set(gcf,'Visible','off')
		title('New Trajectory Distribution')
		width = 0.75;
        
        zoom = [2,50,98,150;...
                0,3,6,9]; % 0%, 3%, 6%, 9.2% strain
        
        for i = 1:length(zoom)
            subplot(2,2,i)
            h = histfit(new_traj(zoom(1,i),1:k/store),10, 'kernel');
            set(h(2),'color','m'); set(h(2),'linewidth',width)
            title([num2str(zoom(2,i)),'% Strain'])
            xlabel('stress (MPa)')
            ylabel('counts')
        end
        print( fullfile(figure_path,'Distribution of Draws from PPD'),'-dpng')
		%print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Distribtuion of Draws from PPD'],'-dpng')
	end