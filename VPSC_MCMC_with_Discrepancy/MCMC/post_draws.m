%Copyright (c) 2020, Denielle Ricciardi
%All rights reserved.                  

function [new_traj, post_draw] = post_draws(N,S,D,curr_theta,curr_Lambda,curr_delta,...
    curr_Delta,strain_inc,FFT_stress,new_traj,post_draw,k,store,figure_path)
    
    %The current parameters at iteration k are approximate samples from the
    %marginal posterior over theta, the random effects precision,
	%the observation error precision, and the discrepancy

%% Posterior predictive draw    
    Kosher = 0;
	while Kosher ==0                
		%sample a random effect from the likelihood given 
		%theta, and the random effect precision
		rand_eff = mvnrnd(curr_theta(S*D + 1:end),curr_Lambda\eye(D));

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
	%new random effect and error precision
	err_pr_mat = diag(repmat(curr_delta,1,N));
	new_traj(:,k/store) = mvnrnd(rand_eff_vpsc + curr_Delta,err_pr_mat\eye(N));        
	 
	figure 
	set(gcf,'Visible','off')
	h1 = plot(strain_inc,FFT_stress(:,:),'r--','Linewidth',.5);
	hold on             
			
	for s = 1:k/store
		xflip = [strain_inc' fliplr(strain_inc')];
		yflip = [new_traj(:,s)' fliplr(new_traj(:,s)')];
		p = patch(xflip,yflip,'b','EdgeAlpha',.2,...
		'FaceColor','none','Linesmoothing','on','Linewidth',1); %edge alpha originally .1
	end

	xlim([-.01 .091])
	ylim([0 200])
	title('Draws from Posterior Predictive Distribution')
	xlabel('strain')
	ylabel('stress (MPa)')
	
	print( fullfile(figure_path,'Draws from Post Pred'),'-dpng')
        
%% Posterior Draw
	figure  
	set(gcf,'Visible','off')
	h1 = plot(strain_inc,FFT_stress(:,:),'r--','Linewidth',.5);
	hold on             
	
	par = curr_theta(S*D+1:end);
	pdraw = VPSC(par(1),par(2),par(3),par(4),strain_inc);
  
	post_draw(:,k/store) = pdraw + curr_Delta;        
	for s = 1:k/store
		xflip = [strain_inc' fliplr(strain_inc')];
		yflip = [post_draw(:,s)' fliplr(post_draw(:,s)')];
		p = patch(xflip,yflip,'b','EdgeAlpha',.2,...
			'FaceColor','none','Linesmoothing','on','Linewidth',1); %edge alpha originally .1
	end
	xlim([-.01 .091])
	ylim([0 200])
	title('Draws from Posterior Distribution')
	xlabel('strain')
	ylabel('stress (MPa)')
	
	print( fullfile(figure_path,'Draws from Posterior'),'-dpng') 
end