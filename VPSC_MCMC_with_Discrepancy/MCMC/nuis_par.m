function nuis_par_dist = nuis_par(strain_inc,FFT_stress,k,store,S,yin,nuis_par_dist,figure_path)

	figure  
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
end