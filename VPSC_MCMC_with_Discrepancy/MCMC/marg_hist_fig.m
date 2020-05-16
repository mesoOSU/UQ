%% Margninal Histogram Figure

function stats = marg_hist_fig(D,S,theta_chain,k,iter_start,thin,figure_path)

stats = zeros(2,S+1);
samples = k/thin;

for d = 1:D
       
	figure
	set(gcf,'Visible','off');
	hold on
	par_index = d:D:(S+1)*D;      
	width = 1.5;

	for s = 1:S+1        
		if s == S + 1
			width = 6;
		end            

		pd = fitdist(theta_chain(iter_start/thin:samples,par_index(s)),'normal');
		stats(:,s) = [pd.mu,pd.sigma];           
		x_values = linspace(pd.mu-4*pd.sigma,pd.mu+4*pd.sigma,10000);
		y = pdf(pd,x_values);       

		if s <=3
			h = plot(x_values,y/max(y),'-o','MarkerIndices',1:200+s:length(y),'MarkerSize',10,'LineWidth',width);
		end

		if 3 < s && s<=6
			h = plot(x_values,y/max(y),'-*','MarkerIndices',1:200+s:length(y),'MarkerSize',10,'LineWidth',width);
		end

		if s > 6
			h = plot(x_values,y/max(y),'LineWidth',width);
		end
		hold on
	end
       
	axis tight 
	xlabel(['\theta_{',num2str(d),'}'],'fontsize',20)

	lgd = legend('\theta^{1}','\theta^{2}','\theta^{3}',...
	'\theta^{4}','\theta^{5}','\theta^{6}',...
	'\theta^{*}');
	set(lgd,'FontSize',18,'FontWeight','bold');
	set(gcf, 'Position',  [100, 100, 1700, 600])

	set(gcf,'Units','inches');
	screenposition = get(gcf,'Position');
	set(gcf,...
		'PaperPosition',[0 0 screenposition(3:4)],...
		'PaperSize',screenposition(3:4));        


	set(gcf,'Units','inches');
	screenposition = get(gcf,'Position');
	set(gcf,...
		'PaperPosition',[0 0 screenposition(3:4)],...
		'PaperSize',screenposition(3:4));

	print(fullfile(figure_path,['Marginal_Histograms_theta_',num2str(d)]),'-dpng','-painters')
end

close all
        
end