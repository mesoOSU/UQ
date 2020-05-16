function [ml_parameters] = diagnostics(par,Delta_chain,log_posterior,k,thin,zz,sigma_obs,strain_inc,final,iter_start,figure_path)

if zz <= 7
    npar = size(par,2);
    name = {'\tau_0','\tau_1','\theta_0','\theta_1',num2str(zz)};
    spr = npar/2; %number of rows in subplot
    spc = npar/2; %number of columns in subplot
end

if zz == 8
    npar = 16;

    name = {'1/\xi_{1}','\nu_{12}','\nu_{13}','\nu_{14}',...
            '\nu_{21}','1/\xi_{2}','\nu_{23}','\nu_{24}',...
            '\nu_{31}','\nu_{32}','1/\xi_{3}','\nu_{34}',...
            '\nu_{41}','\nu_{42}','\nu_{43}','1/\xi_{4}',...
            'Inverse Covariance Matrix'};
    spr = 4;
    spc = 4;
end

if zz == 9
    npar = size(par,2);
    name = {'Error Precision - tau'};
    spr = 1;
    spc = 1;
end


if iter_start + 250 > k
    iter_start = thin;
end

%map estimate
%location in parameter space which has maximum log_posterior value
[~,I] = max(log_posterior(1:k/thin));

ml_parameters = par(I,:);
ml_Delta = Delta_chain(:,I);

%% plot of model fit to data

if zz < 7 
    
    figure
    if final == 0; set(gcf,'Visible','off'); end
        Linethickness = 1;
        h1 = plot(strain_inc,sigma_obs(:,zz),'b','Linewidth',Linethickness);
        axis tight
        hold on
        map_out = VPSC(ml_parameters(1),ml_parameters(2),ml_parameters(3),ml_parameters(4),strain_inc);
        h2 = plot(strain_inc,map_out,'r','Linewidth',2);
        map_out_0 = VPSC(par(1,1),par(1,2),par(1,3),par(1,4),strain_inc);
        h3 = plot(strain_inc,map_out_0,'g--','Linewidth',Linethickness);
        current = VPSC(par(k/thin,1),par(k/thin,2),par(k/thin,3),par(k/thin,4),strain_inc);
        h4 = plot(strain_inc,current,'k--','Linewidth',Linethickness);
        h5 = plot(strain_inc,ml_Delta,'m','Linewidth',Linethickness);
        axis tight
        title(['MAP Parameter Result for data ',num2str(zz)])
        xlabel('Von Mises Strain')
        ylabel('Von Mises Stress (MPa)')
        legend([h1(1), h2(1), h3(1) h4(1) h5(1)],'Ground Truth Data','MAP Parameter Data','Initial Parameter Selection','Current Evaluation','MAP Delta','Location','best')
        xlim([-.001 .1])
        print(fullfile(figure_path,['Model_fit_',num2str(zz)]),'-dpng')
%        print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Model_fit_',num2str(zz)],'-dpng')
end

%% Delta plot

if zz == 1 
    
    figure
    if final == 0; set(gcf,'Visible','off'); end
        Linethickness = 1;
        hold on
        h2 = plot(strain_inc,ml_Delta,'r','Linewidth',2);
        h3 = plot(strain_inc,Delta_chain(:,1),'g--','Linewidth',Linethickness);
        h4 = plot(strain_inc,Delta_chain(:,k/thin),'k--','Linewidth',Linethickness);
        axis tight
        xlim([-.001 .1])
        title(['MAP Result for Delta '])
        xlabel('Von Mises Strain')
        ylabel('Discrepancy (\Delta)')
        legend([h2(1), h3(1) h4(1)],'MAP Delta','Initial Delta','Current Delta','Location','best')

        print(fullfile(figure_path,['Delta']),'-dpng')
end

%% posterior correlation plots
if zz <= 7
    figure
    if final == 0; set(gcf,'Visible','off'); end

    plotmatrix(par(iter_start/thin:k/thin,:))
    axis tight
    title(['Correlation for data ',num2str(zz)])
    
    print(fullfile(figure_path,['Post_corr_plot',num2str(zz)]),'-dpng')
    
    %print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Post_corr_plot_',num2str(zz)],'-dpng')
end

%% trace plots

figure
if final == 0; set(gcf,'Visible','off'); end
    
    if zz <= 7 || zz == 9
    for i = 1:npar
        subplot(spr,spc,i)
        plot(par(iter_start/thin:k/thin,i))
        axis tight
        xlabel('Iteration')
        ylabel('Parameter Value')
        title(name{i})
     end
     end
     
    if zz == 8
    for i = 1:npar
        subplot(spr,spc,i)
        
        for rr = 1:4
            for cc = 1:4
                plot(squeeze(par(rr,cc,iter_start/thin:k/thin)))
                axis tight
                xlabel('Iter')
                ylabel('Par Value')
                title(name{i})
            end
        end
     end
     end
    
    print(fullfile(figure_path,['Trace_plots',name{end}]),'-dpng')
    
    %print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Trace_plots_',name{end}],'-dpng')
    


%% Posterior marginal histograms
figure
if final == 0; set(gcf,'Visible','off'); end

width = 0.75;
    
    if zz <= 7 || zz == 9
        for i = 1:npar
            subplot(spr,spc,i)
            k1 = histfit(par(iter_start/thin:k/thin,i),25,'kernel');
            set(k1(2),'color','m'); set(k1(2),'linewidth',width)
            xlabel('Parameter Value')
            title(name{i})
        end
    end
    
    if zz == 8
    for i = 1:npar
        subplot(spr,spc,i)
        for rr = 1:4
            for cc = 1:4
                k1 = histfit(squeeze(par(rr,cc,iter_start/thin:k/thin)),25,'kernel');
                set(k1(2),'color','m'); set(k1(2),'linewidth',width)
                xlabel('Parameter Value')
                title(name{i})
            end
        end
     end
     end
     
    print(fullfile(figure_path,['Marginal_hist_plots_',name{end}]),'-dpng')
    
    %print(['/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Marginal_hist_plots_',name{end}],'-dpng')

%% log_posterior plot
if zz == 1 %only need to calculate 1x
    figure
    if final == 0; set(gcf,'Visible','off'); end

    plot(iter_start/thin:k/thin+1,log_posterior(iter_start/thin:k/thin+1))
    xlabel('Iteration')
    ylabel('Log Posterior')
    title('Log Posterior')

    print(fullfile(figure_path,'Log_post_plot'),'-dpng')
    
    %print('/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/Log_post_plot','-dpng')
end
    