function trace_figs(S,D,iter_start,k,thin,theta_chain,log_post,block_index,figure_path)

    figure
    set(gcf,'visible','off');
    
    subplot(D+1,1,1)
    plot(log_post(iter_start/thin:k/thin),'r')

    for d = 1:D
        subplot(D+1,1,d+1)
        hold on
        for s = 1:S
            plot(theta_chain(iter_start/thin:k/thin,block_index{s}(d)),'b')
        end

        plot(theta_chain(iter_start/thin:k/thin,block_index{S+1}(d)),'g')
        ylabel(num2str(D))
        xlabel('sample')
    end
    
    
print(fullfile(figure_path,['Traces']),'-dpng')
end