%Copyright (c) 2020, Denielle Ricciardi
%All rights reserved.                  

function [full_cov,e_trans,cov_1,cov_2] = gen_cov(sigma_1,sigma_2,w_1,w_2,nug_1,nug_2,N,...
	strain_inc,el_pl_trans,show_figs,analyze)
	
u = strain_inc;
[~,e_trans] = min(abs(u - el_pl_trans)); %find strain index corresonding to elastic-plastic transition
u_sub_1 = u(1:e_trans); %indices belonging to elastic region
u_sub_2 = u(e_trans + 1:end); %indices belonging to plastic region

%Specify the left limit of the time domain
a = 0;

nugget_1 = eye(e_trans)*nug_1;
nugget_2 = eye(length(u_sub_2))*nug_2;
template = zeros(N,N);

%Non-stationary covariance structure for elastic region
cov_1 = sigma_1.*QQ1d_se(fliplr(u_sub_1),fliplr(u_sub_1),w_1,a);

%Stationary squared exponential covariance structure for plastic region
cov_2 = sq_exp(u_sub_2,u_sub_2,w_2,sigma_2);   

if ~issymmetric(cov_1)
    cov_1_symmetric = nearestSPD(cov_1);
    cov_1 = cov_1_symmetric;
end

if ~issymmetric(cov_2)
    cov_2_symmetric = nearestSPD(cov_2);
    cov_2 = cov_2_symmetric;
end

template(1:e_trans,1:e_trans) = rot90(cov_1,2) + nugget_1;
template(e_trans+1:end,e_trans+1:end) = cov_2 + nugget_2;
cov_2 = cov_2 + nugget_2;
full_cov = template;

mat = cov_1;
eig_mat = eig(mat);
flag = 0;
for i = 1:rank(mat)
  if eig_mat(i) <= 0 
  flag = 1;
  end
end
if flag == 1
   disp('Elastic covariance matrix is not positive definite')
  else
   disp('Elastic covariance matrix is positive definite')
end

mat = cov_2;
eig_mat = eig(mat);
flag = 0;
for i = 1:rank(mat)
  if eig_mat(i) <= 0 
  flag = 1;
  end
end
if flag == 1
   disp('Plastic covariance matrix is not positive definite')
else
   disp('Plastic covariance matrix is positive definite')
end

mat = full_cov;
eig_mat = eig(mat);
flag = 0;
for i = 1:rank(mat)
  if eig_mat(i) <= 0 
  flag = 1;
  end
end
if flag == 1
   disp('Full covariance matrix is not positive definite')
  else
   disp('Full covariance matrix is positive definite')
end

if show_figs == 1
    mu = zeros(N,1);  
    nsims = 10;              
    Delta = zeros(N,nsims);  

    figure
    hold on
    for j = 1:nsims
        Delta(:,j) = mvnrnd(mu,full_cov);
        plot(u,Delta(:,j))
    end   
    xlim([-strain_inc(1)*1.1 strain_inc(end)*1.1])
    
    axis tight
    set(gca, 'FontSize',22,'FontWeight','bold','box','on');
    xlabel('strain','fontsize',24,'fontweight','bold')
    set(gcf, 'Position', [100 100 700 600])
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    filename = 'Delta_full_prior';
    print('-painters',filename,'-dpng')
    
    figure
    hold on
    for j = 1:nsims
        plot(u_sub_2,Delta(e_trans+1:end,j))
    end   
    xlim([-strain_inc(1)*1.1 strain_inc(end)*1.1])
   
    axis tight
    set(gca, 'FontSize',22,'FontWeight','bold','box','on');
    xlabel('strain','fontsize',24,'fontweight','bold')
    set(gcf, 'Position', [100 100 700 600])
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    filename = 'Delta_plastic_prior';
    print('-painters',filename,'-dpng')
    
    
    figure
    hold on 
    for j = 1:nsims
        plot(u_sub_1,Delta(1:e_trans,j))
    end
    
    axis tight
    set(gca, 'FontSize',22,'FontWeight','bold','box','on');
    xlabel('strain','fontsize',24,'fontweight','bold')
    set(gcf, 'Position', [100 100 700 600])
    set(gcf, 'Units', 'inches');
    screenposition = get(gcf,'Position');
    set(gcf,'PaperPosition',[0 0 screenposition(3:4)],...
        'PaperSize',[screenposition(3:4)]);
    filename = 'Delta_elastic_prior';
    print('-painters',filename,'-dpng')
    
end

if analyze == 1    
    disp(['************************************************* '])
    disp(['w 1 = ' num2str(w_1)])
    disp(['nugget 1 = ' num2str(nug_1)])
    disp(['signal variance 1 = ' num2str(sigma_1)]) 
    disp(['The condition of the elastic covariance matrix is ', num2str(cond(full_cov(1:e_trans,1:e_trans)))])   
    
    disp(['************************************************* '])
    disp(['w 2 = ' num2str(w_2)])
    disp(['nugget 2 = ' num2str(nug_2)])
    disp(['signal variance 2 = ' num2str(sigma_2)])
    disp(['The condition of the plastic covariance matrix is ', num2str(cond(full_cov(e_trans+1:end,e_trans+1:end)))])
    
    disp(['************************************************* '])
    disp(['The condition of the covariance matrix is ', num2str(cond(full_cov))])
end
cov_1 = full_cov(1:e_trans,1:e_trans);
cov_2 = full_cov(e_trans+1:end,e_trans+1:end);
end


