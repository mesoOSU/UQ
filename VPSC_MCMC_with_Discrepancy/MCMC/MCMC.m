%Copyright (c) 2020, Denielle Ricciardi
%All rights reserved.                  

%% FOLDER SETUP
[filepath] = fileparts(which('MCMC.m'));
parts = strsplit(filepath, '\');
path = strjoin(parts(1:end-1),'/');

addpath(sprintf('%s/MCMC',path))
addpath(sprintf('%s/vpsc7d_virgin',path))

figure_path = strcat(strjoin(parts(1:end-1),'/'),'/figures');

PATH = getenv('PATH');

setenv('PATH', [PATH sprintf(':%s/MCMC',path),...
                     sprintf(':%s/vpsc7d_virgin',path)]);
                 

%Navigate to correct path to run VPSC code
cd(sprintf('%s//vpsc7d_virgin',path)) 

%% USER DEFENITIONS

%If initial start = 0, chain is initiated from current state of previously saved chain
initial_start = 0; 

%Number of samples to draw in Markov chain
number_iterations = 100000;

%Diagnostics Settings
%Iter_start is the number of samples to be counted as burn-in
iter_start = 1; %Must be a multiple of thin
iter_check = 1000; %Frequency of printing diagnostic plots
num_par = 4; %Number of unknown parameters in VPSC model
thin = 1; %Factor chain is thinned by (should be a factor of num_iter and iter_check) 

%Load data
FFT_data = load(sprintf('%s//FFT_Simulated_Data.mat',path));  
n_data = size(FFT_data.SVM,2); %Number of experiments
store_ppd = 500; %Frequency of storing posterior predictive evaluations  

%Initialize proposal covariance matrix
%One covariance matrix for each block of parameters sampled in a MH step
cov_mat = cell(1,n_data + 2); 

for i = 1:size(cov_mat,2)
    if initial_start == 1 
		%Choose reasonable initial covariance
        cov_mat{i} = eye(num_par)*10;
    else
		import_cov_mat = load(sprintf('%s/MCMC/current_cov_mat.mat',path));
        cov_mat{i} = import_cov_mat.cov_mat{i};
    end    
end     
                 
start_fix = number_iterations;  %Iteration at which to start adaptation of proposal covariance  
stop_fix =  number_iterations; %Iteration at which to end adaptive proposal covariance (only adapt during burn-in)
covariance_check = 2000; %Frequency of adapting porposal variance based on acceptance rate
target_accept = [0.2,0.5]; %Target acceptance rates

v = 1:n_data+2; 
blocks = repelem(v,num_par);%define blocks

%% CREATE BLOCKS
theta_chain_index = 1:1:length(blocks); 
u = unique(blocks);
num_blocks = length(u);
block_index = cell(1,num_blocks);
ml_par = cell(1,num_blocks);

for nn = 1:num_blocks
    block_index{nn} = theta_chain_index(blocks == u(nn));
end

%% LOAD DATA AND SET CODE DEFENITIONS
%Load data
sigma_obs = zeros(size(FFT_data.SVM,1),n_data);

S = size(sigma_obs,2); %Number of random effects
D = 4; %Number of unknown physical model parameters
N = length(sigma_obs); %Number of observations per experiment or random effect

for i = 1:S
    sigma_obs(:,i) = FFT_data.SVM(:,i);
end

strain_inc = FFT_data.EVM;

%Symmetric proposal Distribution
proposal_sample = @(mu,proposal_variance) mvnrnd(mu,proposal_variance);

%Initialize arrays
num_iter = number_iterations; %Number of iterations
theta_chain = zeros(num_iter/thin+1,block_index{S+1}(D)); %Chain of accepted values for theta^1,...,theta^S,theta,t^2
Lambda_chain = zeros(D,D,num_iter/thin+1); %Chain holding random effects precision matrix samples
delta_chain = zeros(num_iter/thin+1,1); %Chain holding error precision samples
Delta_chain = zeros(N,num_iter/thin+1); %Chain holding discrepancy samples

%Parameters for Lambda decomposition
t2_chain = zeros(D,num_iter/thin+1); %Diagonal variances of random effects precision matrix (Lamda)
R_chain = zeros(D,D,num_iter/thin+1); %Off-diagonal correlations of random effects precision matrix (Lambda)

%Array holding 'covariance_check' number of samples for full adjustment 
temp_theta_chain = zeros(covariance_check,block_index{S+2}(D)); 

log_post = zeros(num_iter/thin+1,1); %Array for the log posterior value at each iteration
acceptance = zeros(num_iter,S+2); %Array to track parameter acceptance, last two blocks are theta and t^2
nuis_par_dist = zeros(N,num_iter/store_ppd*S); %Array for nuisance parameter distribution
ppd_draw = zeros(N,num_iter/store_ppd); %Array for posterior predictive distrirbution    
post_draw = zeros(N,num_iter/store_ppd); %Array for VPSC model evaluated at parameter posterior draws

%% INITIALIZE MARKOV CHAIN

%Hyperparameters for R prior (Wishart)
r_not = D + 2; %scalar 
R_not = 1/r_not*eye(D); %DxD positive semidefinite matrix

%Hyperparameters for t (Gamma)
a_t = 1;%1; %SHAPE
b_t = 0.1; %RATE

%Hyperparameters for error precision prior (Gamma)  
a_delta = 1;  %SHAPE parameter
b_delta = .1; %RATE parameter - MATLAB uses a SCALE parameter so inverse must be taken for built-in functions

%Initiation of Markov Chain        
if initial_start == 1
    theta_chain(1,:) = repmat([55 55 700 150],1,S+1);
    R_chain(:,:,1) = r_not.*R_not;
    t2_chain(:,1) = ones(D,1)*.1;
    Lambda_chain(:,:,1) = diag(sqrt(t2_chain(:,1)))*R_chain(:,:,1)*diag(sqrt(t2_chain(:,1)));
    delta_chain(1,:) = 10;
    Delta_chain(:,1) = diag(ones(N));    
else
	import_last_parameters = load(sprintf('%s/MCMC/current_parameters.mat',path));
    theta_chain(1,:) = import_last_parameters.current_parameters{1};
    Lambda_chain(:,:,1) = import_last_parameters.current_parameters{2};
    delta_chain(1,1) = import_last_parameters.current_parameters{3};
    Delta_chain(:,1) = import_last_parameters.current_parameters{4};
    R_chain(:,:,1) = import_last_parameters.current_parameters{5};
    t2_chain(:,1) = import_last_parameters.current_parameters{6};
end

current_theta = theta_chain(1,:);
current_Lambda = Lambda_chain(:,:,1);
current_delta = delta_chain(1,1);
current_t2 = t2_chain(:,1);
current_R = R_chain(:,:,1);
current_Delta = Delta_chain(:,1);

temp_chain(1,1:block_index{S+1}(D)) = theta_chain(1,:);
temp_chain(1,block_index{S+2}) = t2_chain(:,1)';

%% Initialize GP
if initial_start ~= 0
    %Initialize and define covariance structure for Gaussian Process, which is used as a prior for the discrepancy 
    discretization = diff(strain_inc);

    %Length scale
    w_1 = .35*discretization(1); %Elastic region
    w_2 = 50*discretization(2);  %Plastic region

    nug_1 = .1; %Elastic region
    nug_2 = 1E-9; %Plastic region

    signal_variance_1 =  1E14; %Elastic region
    signal_variance_2 = 1E-4; %Plastic region

    el_pl_trans = .0025; %Elastic-plastic transition

    %Flags
    show_figs = 0;
    analyze = 0;

    [Gamma,e_trans,cov_1,cov_2] = gen_cov(signal_variance_1,signal_variance_2,...
        w_1,w_2,nug_1,nug_2,N,strain_inc,el_pl_trans,show_figs,analyze);

    log_det_Gamma = double(vpa(log(det(vpa(Gamma,1000))),1000));
    
%% Calculate the inverse of Gamma
    digitsOld = digits;
    digits(16)
    syms A
    A = vpa(Gamma); 

    disp(['Condition of covariance matrix is ',num2str(cond(double(A)))])

    inv_Gamma = vpa(inv(Gamma));    

    if issymmetric(double(inv_Gamma))
        disp('Inverse is symmetric')
    elseif ~issymmetric(double(inv_Gamma))
        disp('Inverse is not symmetric')
    end

    if ~isposdef(double(inv_Gamma))
       disp('Inverse is not positive definite')
    elseif isposdef(double(inv_Gamma))
       disp('Inverse is positive definite')
    end
    
    %Check quality of inverse
    disp(['|| I - inv(Gamma)*Gamma || = ',num2str(double(norm(eye(N) - inv_Gamma*Gamma,1)))])
    Gamma = vpa(Gamma);  
    save Gamma Gamma inv_Gamma log_det_Gamma 
    
elseif initial_start == 0
    
    load_Gamma = load('Gamma');
    Gamma = load_Gamma.Gamma;
    inv_Gamma = load_Gamma.inv_Gamma;
    log_det_Gamma = load_Gamma.log_det_Gamma;
end

dubIG = double(inv_Gamma);
dubG = double(Gamma);

sigma_vpsc = zeros(N,S); %Initiate array for holding VPSC model evaluations, used in log likelihood calculation 
for i = 1:S    
    par = theta_chain(1,block_index{i});
    sigma_vpsc(:,i) = VPSC(par(1),par(2),par(3),par(4),strain_inc);     
end

%Initial log likelihood
[initial_log_likelihood,vpsc_prop] = log_likelihood(theta_chain(1,:),Lambda_chain(:,:,1),...
    delta_chain(1),Delta_chain(:,1),sigma_obs,strain_inc,sigma_vpsc,1,block_index); 
    
sigma_vpsc = vpsc_prop;

%Initial log posterior
initial_log_posterior = log_posterior(S,N,D,block_index,initial_log_likelihood,theta_chain(1,:),...
    t2_chain(:,1),R_chain(:,:,1),delta_chain(1),Delta_chain(:,1),log_det_Gamma,dubIG,a_delta, ...
	b_delta,R_not,r_not,a_t,b_t);

log_post(1,1) = initial_log_posterior;
current_log_post = log_post(1,1);

%% MCMC LOOP
disp('****************************')
disp('MCMC simulation initiated')
disp('****************************')

tic 
for k = 1:num_iter     
     if mod(k,covariance_check) == 0         
         if k <= stop_fix && k >= start_fix         
             accept_rate = (1 + sum(acceptance(k-covariance_check + 1:k,:)))./covariance_check;                
             for aa = 1:num_blocks            
                if accept_rate(aa) < target_accept(1) || accept_rate(aa) > target_accept(2)                   
                   C = cov(temp_chain(:,block_index{aa}));
                   cov_mat{aa}(logical(eye(D))) = diag(C).*accept_rate(aa)./(target_accept(1) + range(target_accept)/2);
                   cov_mat{aa} = full_cov_mat_adj(temp_chain(:,blocks == aa),D,cov_mat{aa});                    
                end   
             end 
         end 
         temp_chain = zeros(covariance_check,block_index{S+2}(D)); %refresh for next full-adjustment
     end
        
    proposed_theta = current_theta;    
    for b = 1:length(unique(blocks)) - 1 
         
        %Check that proposed Voce parameters meet 'kosher' conditions
        Kosher = 0;        
        while Kosher ==0            
            %Draw proposed parameters from the symmetric candidate distribution
            %Centered at the current parameter
            proposed_theta(block_index{b}) = proposal_sample(current_theta(block_index{b}),cov_mat{b});

            %Check that voce parameters meet 'Kosher' requriements            
            if proposed_theta(b*D) >= proposed_theta(b*D-1)
                continue
            end            
            if proposed_theta(b*D) < 0 || proposed_theta(b*D-3) < 0 || proposed_theta(b*D-2) < 0
                continue
            end            
            Kosher = 1;
        end        
        
        try            
		u = unifrnd(0,1);

		%Proposed log likelihood
		[proposed_log_likelihood,vpsc_prop] = log_likelihood(proposed_theta,current_Lambda,...
		current_delta,current_Delta,sigma_obs,strain_inc,sigma_vpsc,b,block_index); 

		%Proposed log posterior
		proposed_log_posterior = log_posterior(S,N,D,block_index,proposed_log_likelihood,proposed_theta,...
		current_t2,current_R,current_delta,current_Delta,log_det_Gamma,dubIG,a_delta,b_delta,R_not,r_not,a_t,b_t);    

		if log(u) < proposed_log_posterior-current_log_post
			current_theta(block_index{b}) = proposed_theta(block_index{b}); 
			current_log_post = proposed_log_posterior;                
			acceptance(k,b) = 1;
			sigma_vpsc = vpsc_prop;                         
		else
			acceptance(k,b) = 0;
			proposed_theta(block_index{b}) = current_theta(block_index{b}); 
		end

		catch
			sprintf('Error in MH step iteration (%i,%i)',k,b)
			kk = k+1;
			sprintf('\n Continuing to next iteration with parameter(%i,%i) = parameter(%i,%i)',kk,b,k,b)
			proposed_theta(block_index{b}) = current_theta(block_index{b}); 
			acceptance(k,b) = 0;            
			continue            
        end        
    end  
    
      
%% Sample t^2 in MH step
             
    b = S+2;
    pos_t2 = 0;
    while pos_t2 == 0
        proposed_t2 = proposal_sample(current_t2,cov_mat{b});
        if proposed_t2 > 0
            pos_t2 = 1;
        end
    end
    
    proposed_lambda = diag(sqrt(proposed_t2))*current_R*diag(sqrt(proposed_t2));
    u = unifrnd(0,1);      
    try
    [proposed_log_likelihood,~] = log_likelihood(current_theta,proposed_lambda,...
        current_delta,current_Delta,sigma_obs,strain_inc,sigma_vpsc,b,block_index);

    proposed_log_posterior = log_posterior(S,N,D,block_index,proposed_log_likelihood,current_theta,...
        proposed_t2,current_R,current_delta,current_Delta,log_det_Gamma,dubIG,a_delta,b_delta,R_not,r_not,a_t,b_t);

    if log(u) < proposed_log_posterior-current_log_post
        current_t2 = proposed_t2;
        current_log_post = proposed_log_posterior;                
        acceptance(k,b) = 1;               
    else       
		acceptance(k,b) = 0;
    end

    catch
        sprintf('Error in MH step iteration (%i,%i)',k,b)
        kk = k+1;
        sprintf('\n Continuing to next iteration with parameter(%i,%i) = parameter(%i,%i)',kk,b,k,b)    
        acceptance(k,b) = 0;        
    end
	
 %% Gibbs updates 
 
 %Parameters Psi (delta), R, and Delta are sampled directly from their full-condition posterior distributions
 %Always accepted
 
 %% Sample R
    try
    sq_diff_pars  = 0;
    for i = 1:S
        sq_diff_pars = sq_diff_pars + diag(sqrt(current_t2))*(current_theta(block_index{i}) ...
        - current_theta(block_index{S+1}))'*(current_theta(block_index{i})...
        - current_theta(block_index{S+1}))*diag(sqrt(current_t2));
    end     
    current_R = wishrnd((R_not\eye(D) + sq_diff_pars)\eye(D),r_not + S);         
    catch
        %disp(['Numerical error in Gibbs step for R at iteration',num2str(k)])        
    end    
    current_Lambda = diag(sqrt(current_t2))*current_R*diag(sqrt(current_t2));
    
%% Sample Psi (delta)- the error precision       
                          
	sq_diff = (sigma_obs - sigma_vpsc - current_Delta).^2;         
	current_delta = gamrnd(S*N/2 + a_delta,...
		1/(b_delta + 1/2*(sum(sum(sq_diff))))); 
        
%% Sample Delta - the discrepancy 

	err_pr_mat = diag(repmat(current_delta,1,N));
	sample_mean = mean((sigma_obs - sigma_vpsc),2);
	fc_var = (dubIG + N.*err_pr_mat)\eye(N);
	fc_mean = fc_var*(N.*err_pr_mat)*sample_mean;          
	current_Delta = mvnrnd(fc_mean',fc_var)';        
        
%% Calculate log likelihood after gibbs step
	[Gibbs_log_likelihood,vpsc_prop] = log_likelihood(current_theta,...
	current_Lambda,current_delta,current_Delta,sigma_obs,strain_inc,sigma_vpsc,b,block_index); 

	Gibbs_log_posterior = log_posterior(S,N,D,block_index,Gibbs_log_likelihood,...
	current_theta,current_t2,current_R,current_delta,current_Delta,...
	log_det_Gamma,dubIG,a_delta,b_delta,R_not,r_not,a_t,b_t);

	current_log_post = Gibbs_log_posterior;        

    %Save parameters by increment defined by 'thin'
	if mod(k,thin) == 0
		theta_chain(k/thin + 1,:) = current_theta;
		t2_chain(:,k/thin+1) = current_t2;
		R_chain(:,:,k/thin+1) = current_R;
		Lambda_chain(:,:,k/thin + 1) = current_Lambda;
		delta_chain(k/thin + 1,1) = current_delta;
		Delta_chain(:,k/thin+1) = current_Delta;
		log_post(k/thin+1,1) = current_log_post; 
	end
	
    %After each iteration save current state of all parameters sampled in 
    %MH steps in temporary array to be used in covariance adjustment        
	temp_index = rem(k,covariance_check);
	temp_chain(temp_index+1,1:block_index{S+1}(D)) = current_theta; 
	temp_chain(temp_index+1,block_index{S+2}) = current_t2';  
	  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%The rest of the script is for generating diagnostic plots%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
	if mod(k,store_ppd) == 0
		%Generate nuisance parameter distribution
		nuis_par_dist = nuis_par(strain_inc,sigma_obs,k,store_ppd,S,sigma_vpsc,nuis_par_dist,figure_path);        

		%Genereate Posterior Predictive Distribution
		try
		[ppd_draw, post_draw] = post_draws(N,S,D,current_theta,current_Lambda,current_delta,...
		current_Delta,strain_inc,sigma_obs,ppd_draw,post_draw,k,store_ppd,figure_path);      
		catch
		end   
	end

%% Generate diagnostics    
    if mod(k,iter_check) == 0
   
        toc;
        try 
         
        sum_acceptance = sum(acceptance(k-iter_check+1:k,:));
        acceptance_ratio = sum_acceptance/iter_check;        
        str = ["-----------------------------------------";...
        sprintf("Number of iterations so far = %i/%i",k,num_iter);...
        sprintf("Time for last %i iterations = %.2f minutes",iter_check,toc/60);...
        sprintf("Time of last diagnostic = %s", datestr(now, 'dd/mm/yy HH:MM'));...
        sprintf("Acceptance rates in past %i iterations:",iter_check);...
        sprintf("%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f",acceptance_ratio)];    
        fid = fopen(fullfile(figure_path,'stats.txt'), 'a+'); 
        fprintf(fid, '%s\n', str);        
        fclose(fid);
        
        disp('-----------------------------------------')
        disp(['Number of iterations so far = ',num2str(k),'/',num2str(num_iter)])
        disp(['Time for last ',num2str(iter_check),' iterations = ',num2str(toc/60),' minutes'])
        disp(['Time of last diagnostic = ', datestr(now, 'dd/mm/yy HH:MM')])
        sum_acceptance = sum(acceptance(k-iter_check+1:k,:));
        acceptance_ratio = sum_acceptance/iter_check;
        disp(['Acceptance rates in past ',num2str(iter_check),' iterations:'])
        disp(num2str(acceptance_ratio))        
         
        for zz = 1:num_blocks - 1
			[ml_voce] = diagnostics(theta_chain(1:k/thin + 1,block_index{zz}),...
				Delta_chain,log_post,k,thin,zz,sigma_obs,strain_inc,0,iter_start,figure_path);
			ml_par{zz} = ml_voce;
        end
        
        try
            stats = marg_hist_fig(D,S,theta_chain,k,iter_start,thin,figure_path);
        catch
        end
        
        figure
        set(gcf,'Visible','off')
        hold on
        for s = iter_start/thin:100:k/thin
            xflip = [strain_inc' fliplr(strain_inc')];
            yflip = [Delta_chain(:,s)' fliplr(Delta_chain(:,s)')];
            patch(xflip,yflip,'b','EdgeAlpha',.1,'FaceColor','none','Linesmoothing','on','Linewidth',1)
        end
        title('Delta Posterior')
        xlabel('strain')
        print( fullfile(figure_path,'Delta Posterior'),'-dpng')        
        
        trace_figs(S,D,iter_start,k,thin,theta_chain,log_post,block_index,figure_path);
            
        disp('Diagnostic plots have been created in working directory')        
        save(fullfile(figure_path,'Workspace.mat'))
        save(fullfile(figure_path,'current_cov_mat.mat'), 'cov_mat')
        
        current_parameters{1} = current_theta;
        current_parameters{2} = current_Lambda;
        current_parameters{3} = current_delta;
        current_parameters{4} = current_Delta;
        current_parameters{5} = current_R;
        current_parameters{6} = current_t2;
        
        save( fullfile(figure_path,'current_parameters.mat'),'current_parameters')   
                 
        catch 
			sprintf('Invalid file identifier in iteration (%i,%i), no diagnostics created. Continuing to next iteration',k,b)
			save(fullfile(figure_path,'100K.mat'))
			save(fullfile(figure_path,'current_cov_mat.mat'), 'cov_mat')

			current_parameters{1} = current_theta;
			current_parameters{2} = current_Lambda;
			current_parameters{3} = current_delta;
			current_parameters{4} = current_Delta;
			current_parameters{5} = current_R;
			current_parameters{6} = current_t2;			
			save( fullfile(figure_path,'current_parameters.mat'),'current_parameters') 
        end     
        tic
    end
    close all
end