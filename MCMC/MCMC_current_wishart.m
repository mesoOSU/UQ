%This is the same as MCMC_random_variance except the Voce parameters
%are not rounded and no variance limit is set
%parameters are initiated at past MAP parameter values

% ADAPTIVE
% BLOCKS                                     

%% FOLDER SETUP
addpath(pwd)
addpath(':/users/PAS1064/osu8614/Bayesian Inference/VPSC_MCMC_Gamma_06/Psi_0_6RandEff/MCMC')

PATH = getenv('PATH');
setenv('PATH', [PATH ':/users/PAS1064/osu8614/Bayesian Inference/VPSC_MCMC/VPSC_MCMC_Gamma_06/Psi_0_6RandEff/vpsc7d_virgin',...
                     ':/users/PAS1064/osu8614/Bayesian Inference/VPSC_MCMC/VPSC_MCMC_Gamma_06/Psi_0_6RandEff/MCMC']);

[filepath] = fileparts(which('MCMC_current_wishart.m'));
parts = strsplit(filepath, '/');
path = strjoin(parts(1:end-1),'/');

figure_path = strcat(strjoin(parts(1:end-1),'/'),'/figures');

%Navigate to correct path to run VPSC code
cd(sprintf('%s//vpsc7d_virgin',path)) 

%% USER DEFENITIONS

number_iterations = 100000;

%Diagnostics
iter_start = 1; 
iter_check = 1000; 

%load data

FFT_data = load(sprintf('%s//FFT_Plastic_Data_06.mat',path)); %150x6 array
n_data = size(FFT_data.SVM,2);

store = 500; %frequency of storing posterior predictive evaluations 
                 %for uncertainty plot


import_cov_mat = load(sprintf('%s/current_cov_mat.mat',path));
cov_mat = cell(1,n_data + 1); %one block per data set, 1 proposal cov matrix per block 

 for i = 1:size(cov_mat,2)
     cov_mat{i} = import_cov_mat.cov_mat{i};
 end    
                 
start_fix = number_iterations;%100000;  %start adaptive proposal variance at this iteration
stop_fix =  number_iterations;%100000; %end adaptive proposal variance at this iteration
covariance_check = 1000;%100000; %frequency of adapting porposal variance based on acceptance rate
full_adj = 1000;%100000; %frequency of adjusting prop_var with cov ratio
%compare = zeros(length(proposal_variance),4); %store information on how proposal variance ratio is adjusted
target_accept = [0.2,0.5];

blocks = [1,1,1,1,...
          2,2,2,2,...
          3,3,3,3,...
          4,4,4,4,...
          5,5,5,5,...
          6,6,6,6,...
          7,7,7,7]; %define blocks
      
%% CREATE BLOCKS
 
% parameter_index = [1,2,3,4,5,6,7,8]; %tau0, tau1, theta0, theta1, Theta1, Theta2, Theta3, Theta4
theta_chain_index = 1:1:(n_data + 1)*4; %7*voce
u = unique(blocks);
num_blocks = length(u);
block_index = cell(1,num_blocks);
ml_par = cell(1,num_blocks);
% variance_index = cell(1,num_blocks);
for nn = 1:num_blocks
    block_index{nn} = theta_chain_index(blocks == u(nn));
%     variance_index{nn} = proposal_variance(blocks == u(nn)); 
end

%% LOAD DATA AND SET CODE DEFENITIONS

%Load data
sigma_obs = zeros(size(FFT_data.SVM,1),n_data);

S = size(sigma_obs,2);
D = 4;
N = length(sigma_obs);

for i = 1:S
    sigma_obs(:,i) = FFT_data.SVM(:,i);
end

strain_inc = FFT_data.EVM;

%define hyper-parameters for precision gamma prior distributions 
a_tau = 1;      %.01; %2.5; %SHAPE parameter
b_tau = 0.1;   %.01; %50; %RATE parameter - MATLAB uses a SCALE parameter so inverse must be taken
a_xi = 1;     %.01;%.01; %1.1; %SHAPE parameter
b_xi = 0.1;       %.01; %7; %RATE parameter - MATLAB uses a SCALE parameter so inverse must be taken

%Proposal Distribution
proposal_sample = @(mu,proposal_variance) mvnrnd(mu,proposal_variance);

% Initialize vectors
num_iter = number_iterations; %number of iterations

theta_chain = zeros(num_iter+1,(S + 1)*D); %chain holding accepted values for theta^1,...,theta^S,theta*
Delta_chain = zeros(D,D,num_iter+1); %chain holding inverse covariance matrix
delta_chain = zeros(num_iter+1,1); %chain holding error precision

log_post = zeros(num_iter+1,1); %initiate array for the log posterior value at each iteration
acceptance = zeros(num_iter,S+1); %initiate array for parameter acceptance, last two blocks are not accepted/rejected
nuis_par_dist = zeros(N,num_iter/store*S); %initiate array for nuisance parameter distribution
ppd_draw = zeros(N,num_iter/store); %initiate array for model evaluations from posterior predictive distrirbution     

%% INITIALIZE MARKOV CHAIN

% initiation of Markov Chain
import_last_parameters = load(sprintf('%s/current_parameters.mat',path));

theta_chain(1,:) = import_last_parameters.current_parameters{1};
%theta_chain(1,S*D + 1:end) = [55 55 700 220];

Delta_chain(:,:,1) = import_last_parameters.current_parameters{2};

delta_chain(1,1) = import_last_parameters.current_parameters{3};

%scalar
v_not = D;

%V_not = import_last_parameters.current_parameters{4};
V_not = eye(D);

sigma_vpsc = zeros(N,S); %initiate array for holding VPSC_stress, used in log likelihood evaluation 

tic
for i = 1:S
    sigma_vpsc(:,i) = VPSC(theta_chain(1,i*D-3),theta_chain(1,i*D-2),theta_chain(1,i*D-1),theta_chain(1,i*D),strain_inc);
end
toc

[initial_log_likelihood,vpsc_prop] = log_likelihood(theta_chain(1,:),Delta_chain(:,:,1),delta_chain(1),sigma_obs,strain_inc,sigma_vpsc,1); %likelihood of initiation parameters, 1is for block 1
sigma_vpsc = vpsc_prop;

initial_log_posterior = log_posterior(initial_log_likelihood,theta_chain(1,:),Delta_chain(:,:,1),delta_chain(1),a_tau,b_tau,V_not,v_not,D);

log_post(1,1) = initial_log_posterior;
previous_log_posterior = log_post(1,1);

%% Run this section to start chain from a previous iteration
% k =91325; %iteration to restart chain
% for i = 1:size(yin,2)
%     sigma_vpsc(:,i) = VPSC(theta_chain(k,i*4-3),theta_chain(k,i*4-2),...
%        theta_chain(k,i*4-1),theta_chain(k,i*4),strain_inc);
% end
% 
% % 
% logposterior_previous = log_post(k,1);
%
% proposal_variance = pvar_store_temp(floor(k/covariance_check)+1,:);
%%% proposal_variance = pvar_store(floor(k/covariance_check)+1,:);
%% MCMC LOOP
   
tic
%ktic = tic;
for k = 1:num_iter

   % if mod(k,100) == 0
%        disp(['iteration = ',num2str(k),'  time for last 100 = ',num2str(toc(ktic))])
%        ktic = tic;
%    end
     
     if mod(k,covariance_check) == 0
         
%          pvar_store(k/covariance_check+1,:) = proposal_variance;
%            pvar_store_temp(k/covariance_check + 1,:) = cov_mat;
         if k <= stop_fix && k >= start_fix
         
             accept_rate = (1 + sum(acceptance(k-covariance_check + 1:k,:)))./covariance_check;
                
             for aa = 1:num_blocks          
            
                if accept_rate(aa) < target_accept(1) || accept_rate(aa) > target_accept(2)
                    
                    cov_mat{aa}(logical(eye(D))) = cov_mat{aa}(logical(eye(D))).*accept_rate(aa)./(target_accept(1) + range(target_accept)/2);
                    cov_mat_adj = full_cov_mat_adj(theta_chain(k-covariance_check+1:k,blocks == aa),D,cov_mat{aa});
                    cov_mat{aa} = cov_mat_adj;   
                    
                end   
             end 
         end                       
     end
     
    proposed_parameter = theta_chain(k,:);    
    
    for b = 1:length(unique(blocks))
      
        % the following constraints only apply to voce parameters
            
         Kosher = 0;
        
         while Kosher ==0
            
            %draw proposed parameters from the symmetric candidate distribution
            %centered at the current parameter
            proposed_parameter(block_index{b}) = proposal_sample(theta_chain(k,block_index{b}),cov_mat{b});
            
            %check that voce parameters meet 'Kosher' requriements
            
            if proposed_parameter(b*D) >= proposed_parameter(b*D-1)
                continue
            end
            
            if proposed_parameter(b*D) < 0 || proposed_parameter(b*D-3) < 0 || proposed_parameter(b*D-2) < 0
                continue
            end
            
            Kosher = 1;
         end      

        try
            
            u = unifrnd(0,1);
   
            [proposed_log_likelihood,vpsc_prop] = log_likelihood(proposed_parameter,Delta_chain(:,:,k),delta_chain(k),sigma_obs,strain_inc,sigma_vpsc,b); %likelihood of proposed (k)
            %yout updates the column of yin which corresponds to the
            %current block. VPSC is evaluated at proposed parameters and saved in yin(block)
            %if the parameters are accepted yin = yout, if not, yin not
            %updated
      
            proposed_log_posterior = log_posterior(proposed_log_likelihood,proposed_parameter,Delta_chain(:,:,k),delta_chain(k),a_tau,b_tau,V_not,v_not,D);
    
            if log(u) < proposed_log_posterior-previous_log_posterior
                theta_chain(k+1,block_index{b}) = proposed_parameter(block_index{b}); %save parameter to next iteration in chain & parameters stay in proposal
                previous_log_posterior = proposed_log_posterior;                
                acceptance(k,b) = 1;
                sigma_vpsc = vpsc_prop;
                       
            else
                theta_chain(k+1,block_index{b}) = theta_chain(k,block_index{b});
                acceptance(k,b) = 0;
                proposed_parameter(block_index{b}) = theta_chain(k,block_index{b}); %if proposed parameter denied, return to last accepted value
            end
        
        catch
            sprintf('Error in iteration (%i,%i)',k,b)
            kk = k+1;
            sprintf('\n Continuing to next iteration with parameter(%i,%i) = parameter(%i,%i)',kk,b,k,b)
            theta_chain(k+1,block_index{b}) = theta_chain(k,block_index{b}); %if an error: return to last accepted value
            proposed_parameter(block_index{b}) = theta_chain(k,block_index{b}); %return proposal to last accepted value
            acceptance(k,b) = 0;
            
            continue
            
        end
        
    end
    
    %sample delta and Lambda from full conditional (gibbs sampler)
        %samples directly from posterior, always accepted
        
        %sample xi (the sample precision)
       
        sq_diff_pars = zeros(D,D);
             for i = 1:S
                 sq_diff_pars = sq_diff_pars + (theta_chain(k+1,D*i-3:D*i) - theta_chain(k+1,S*D+1:end))'*(theta_chain(k+1,D*i-3:D*i) - theta_chain(k+1,S*D+1:end));
             end
          
%          C = sq_diff_pars/5; %sample covariance   
%          for i = 1:D
%             for j = 1:D
%                 Correlation(i,j) = C(i,j)/(sqrt(C(i,i))*sqrt(C(j,j)));
%             end
%          end
         

        % W = wishrnd(Sigma,df) generates a random matrix W having the Wishart distribution with covariance matrix Sigma and with df degrees of freedom
        Delta_chain(:,:,k+1) = wishrnd((V_not\eye(D) + sq_diff_pars)\eye(D),v_not + S);
                
        
%       sample tau (the error precision)
          
                         
            sq_diff = (sigma_obs-sigma_vpsc).^2;
            
           delta_chain(k+1) = gamrnd(S*N/2 + a_tau,...
                1/(b_tau + 1/2*(sum(sum(sq_diff))))); 
        
       %calculate log likelihood after gibbs step
       [proposed_log_likelihood,vpsc_prop] = log_likelihood(theta_chain(k+1,:),Delta_chain(:,:,k+1),delta_chain(k+1),sigma_obs,strain_inc,sigma_vpsc,b); 
            
       %calculate log posterior after gibbs step
       proposed_log_posterior = log_posterior(proposed_log_likelihood,theta_chain(k+1,:),Delta_chain(:,:,k+1),delta_chain(k+1),a_tau,b_tau,V_not,v_not,D);
                
       previous_log_posterior = proposed_log_posterior; 
        
    
        
    log_post(k+1) = previous_log_posterior; %only save log posterior at the end of the block cycle (once per iteration)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%The rest of the script is for generating diagnostic plots%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if mod(k,store) == 0
        
        %generate nuisance parameter distribution
        nuis_par_dist = nuis_par(strain_inc,sigma_obs,k,store,S,sigma_vpsc,nuis_par_dist,figure_path);
        
        
        % Genereate Posterior Predictive Distribution
         ppd_draw = post_pred_draw(N,S,D,theta_chain(k+1,:),Delta_chain(:,:,k+1),delta_chain(k+1),...
            block_index,strain_inc,sigma_obs,ppd_draw,k,store,figure_path);      
        
    end
    
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
        
        %fid = fopen('/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/iterations50_100.txt', 'a+');
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
        
        for zz = 1:num_blocks 
            [ml_voce] = diagnostic_test(theta_chain(1:k,block_index{zz}),log_post,k,zz,sigma_obs,strain_inc,0,iter_start,figure_path);
            ml_par{zz} = ml_voce;
        end
        
            
        disp('Diagnostic plots have been created in working directory')
        
        save(fullfile(figure_path,'100K.mat'))
        save(fullfile(figure_path,'current_cov_mat.mat'), 'cov_mat')
        
        %save('/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/iterations50_100.mat')
        %save('/users/PAS1064/osu8614/Bayesian Inference/Wishart_Gibbs_Step/figures/current_cov_mat.mat', 'cov_mat')
        current_parameters{1} = theta_chain(k,:);
        current_parameters{2} = Delta_chain(:,:,k);
        current_parameters{3} = delta_chain(k);
        current_parameters{4} = V_not;
        
        save( fullfile(figure_path,'current_parameters.mat'),'current_parameters')   
    
             
        catch 
        sprintf('Invalid file identifier in iteration (%i,%i), no diagnostics created. Continuing to next iteration',k,b)
        save(fullfile(figure_path,'100K.mat'))
        save(fullfile(figure_path,'current_cov_mat.mat'), 'cov_mat')
        current_parameters{1} = theta_chain(k,:);
        current_parameters{2} = Delta_chain(:,:,k);
        current_parameters{3} = delta_chain(k);
        current_parameters{4} = V_not;
        save( fullfile(figure_path,'current_parameters.mat'),'current_parameters')
 
        end     
        tic
    end
    close all
end