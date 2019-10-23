 function [vpsc_out] = VPSC(tau0,tau1,theta0,theta1,strain_inc)

   %% Generate new crystal file
    rate = 10;
    voce = {tau0,tau1,theta0,theta1};
       
    str = fileread('template_FCC_cu.txt');
    fid = fopen('old_FCC_cu.sx','wt');
    fprintf(fid,str,rate,voce{:});
    fclose(fid);

    %remove extra spaces in FCC_cu file
    filecontent = fileread('old_FCC_cu.sx');
    newcontent = regexprep(filecontent, '\r\n', '\n');
    fid = fopen('FCC_cu.sx', 'w');
    fwrite(fid, newcontent);
    fclose(fid);
 
    %% Run VPSC code
    % ATTENTION!!!!!!!
    % the '&' is needed at the end of the call to the .exe file in LINUX
    
     evalc('!vpsc7');        
%      !vpsc7 

     
           
     ES = importdata('STR_STR.OUT');
     STR_STR = ES.data;
            
     %define sample points for interpolation
     x = STR_STR(:,1);
     v = STR_STR(:,2);
     
%     figure
%     plot(x,v)
     
%     %query points
     xq = strain_inc;
     
%     %interpolate VPSC dtat at query points and plot
     vq = interp1(x,v,xq);
     vpsc_out = vq; %simulated data: VPSC results for query points
          
end