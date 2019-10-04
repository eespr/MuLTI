%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab script runs a transdimensional MCMC 
% inversion for S-Wave data 

% Author: Phil Livermore / Siobhan Killingbeck
% It is based on Matlab code written by Thomas Bodin.

% The physical model consists of internal layers, each with an associated shear wave
% speed vs defined by Voronoi nuclei.

% The domain is divided into a number of layers, num_layers (which could be one)
% each with its own prior distribution on Vs.

% Each layer has a special nuclei that cannot leave its layer ("confined" nuclei).
% npt is the number of "floating" nuclei that can change layer.

% The total number of nuclei is therefore npt_max + num_layers

clear all % clear all the variables previously allocated
close all % close all the figures

%%%%%%%%%%%%%%%%%%%%%%%%%LOAD DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Loading dispersion curves picked from 1D shot gathers%%%%%%%%%%%%%

load('4.real_data_picks.mat'); %load picked dispersion curves 
freq = data(:,1);% set frequency
data = data(:,2); % set observed data
nd =numel(data); % nd is the number of data points

%determine fitting error of picking dispersion curve in m/s
fitting_error = half_width(:,2); % set fitting error to half width of waveform's dispersion image

running_mode = 1; %1 -> find posterior; 0 -> Find priors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     VERY IMPORTANT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%th

burn_in=10000; % burn-in period
nsample=1000000; % total number of samples

%Uniform prior on depths for nuclei
priors.depth_min=0; % Cannot be changed, always from the surface.
priors.depth_max=40;
priors.npt_max = 30;
priors.npt_min = 0; 

num_layers = 3;

priors.layer_depths = zeros(num_layers-1,1); %the last layer depth is infinity (and is not defined).
priors.vsmin = zeros(num_layers,1);
priors.vsmax = zeros(num_layers,1);
priors.vp = zeros(num_layers,1);
priors.density = zeros(num_layers,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the priors and layer geometry
if num_layers == 3
    priors.layer_depths(1) = 2 ; % snow depth.
    priors.layer_depths(2) = 25.5 ; % ice depth.
    priors.vsmin(1) = 500; % Snow
    priors.vsmax(1) = 1700; 
    priors.vsmin(2) = 1700; % Ice
    priors.vsmax(2) = 1950;
    priors.vsmin(3) = 200; % rock
    priors.vsmax(3) = 2800;
    priors.vp(1) = 2500; % snow fixed vp
    priors.vp(2) = 3810; % ice fixed Vp
    priors.vp(3) = 4000; % rock fixed Vp
    priors.density(1) = 0.47; % snow fixed density
    priors.density(2) = 0.92; % ice fixed density
    priors.density(3) = 2.5; % rock fixed density
elseif num_layers == 2
    priors.layer_depths(1) = 2 ; % snow depth.
    priors.vsmin(1) = 500; % Snow
    priors.vsmax(1) = 1700;
    priors.vsmin(2) = 200; % rock
    priors.vsmax(2) = 2800;
    priors.vp(1) = 2500; % snow fixed vp    
    priors.vp(2) = 4000; % rock fixed Vp
    priors.density(1) = 0.47; % snow fixed density    
    priors.density(2) = 2.5; % rock fixed density
else % num_layers = 1 
    priors.vsmin(1) = 200; % wide range of constraints 
    priors.vsmax(1) = 2800;
    priors.vp(1) = 3500; % fixed vp   
    priors.density(1) = 2; % fixed density  
end

npt_init=1; % initial number of floating nuclei

sigma_change_vs=20; % std deviation of Gaussian proposal on Change vs value
sigma_move_depth = 1; % std deviation of Gaussian proposal on MOVE (change depth)
sigma_birth_vs = 400; % std deviation of Gaussian proposal on BIRTH
% (this number is also present in the DEATH
% acceptance term when taking in acount the reverse jump )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LESS IMPORTANT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default');
rng(4);
% You can change the 1 to any other number to change the seed.

% Define the limits of your model
x_min = priors.depth_min;
x_max = priors.depth_max;
dis=80; % steps to discretize the model. Trade-off between computational time and accuracy. 

show=10000; % show statistics of the chain every "show" samples
thin = 100; % thining

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
x =linspace(x_min,x_max,dis); % discretize the model
y = linspace(min(priors.vsmin), max(priors.vsmax), dis); %y limits are Vs priors limits
num=ceil((nsample-burn_in)*0.025/thin); % number of collected samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Preallocation of variables

b=0;
bb=0;
AV=zeros(dis,1);

AB=0;
AD=0;
PB=0;
PD=0;

AcV=0;
PV=0;
AP=0;
PP=0;

errors_gpdc = nan(nsample,(1 + (priors.npt_max+num_layers)*4));
best=zeros(dis,1);
val_min=zeros(dis,1);
val_max=zeros(dis,1);
ind_min_vs=zeros(dis,1);
ind_max_vs=zeros(dis,1);
hist_vs=zeros(dis,1);
hist_density=zeros(dis,1);
hierhist=zeros(nsample,1);
change_points=zeros(nsample*(priors.npt_max+num_layers),1);
cov=zeros(nsample,1);
nnuclei=zeros(nsample,1);
sup=zeros(dis,1);
inf=zeros(dis,1);
MINI_vs=zeros(dis,num);
MAXI_vs=zeros(dis,num);
CI_density=zeros((dis-1),(dis-1));

nnucleihist=zeros(priors.npt_max+num_layers,1);
nuclei_depths=zeros(priors.npt_max+num_layers,1);
nuclei_vs = zeros(priors.npt_max+num_layers,1);
nuclei_vp = zeros(priors.npt_max+num_layers,1);
nuclei_density = zeros(priors.npt_max+num_layers,1);
thickness = zeros(priors.npt_max+num_layers,1);
density = zeros(priors.npt_max+num_layers,1);
vs = zeros(priors.npt_max+num_layers,1);
vp = zeros(priors.npt_max+num_layers,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize - Define randomly the first model of the chain. Make sure
% that it satisfies the priors.

npt=npt_init;

    for i=1:npt+num_layers
% define the layer depth. The first num_layer nuclei are special: the ith
% nuclei must reside in the ith layer. 

        if i <= num_layers
            if i == 1
                top_of_layer = priors.depth_min;
                    if num_layers > 1
                        bottom_of_layer = priors.layer_depths(1);
                    else
                        bottom_of_layer = priors.depth_max;
                    end
                
            elseif i < num_layers
                top_of_layer = priors.layer_depths(i-1);
                bottom_of_layer = priors.layer_depths(i);
            else 
                top_of_layer = priors.layer_depths(i-1);
                bottom_of_layer = priors.depth_max;
            end
                
            nuclei_depths(i)= (bottom_of_layer + top_of_layer) / 2; % fix nuclei to be in middle of layer
        else    
            nuclei_depths(i)=priors.depth_min+rand*(priors.depth_max-priors.depth_min); % position of floating nuclei
        end
       
    % For each nuclei, find out which layer it is in:
        layer = num_layers;
        for j = 1:num_layers - 1
            if nuclei_depths(i) <= priors.layer_depths(j)
                layer = j;
                break
            end
        end              
      
     % the variable 'layer' is the layer of the nuclei:
       nuclei_vs(i)=priors.vsmin(layer)+rand*(priors.vsmax(layer)-priors.vsmin(layer)); % vs 
       nuclei_density(i)=priors.density(layer);
       nuclei_vp(i)=priors.vp(layer);
   
    end 

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE INITIAL MISFIT
% (Here, 'like' is the misfit )
like=0;
[thickness, density, vp, vs, priors_OK] = thicknesses_and_priors(nuclei_depths, nuclei_density, nuclei_vp, nuclei_vs, npt, num_layers, priors); 
    
%thickness(npt) = []; %remove last thickness for mat_disperse.  
if running_mode == 1
    forward_model = [freq zeros(length(data),5)];
    try
    forward_model = gpdc(thickness, vp.', vs.', (density*1000).', 'fV', freq.');
    forward_model = [forward_model(:, 1), rdivide(1, forward_model(:, 2:end))]; %convert forward model into velocity (m/s)
    catch %creating error file to save varibles which do not run through the gpdc code called errors_gpdc
        errors_gpdc(1,1:(1+length(thickness)+length(vp.')+length(vs.')+length((density*1000).'))) = [0, thickness, vp.', vs.', (density*1000).'];
    end %end
    
    %%%%%%%%%% computing multimodal misfit %%%%%%%%%%%
    min_misfit_all_modes = NaN(length(freq),1);

    for i = 1:length(freq) % frequency samples, this should match the frequency samples.
    misfit_fm = abs(data(i) - forward_model(i,2));
    misfit_m1 = abs(data(i) - forward_model(i,3));
    misfit_m2 = abs(data(i) - forward_model(i,4));
    %misfit_m3 = abs(data(i) - forward_model(i,5));
    %misfit_m4 = abs(data(i) - forward_model(i,6));
        if i == 17 %set first freq picked to the fundamental mode
            misfit_all_modes = [misfit_fm];
        else
            misfit_all_modes = [misfit_fm misfit_m1 misfit_m2];
        end
    min_misfit_all_modes(i,1) = min(misfit_all_modes);

    end %end multimodal misfit
    
    like = nansum( (min_misfit_all_modes).^2 ./(2 * fitting_error.^2) );
else
like = 1;
end

like_best=1e99;
like_init=like;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%% START RJ-MCMC SAMPLING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Total number of samples %i\n',nsample);
fprintf('Acceptance rates:\n');
fprintf('Iteration  Change Vs    Move Depth    Birth     Death\n');
for s=1:nsample
    out=1;
    % Print statistics of the chain, The best is to "tune" these
    % ratios to 44 %. (Rosental 2000).
    if (mod(s,show)==0)
        number_of_samples = s;
        number_of_nuclei = npt;
        if (s>burn_in)
            fprintf('%7i     %5.2f         %5.2f       %5.2f     %5.2f\n',s, 100*AcV/PV, 100*AP/PP, 100*AB/PB,100*AD/PD);
        end
    end
     
    
    birth=0;
    move=0;
    death=0;
    
    nuclei_vs_prop=nuclei_vs;
    nuclei_vp_prop = nuclei_vp;
    nuclei_density_prop = nuclei_density;
    nuclei_depths_prop = nuclei_depths;
    
    like_prop = like;
    %----------------------------------------------------------------------
    % Every even iteration, propose a new changed parameter value
    if (mod(s,2)==0) % Change Value
        if (s>burn_in)
            PV=PV+1;
        end
        npt_prop = npt;
        ind=ceil(rand*(npt+num_layers));
        nuclei_vs_prop(ind) = nuclei_vs(ind) + randn * sigma_change_vs;

        %-----------------------------------------------------------------------
        % Every odd iteration change the nuclei tesselation
    else % Change position

        %u=1; % turning off birth/death
        u=rand; % Chose randomly between 3 different types of moves
        if (u<0.333) % BIRTH ++++++++++++++++++++++++++++++++++++++
            birth=1;
            if (s>burn_in)
                PB=PB+1;
            end
            npt_prop = npt+1;
            nuclei_depths_prop(1:npt+num_layers) = nuclei_depths(1:npt+num_layers);
            nuclei_depths_prop(npt+num_layers+1) = priors.depth_min+rand*(priors.depth_max-priors.depth_min);
            ind=whichnuclei(nuclei_depths(1:npt+num_layers),nuclei_depths_prop(npt+num_layers+1), num_layers, priors);

      
            nuclei_density_prop(npt+num_layers+1)=nuclei_density(ind);
            nuclei_vp_prop(npt+num_layers+1)=nuclei_vp(ind);
            nuclei_vs_prop(npt+num_layers+1)=nuclei_vs(ind)+randn*sigma_birth_vs;
            
 % find which layer it's in and find the product of Priors:
 
        layer = num_layers;
        Prod_delta_prior = priors.vsmax(layer)-priors.vsmin(layer);
        
        for j = 1:num_layers - 1
            if nuclei_depths_prop(npt+num_layers+1) <= priors.layer_depths(j)
                layer = j;
                Prod_delta_prior = (priors.vsmax(j)-priors.vsmin(j) );
                break
            end
        end  
               
       prob = 1.0 / (sigma_birth_vs*sqrt(2*pi)) * exp(-( nuclei_vs_prop(num_layers+npt+1) - nuclei_vs(ind) )^2/(2*sigma_birth_vs^2));
            
        elseif (u<0.666) % DEATH +++++++++++++++++++++++++++++++++++++++++
            death=1;
            if (s>burn_in)
                PD=PD+1;
            end
            
            npt_prop = npt-1;   
  % choose a floating nuclei to remove
            ind=ceil(rand*npt)+num_layers;
            
            nuclei_depths_prop(1:num_layers+npt-1) = [nuclei_depths(1:ind-1) ; nuclei_depths(ind+1:num_layers+npt)];
            nuclei_density_prop(1:num_layers + npt-1)= [nuclei_density(1:ind-1) ; nuclei_density(ind+1:num_layers+npt)];
            nuclei_vs_prop(1:num_layers + npt-1)= [nuclei_vs(1:ind-1) ; nuclei_vs(ind+1:num_layers+npt)];
            nuclei_vp_prop(1:num_layers + npt-1)= [nuclei_vp(1:ind-1) ; nuclei_vp(ind+1:num_layers+npt)];
           
                death_pt_density = nuclei_density(ind);
                death_pt_vs = nuclei_vs(ind);
                death_pt_vp = nuclei_vp(ind);
                death_pt_depth = nuclei_depths(ind);
                
                        
    %GET prob
                node=whichnuclei(nuclei_depths_prop(1:npt_prop+num_layers),death_pt_depth, num_layers,priors);
                %prob=(1/(sigmav*sqrt(2*pi)))*exp(-(pt(ind,2)-pt_prop(node,2))^2/(2*sigmav^2));
              
 % find which layer it's in and find the product of Priors:
        layer = num_layers;
        Prod_delta_prior = priors.vsmax(layer)-priors.vsmin(layer);
        for j = 1:num_layers - 1
            if death_pt_depth <= priors.layer_depths(j)
                layer = j;
                Prod_delta_prior = priors.vsmax(j)-priors.vsmin(j) ;
                break
            end
        end  
                
        % the case of npt_prop = -1 is a rather special case but it results in an error. 
        % when num_layers = 1. It is never excepted.
        if npt_prop == -1 
            prob = 1; %set to anything.
        else
            prob = 1.0 / (sigma_birth_vs*sqrt(2*pi)) * exp(-( nuclei_vs_prop(node) - death_pt_vs )^2/(2*sigma_birth_vs^2));
        end
            
            
            
        else % MOVE +++++++++++++++++++++++++++++++++++++++++++++++++++++++
            if (s>burn_in)
                PP=PP+1;
            end
            move=1;
            npt_prop = npt;
 % choose the nuclei to move
            ind=ceil(rand*(npt+num_layers));
 
            if num_layers == 1 || ind > num_layers %if ind is a 'floating' nuclei or depth constraints are not applied, move nuclei randomly using sigma_move_depth
                nuclei_depths_prop(ind) = nuclei_depths(ind)+randn*sigma_move_depth;
            else %if ind is a 'confined' nuclei, move nuclei randomly within the range of the layer depths
                if ind == 1
                top_of_layer = priors.depth_min;
                bottom_of_layer = priors.layer_depths(1);
                elseif ind < num_layers
                top_of_layer = priors.layer_depths(ind-1);
                bottom_of_layer = priors.layer_depths(ind);
                else
                top_of_layer = priors.layer_depths(ind-1);
                bottom_of_layer = priors.depth_max;
                end
                nuclei_depths_prop(ind) = (bottom_of_layer-top_of_layer).*rand(1) + top_of_layer;
            end
                           
% Find move probability

% find which layer the nuclei is currently in:
        layer = num_layers;
        move_prob1 = priors.vsmax(layer)-priors.vsmin(layer);
        for j = 1:num_layers - 1
            if nuclei_depths(ind) <= priors.layer_depths(j)
                layer = j;
                move_prob1 = priors.vsmax(j)-priors.vsmin(j) ;
                break
            end
        end  
        
% find which layer the nuclei will move to:
        layer = num_layers;
        move_prob2 = priors.vsmax(layer)-priors.vsmin(layer);
        for j = 1:num_layers - 1
            if nuclei_depths_prop(ind) <= priors.layer_depths(j)
                layer = j;
                move_prob2 = priors.vsmax(j)-priors.vsmin(j) ;
                break
            end
        end  
        move_prob = move_prob1 / move_prob2;
        
        end   
          
    end % Change the position
    %----------------------------------------------------------------------
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE MISFIT OF THE PROPOSED MODEL
    % If the proposed model is not outside the bounds of the uniform prior,
    % compute its misfit : "like_prop"
    if out==1
        like_prop=0;
        [thickness, density, vp, vs, priors_OK] = thicknesses_and_priors(nuclei_depths_prop, nuclei_density_prop, nuclei_vp_prop, nuclei_vs_prop, npt_prop, num_layers, priors);

        if priors_OK == 0 
            out = 0;
            like_prop = 0;
        else
        
            if running_mode == 1
                forward_model = [freq zeros(length(data),5)];
                    try
                    forward_model = gpdc(thickness, vp.', vs.', (density*1000).', 'fV', freq.');
                    forward_model = [forward_model(:, 1), rdivide(1, forward_model(:, 2:end))]; %convert forward model into velocity (m/s)
                    catch
                       errors_gpdc(s+1,1:(1+length(thickness)+length(vp.')+length(vs.')+length((density*1000).'))) = [s, thickness, vp.', vs.', (density*1000).'];
                    end %end
                    
                    %%%%%%%%%%% computing multimodal misfit %%%%%%%%%%%
                    min_misfit_all_modes = NaN(length(freq),1);

                    for i = 1:length(freq) % frequency samples, this should match the frequency samples.

                    misfit_fm = abs(data(i) - forward_model(i,2));
                    misfit_m1 = abs(data(i) - forward_model(i,3));
                    misfit_m2 = abs(data(i) - forward_model(i,4));
                    %misfit_m3 = abs(data(i) - forward_model(i,5));
                    %misfit_m4 = abs(data(i) - forward_model(i,6));
                        if i == 17 %set first freq picked to the fundamental mode
                            misfit_all_modes = [misfit_fm];
                        else
                            misfit_all_modes = [misfit_fm misfit_m1 misfit_m2];
                        end
                    min_misfit_all_modes(i,1) = min(misfit_all_modes);
                    end
                    %end
                    like_prop = nansum( (min_misfit_all_modes).^2 ./(2 * fitting_error.^2) );
            else
            like_prop = 1;
            end
        end
        
    end %if (out==1)
    
       
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SEE WHETHER MODEL IS ACCEPTED
    
    accept=0;
 % This avoids runtime errors, as if out=0 Matlab insists on evaluating
 % Prod_delta_prior even though it never affects the calculation if out==0
 % (for the result is always 0).
 
    if out == 0 
    Prod_delta_prior = 1;
    prob = 1;
    end
    
    % THe acceptance term takes different
    % values according the the proposal that has been made.
    
    if (birth==1)        
        if (rand<((1/(Prod_delta_prior*prob))*exp(log(out)-like_prop+like)))
            accept=1;
            if (s>burn_in)
                AB=AB+1;
            end
        end
    elseif (death==1)
        
        if (rand<(Prod_delta_prior*prob*exp(log(out)-like_prop+like)))
            accept=1;
            if (s>burn_in)
                AD=AD+1;
            end
        end
        
    elseif (move == 1) % NO JUMP, i.e no change in dimension
        
        if (rand<(move_prob * exp(log(out)-like_prop+like)))
            accept=1;
            if (s>burn_in)
            AP=AP+1;       
            end %if (s>burn_in)
        end
        
    else %change v_s
       if  (rand<exp(log(out)-like_prop+like))
        accept=1;
            if (s>burn_in)
            AcV=AcV+1;       
            end %if (s>burn_in)      
    end
    end
    
    % If accept, update the values
    if (accept==1)
        npt=npt_prop;
        nuclei_depths = nuclei_depths_prop;
        nuclei_density = nuclei_density_prop;
        nuclei_vs = nuclei_vs_prop;
        nuclei_vp = nuclei_vp_prop;
        like=like_prop;
        for i=1:dis
            ind=whichnuclei(nuclei_depths(1:npt+num_layers),x(i),num_layers,priors);
            hist_vs(i)=nuclei_vs(ind);
        end
        [N]=histcounts2(x,hist_vs',x,y);
        vs_edge=y(1:(dis-1));
        depth_edge=x(1:(dis-1));
    end
    
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We collect the samples for the ensemble solution
    
    if (s>burn_in)
        if (mod(s,thin)==0)
            b=b+1;
            % DO THE AVERAGE
            
            for i=1:dis
                ind=whichnuclei(nuclei_depths,x(i), num_layers,priors);
                
                AV(i,1)=AV(i,1)+nuclei_vs(ind);               
                
                % Do the 95% credible interval for vs
                if (b<=num)
                    MINI_vs(i,b)=nuclei_vs(ind);
                    MAXI_vs(i,b)=nuclei_vs(ind);
                    if (b==num)
                        [val_min_vs(i) ind_min_vs(i)]=min(MAXI_vs(i,:));
                        [val_max_vs(i) ind_max_vs(i)]=max(MINI_vs(i,:));
                    end
                    
                else
                    if (nuclei_vs(ind)>val_min_vs(i))
                        MAXI_vs(i,ind_min_vs(i))=nuclei_vs(ind);
                        [val_min_vs(i) ind_min_vs(i)]=min(MAXI_vs(i,:));
                    end
                    if (nuclei_vs(ind)<val_max_vs(i))
                        MINI_vs(i,ind_max_vs(i))=nuclei_vs(ind);
                        [val_max_vs(i) ind_max_vs(i)]=max(MINI_vs(i,:));
                    end
                end
                
            end
            nnucleihist(npt+num_layers,1)=nnucleihist(npt+num_layers,1)+1;
            %Do the histogram on change points
            nuclei=nuclei_depths(1:npt+num_layers);
            nuclei=sort(nuclei);
            for i = 1:npt-1+num_layers
                bb=bb+1;
                cp= (nuclei(i+1)+nuclei(i))/2;
                change_points(bb)=cp;
            end
        end
        
        if accept==1 %do histergram on CI density Vs accepted
            CI_density = CI_density + N;
        end
        
    end %if burn-in
    
    cov(s)=like; % Convergence of the misfit
    nnuclei(s)=npt+num_layers; % Convergence of number of nuclei
            
    % Get the best model
    if priors_OK ==1 && (like<like_best) 
        depths_best = nuclei_depths;
        density_best = nuclei_density;
        vs_best = nuclei_vs;
        vp_best = nuclei_vp;
        npt_best = npt;
        like_best = like;
        if running_mode ==1
        forward_model_best = forward_model;
        end
    end   
   
    
end% the Sampling of the mcmc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%% Take average and credible intervals %%%%%%%%%%%
AV=AV./b;
%
best_model_vs = zeros(1,dis);

for i=1:dis
    ind=whichnuclei(depths_best(1:num_layers+npt_best),x(i),num_layers, priors);
    best_model_vs(i)=vs_best(ind);
end

%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the best, average and mode solutions

% Vs
for i=1:dis
    [val_min ind_min]=min(MAXI_vs(i,:));
    [val_max ind_max]=max(MINI_vs(i,:));
    sup(i)=val_min;
    inf(i)=val_max;
end
%
x2 = [x'; flipud(x')];
inbetween = [inf; flipud(sup)];
%%%%%%%%%%%%%%%%%%%%%%
N = histcounts2(x2,inbetween,x,y);
for j=1:(dis-1)
    f=find(N(j,:),1,'first');
    l=find(N(j,:),1,'last');
    for i = f:l
        N(j,i)=1;
    end
end
CI_density_limit=CI_density.*N;
%%%%%%%%%%%%%%%%%%%%%%

figure
imagesc(vs_edge,depth_edge,CI_density_limit)
hold on
plot(inbetween, x2,'g');
hold on
plot(best_model_vs,x,'k','LineWidth',2);
hold on
plot(AV(:,1),x,'r','LineWidth',2);
title('Vs density plot');
h=legend('95% CI','best model sampled', 'true model', 'average solution');
%set(gca,'Ydir','normal')
set(h,'Interpreter','none','Location','NorthOutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalising density plot
sizeCI = size(CI_density_limit);
len = sizeCI(1,1);
wid = sizeCI(1,2);
CI_density_limitN = zeros(len,wid);
for i = 1:len
    m = sum(CI_density_limit(i,:));
    for j = 1:wid
        CI_density_limitN(i,j) = CI_density_limit(i,j)/m;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find mode solution 
%calculating modal solution from CI_density plot
    mode = zeros((dis-1),(dis-1));
    for i=1:(dis-1)
        [M,I]=max(CI_density_limit(i,:));
        mode(i,I)=1;
    end
    vs_mode=zeros((dis-1),(dis-1));
    for i=1:(dis-1)
    vs_mode(i,:) = mode(i,:).* vs_edge;
    end   
%removing the zeros in the solution
    vs_mode_n0 = nan((dis-1),1);
    for i=1:(dis-1)
    vs_mode_n0(i,1) = max(vs_mode(i,:));
    end
    vs_mode_n0(dis,1) = vs_mode_n0(dis-1,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VS plot (Normalised) 
figure
imagesc(vs_edge,depth_edge,CI_density_limitN)
hold on
plot(best_model_vs,x,'k','LineWidth',2);
%plot(inbetween, x2,'g');
%hold on
%plot(vs_mode_n0,x,'k','LineWidth',2);
%hold on
%plot(AV(:,1),x,'r','LineWidth',2);
title('Shear Wave Velocity Plot');
% Create ylabel
ylabel('Depth (m)');
% Create xlabel
xlabel('Vs (m/s)');
colormap jet
colorbar
caxis([0 0.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot best model modal dispersion curves and observed data 
if running_mode == 1
figure;
plot( freq, data, '*');
hold on;
plot( freq, forward_model_best(:,2:6), 'r');
% Create ylabel
ylabel('Phase velocity (m/s)');
% Create xlabel
xlabel('Frequency (Hz)');
end


%%%%%%%%%%%%%%% Plot Statistics of the chain %%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
semilogy(cov);
hold on
%line([burn_in burn_in],[0 cov(1)],'LineWidth',4,'Color',[1 0 0]);
xlabel('iterations','FontSize',14)


title('Data Misfit','FontSize',16)
subplot(2,1,2);
line([burn_in burn_in],[0 priors.npt_max+num_layers],'LineWidth',4,'Color',[1 0 0]);
hold on
plot(nnuclei)
xlabel('iterations','FontSize',14)
title('Number of nuclei','FontSize',16)


%%%%%%%%%%%%%% Plot Marginal posteriors %%%%%%%%%%%%%%%%%%%

% Plot histogram on change points
figure
subplot(2,1,1)
hist(change_points(1:bb),500)
title('Probability of change points','FontSize',14)

% Plot histogram on number of nuclei
subplot(2,1,2);
bar(nnucleihist)
title('Posterior Distribution on number of nuclei ','FontSize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%SAVING WORKSPACE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input file name remember to change name for each inversion
save('4.Constrained_Real_data.mat') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      THE END                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

