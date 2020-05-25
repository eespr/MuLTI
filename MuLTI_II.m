%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  MuLTI                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab script runs a transdimensional MCMC 
% inversion for S-Wave data 

% Author: Phil Livermore / Siobhan Killingbeck
% School of Earth and Environemnt, The University of Leeds
% It is based on Matlab code written by Thomas Bodin.

% The physical model consists of internal layers, each with an associated Vs and 
% Vp defined by Voronoi nuclei.

% The domain is divided into a number of layers, num_layers (which could be one)
% each with its own prior distribution on Vs and Vp.

% Each layer has a special nuclei that cannot leave its layer ("confined" nuclei).
% npt is the number of "floating" nuclei that can change layer.

% The total number of nuclei is therefore npt_max + num_layers

clear all % clear all the variables previously allocated
close all % close all the figures

%%%%%%%%%%%%%%%%%%%%%%%%%LOAD DATA%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%Loading dispersion curves picked from 1D shot gathers%%%%%%%%%%%%%

load('input_data.mat'); %load picked dispersion curves
load('vpdata.mat'); %load vp priors as variable named "vpdata" and 
% column 1:vpdepths, column 2:vpmean, column 3:vpsd
freq = data(:,1);% set frequency
data = data(:,2); % set observed data
nd =numel(data); % nd is the number of data points

%determine fitting error of picking dispersion curve in m/s
fitting_error = linspace(10,40,nd).'; % set fitting error to half width of waveform's dispersion image

running_mode = 1; %1 -> find posterior; 0 -> Find priors.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     VERY IMPORTANT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%th

burn_in=10000; % burn-in period
nsample=500000; % total number of samples

%Uniform prior on depths for nuclei
priors.depth_min=0; % Cannot be changed, always from the surface.
priors.depth_max=50;
priors.npt_max = 50;
priors.npt_min = 0;
dis=100; % steps to discretize the model. Trade-off between computational time and accuracy.
num_layers = 3; % depth constraining parameter, if set to 1 no depth constraints will be applied.

% Define the limits of your model
x_min = priors.depth_min;
x_max = priors.depth_max;
x =linspace(x_min,x_max,dis); % discretize the model
priors.layer_depths = zeros(num_layers-1,1); %the last layer depth is infinity (and is not defined).
priors.vsmin = zeros(num_layers,1);
priors.vsmax = zeros(num_layers,1);
priors.vp = zeros(num_layers,1);
priors.density = zeros(num_layers,1);
priors.vpmin = min(vpdata(:,2));%define Vp grid size for pdf
priors.vpmax = max(vpdata(:,2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the priors and layer geometry
if num_layers == 3
    priors.layer_depths(1) = 18 ; % layer 1 depth.    
    priors.layer_depths(2) = 27 ; % layer 2 depth.
    
    priors.vsmin(1) = 500; % layer 1
    priors.vsmax(1) = 1800; 
    priors.vsmin(2) = 800; % layer 2
    priors.vsmax(2) = 1800;
    priors.vsmin(3) = 1500; % layer 3
    priors.vsmax(3) = 2100;   
    
    priors.density(1) = 0.716; % layer 1 fixed density
    priors.density(2) = 0.882; % layer 2 fixed density
    priors.density(3) = 0.893; % layer 3 fixed density
elseif num_layers == 2
    priors.layer_depths(1) = 18 ; % layer 1 depth.    
    
    priors.vsmin(1) = 500; % layer 1
    priors.vsmax(1) = 1800; 
    priors.vsmin(2) = 800; % layer 2
    priors.vsmax(2) = 2100;
    
    priors.density(1) = 0.716; % layer 1 fixed density
    priors.density(2) = 0.885; % layer 2 fixed density
else 
    priors.vsmin(1) = 500; % layer 1
    priors.vsmax(1) = 2100;  
    
    priors.density(1) = 0.885; % layer 1 fixed density
end
    priors.PRmin(1:num_layers) = 0.1;% PR limits
    priors.PRmax(1:num_layers) = 0.49;
    
    priors.vpdepth = x.';% Vp depths 
    priors.vpmean = interp1(vpdata(:,1), vpdata(:,2), priors.vpdepth);% Vp mean values at each depth
    priors.vpsd = interp1(vpdata(:,1), vpdata(:,3), priors.vpdepth);% Vp standard deviation 

npt_init=1; % initial number of floating nuclei

sigma_change_vs=200; % std deviation of Gaussian proposal on Change vs value
sigma_move_depth = 4; % std deviation of Gaussian proposal on MOVE (change depth)
sigma_birth_vs = 200; % std deviation of Gaussian proposal on BIRTH
% (this number is also present in the DEATH acceptance term when taking in acount the reverse jump )
sigma_change_vp = 400;
sigma_birth_vp= 200;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LESS IMPORTANT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rng('default');
rng(1);
% You can change the 1 to any other number to change the seed.

show=1000; % show statistics of the chain every "show" samples
thin = 1000; % thining

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
y = linspace(min(priors.vsmin), max(priors.vsmax), dis); %y limits are Vs priors limits
yp = linspace(priors.vpmin, priors.vpmax, dis); %yp limits are Vp grid size defined by user
yVpVs = linspace(0, 5, dis);
yPR = linspace(0, 0.5, dis);
num=ceil((nsample-burn_in)*0.025/thin); % number of collected samples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      Preallocation of variables

b=0;
bb=0;
AV=zeros(dis,1);
AV_vp=zeros(dis,1);

AB=0;
AD=0;
PB=0;
PD=0;

AcV=0;
AcVp=0;
PV=0;
PVp=0;
AP=0;
PP=0;

errors_gpdc = nan(nsample,(1 + (priors.npt_max+num_layers)*4));
best=zeros(dis,1);
val_min=zeros(dis,1);
val_max=zeros(dis,1);
ind_min_vs=zeros(dis,1);
ind_max_vs=zeros(dis,1);
ind_min_vp=zeros(dis,1);
ind_max_vp=zeros(dis,1);
hist_vs=zeros(dis,1);
hist_vp=zeros(dis,1);
hist_density=zeros(dis,1);
hist_VpVs=zeros(dis,1);
hist_PR=zeros(dis,1);
hierhist=zeros(nsample,1);
change_points=zeros(nsample*(priors.npt_max+num_layers),1);
cov=zeros(nsample,1);
nnuclei=zeros(nsample,1);
sup=zeros(dis,1);
inf=zeros(dis,1);
sup_vp=zeros(dis,1);
inf_vp=zeros(dis,1);
MINI_vs=zeros(dis,num);
MAXI_vs=zeros(dis,num);
MINI_vp=zeros(dis,num);
MAXI_vp=zeros(dis,num);
CI_density=zeros((dis-1),(dis-1));
CI_density_vp=zeros((dis-1),(dis-1));
CI_density_VpVs=zeros((dis-1),(dis-1));
CI_density_PR=zeros((dis-1),(dis-1));

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
       nuclei_vp(i) = interp1(vpdata(:,1), vpdata(:,2), nuclei_depths(i));
   
    end 

    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE INITIAL MISFIT
% (Here, 'like' is the misfit )
like=0;
[thickness, density, vp, vs, priors_OK] = thicknesses_and_priors_III(nuclei_depths, nuclei_density, nuclei_vp, nuclei_vs, npt, num_layers, priors); 
     
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
        if i == 1 %set first freq picked to the fundamental mode
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
fprintf('Iteration  Change Vs   Change Vp    Move Depth    Birth     Death\n');
for s=1:nsample
    out=1;
    % Print statistics of the chain, The best is to "tune" these
    % ratios to 44 %. (Rosental 2000).
    if (mod(s,show)==0)
        number_of_samples = s;
        number_of_nuclei = npt;
        if (s>burn_in)
            fprintf('%7i     %5.2f      %5.2f       %5.2f       %5.2f     %5.2f\n',s, 100*AcV/PV, 100*AcVp/PVp, 100*AP/PP, 100*AB/PB,100*AD/PD);
        end
    end
     
    
    birth=0;
    move=0;
    death=0;
    change_Vs = 0;
    change_Vp = 0;
    
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
        u = rand; % choose randomly between changing Vp (prob 1/6) and Vs (prob 1/3)
        if (u<0.333) %change Vp
            change_Vp = 1;
            if (s>burn_in)
            PVp=PVp+1;
            end
            nuclei_vp_prop(ind) = nuclei_vp(ind) + randn * sigma_change_vp;
            
            mean_vp = interp1(vpdata(:,1), vpdata(:,2), nuclei_depths_prop(ind)); 
            sigma_vp = interp1(vpdata(:,1), vpdata(:,3), nuclei_depths_prop(ind));% sample vp at depth nuclei_depths_prop(ind) 
            Prod_delta_prior_prop = 1.0 / (sigma_vp*sqrt(2*pi)) * exp(-( nuclei_vp_prop(ind) - mean_vp )^2/(2*sigma_vp^2));
            Prod_delta_prior_current = 1.0 / (sigma_vp*sqrt(2*pi)) * exp(-( nuclei_vp(ind) - mean_vp )^2/(2*sigma_vp^2));
            Prod_delta_prior = Prod_delta_prior_prop / Prod_delta_prior_current;
        else %change Vs
             change_Vs = 1;
             nuclei_vs_prop(ind) = nuclei_vs(ind) + randn * sigma_change_vs;
             Prod_delta_prior = 1.0;
        end

        %-----------------------------------------------------------------------
        % Every odd iteration change the nuclei tesselation
    else % Change position

        u=rand; % Chose randomly between 3 different types of moves
        %u=1;%fix number layers to initial 
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
            nuclei_vp_prop(npt+num_layers+1)=nuclei_vp(ind)+randn*sigma_birth_vp;
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
        mean_vp = interp1(vpdata(:,1), vpdata(:,2), nuclei_depths_prop(npt+num_layers+1)); 
        sigma_vp = interp1(vpdata(:,1), vpdata(:,3), nuclei_depths_prop(npt+num_layers+1));% sample vp at depth nuclei_depths_prop(npt+num_layers+1) 
        Prod_delta_prior_vp = 1.0 / (sigma_vp*sqrt(2*pi)) * exp(-( nuclei_vp_prop(num_layers+npt+1) - mean_vp )^2/(2*sigma_vp^2));
        Prod_delta_prior = Prod_delta_prior / Prod_delta_prior_vp;
        prob_vs = 1.0 / (sigma_birth_vs*sqrt(2*pi)) * exp(-( nuclei_vs_prop(num_layers+npt+1) - nuclei_vs(ind) )^2/(2*sigma_birth_vs^2));
        prob_vp = 1.0 / (sigma_birth_vp*sqrt(2*pi)) * exp(-( nuclei_vp_prop(num_layers+npt+1) - nuclei_vp(ind) )^2/(2*sigma_birth_vp^2));
        prob = prob_vs*prob_vp;
            
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
        
        mean_vp = interp1(vpdata(:,1), vpdata(:,2), death_pt_depth); 
        sigma_vp = interp1(vpdata(:,1), vpdata(:,3), death_pt_depth);% sample vp at depth death_pt_depth
        Prod_delta_prior_vp = 1.0 / (sigma_vp*sqrt(2*pi)) * exp(-( death_pt_vp - mean_vp )^2/(2*sigma_vp^2));
        Prod_delta_prior = Prod_delta_prior / Prod_delta_prior_vp;
        
        % the case of npt_prop = -1 is a rather special case but it results in an error. 
        % when num_layers = 1. It is never excepted.
        if npt_prop == -1 
            prob = 1; %set to anything.
        else
            prob_vs = 1.0 / (sigma_birth_vs*sqrt(2*pi)) * exp(-( nuclei_vs_prop(node) - death_pt_vs )^2/(2*sigma_birth_vs^2));
            prob_vp= 1.0 / (sigma_birth_vp*sqrt(2*pi)) * exp(-( nuclei_vp_prop(node) - death_pt_vp )^2/(2*sigma_birth_vp^2));
            prob = prob_vs * prob_vp;
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
        
        mean_vp = interp1(vpdata(:,1), vpdata(:,2), nuclei_depths(ind));
        sigma_vp = interp1(vpdata(:,1), vpdata(:,3), nuclei_depths(ind));% sample vp at depth nuclei_depths(ind)
        prob_vp1 = 1.0 / (sigma_vp*sqrt(2*pi)) * exp(-( nuclei_vp(ind) - mean_vp )^2/(2* sigma_vp ^2));
        
% find which layer the nuclei will move to:
        layer = num_layers;
        move_prob2 = priors.vsmax(layer)-priors.vsmin(layer);
        
        for j = 1:num_layers - 1
            if nuclei_depths_prop(ind) <= priors.layer_depths(j)
                layer = j;
                move_prob2 = priors.vsmax(j)-priors.vsmin(j);
                break
            end
        end
        
        mean_vp = interp1(vpdata(:,1), vpdata(:,2), nuclei_depths_prop(ind));
        sigma_vp = interp1(vpdata(:,1), vpdata(:,3), nuclei_depths_prop(ind));% sample vp at depth nuclei_depths_prop(ind)
        prob_vp2 = 1.0 / (sigma_vp*sqrt(2*pi)) * exp(-( nuclei_vp_prop(ind) - mean_vp )^2/(2* sigma_vp ^2));

        move_prob = move_prob1/move_prob2*prob_vp2/prob_vp1;
        
        end   
          
    end % Change the position
    %----------------------------------------------------------------------
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % COMPUTE MISFIT OF THE PROPOSED MODEL
    % If the proposed model is not outside the bounds of the uniform prior,
    % compute its misfit : "like_prop"
    if out==1
        like_prop=0;
        [thickness, density, vp, vs, priors_OK] = thicknesses_and_priors_III(nuclei_depths_prop, nuclei_density_prop, nuclei_vp_prop, nuclei_vs_prop, npt_prop, num_layers, priors);

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
                        if i == 1 %set first freq picked to the fundamental mode
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
       if  (rand<Prod_delta_prior * (exp(log(out)-like_prop+like)))
        accept=1;
        if (s>burn_in)
                if(change_Vs == 1)
                    AcV=AcV+1;   
                else %( change_Vp == 1)
                    AcVp = AcVp + 1;
                end   
        end
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
    end
    
    %calculate elastic moduli
    VpVs = nuclei_vp./nuclei_vs;
    PR = 0.5*(((VpVs.^2)-2)./((VpVs.^2)-1));
    
        for i=1:dis
            ind=whichnuclei(nuclei_depths(1:npt+num_layers),x(i),num_layers,priors);
            hist_vs(i)=nuclei_vs(ind);
            hist_vp(i)=nuclei_vp(ind);
            hist_VpVs(i)=VpVs(ind);
            hist_PR(i)=PR(ind);
        end
        
        [N]=histcounts2(x,hist_vs',x,y);
        [M]=histcounts2(x,hist_vp',x,yp);
        [V]=histcounts2(x,hist_VpVs',x,yVpVs);
        [P]=histcounts2(x,hist_PR',x,yPR);
        
        vs_edge=y(1:(dis-1));
        vp_edge=yp(1:(dis-1));
        VpVs_edge=yVpVs(1:(dis-1));
        PR_edge=yPR(1:(dis-1));
        depth_edge=x(1:(dis-1));
           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % We collect the samples for the ensemble solution
    
    if (s>burn_in)
        if (mod(s,thin)==0)
            b=b+1;
            % DO THE AVERAGE
            
            for i=1:dis
                ind=whichnuclei(nuclei_depths,x(i), num_layers,priors);
                
                AV(i,1)=AV(i,1)+nuclei_vs(ind);
                AV_vp(i,1)=AV_vp(i,1)+nuclei_vp(ind); 
                
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
                
                % Do the 95% credible interval for vp
                if (b<=num)
                    MINI_vp(i,b)=nuclei_vp(ind);
                    MAXI_vp(i,b)=nuclei_vp(ind);
                    if (b==num)
                        [val_min_vp(i) ind_min_vp(i)]=min(MAXI_vp(i,:));
                        [val_max_vp(i) ind_max_vp(i)]=max(MINI_vp(i,:));
                    end
                    
                else
                    if (nuclei_vp(ind)>val_min_vp(i))
                        MAXI_vp(i,ind_min_vp(i))=nuclei_vp(ind);
                        [val_min_vp(i) ind_min_vp(i)]=min(MAXI_vp(i,:));
                    end
                    if (nuclei_vp(ind)<val_max_vp(i))
                        MINI_vp(i,ind_max_vp(i))=nuclei_vp(ind);
                        [val_max_vp(i) ind_max_vp(i)]=max(MINI_vp(i,:));
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
        
            CI_density = CI_density + N;
            CI_density_vp = CI_density_vp + M;
            CI_density_VpVs = CI_density_VpVs + V;
            CI_density_PR = CI_density_PR + P;
        
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

%%%%%%%%%%% Take average and credible intervals %%%%%%%%%%%
AV=AV./b;
AV_vp=AV_vp./b;
%
best_model_vs = zeros(1,dis);

for i=1:dis
    ind=whichnuclei(depths_best(1:num_layers+npt_best),x(i),num_layers, priors);
    best_model_vs(i)=vs_best(ind);
end

%%%%%%%%%%%%%%%%%%%%%%% PLOT RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:dis
    [val_min ind_min]=min(MAXI_vs(i,:));
    [val_max ind_max]=max(MINI_vs(i,:));
    sup(i)=val_min;
    inf(i)=val_max;
end
%
x2 = [x'; flipud(x')];
inbetween = [inf; flipud(sup)];
%%%%%%%%%%%%%%%%%%%%%
% Plot only points within the given credible interval
CI_density_limit=CI_density;
for i = 1:dis
    for j = 1:dis
        if(  y(j) < inf(i) ||  y(j) > sup(i) )
            CI_density_limit(i,j) = 0;
        end
    end
end
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
imagesc(vs_edge,depth_edge,CI_density_limitN) %pdf plot
hold on
plot(vs_mode_n0,x,'k','LineWidth',2); % mode solution
hold on
plot(AV,x,'r','LineWidth',2); %average solution
title('Shear Wave Velocity Plot');
% Create ylabel
ylabel('Depth (m)');
% Create xlabel
xlabel('Vs (m/s)');
colormap jet
colorbar
caxis([0 0.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:dis
    [val_min_vp ind_min_vp]=min(MAXI_vp(i,:));
    [val_max_vp ind_max_vp]=max(MINI_vp(i,:));
    sup_vp(i)=val_min_vp;
    inf_vp(i)=val_max_vp;
end
%
x2 = [x'; flipud(x')];
inbetween_vp = [inf_vp; flipud(sup_vp)];
%%%%%%%%%%%%%%%%%%%%%%
% Plot only points within the given credible interval
CI_density_vp_limit=CI_density_vp;
for i = 1:dis
    for j = 1:dis
        if(  yp(j) < inf_vp(i) ||  yp(j) > sup_vp(i) )
            CI_density_vp_limit(i,j) = 0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalising density plot
sizeCI = size(CI_density_vp_limit);
len = sizeCI(1,1);
wid = sizeCI(1,2);
CI_density_vp_limitN = zeros(len,wid);
for i = 1:len
    m = sum(CI_density_vp_limit(i,:));
    for j = 1:wid
        CI_density_vp_limitN(i,j) = CI_density_vp_limit(i,j)/m;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find mode solution 
%calculating modal solution from CI_density plot
    mode_vp = zeros((dis-1),(dis-1));
    for i=1:(dis-1)
        [M,I]=max(CI_density_vp_limit(i,:));
        mode_vp(i,I)=1;
    end
    vp_mode=zeros((dis-1),(dis-1));
    for i=1:(dis-1)
    vp_mode(i,:) = mode_vp(i,:).* vp_edge;
    end   
%removing the zeros in the solution
    vp_mode_n0 = nan((dis-1),1);
    for i=1:(dis-1)
    vp_mode_n0(i,1) = max(vp_mode(i,:));
    end
    vp_mode_n0(dis,1) = vp_mode_n0(dis-1,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Vp plot (Normalised)
figure
imagesc(vp_edge,depth_edge,CI_density_vp_limitN) %pdf plot
hold on
plot(vp_mode_n0,x,'r','LineWidth',2); % mode solution
hold on
plot(priors.vpmean,x,'black','LineWidth',2); %input mean vp
hold on 
plot((priors.vpmean+priors.vpsd),x,'g','LineWidth',2); %input sd vp
hold on 
plot((priors.vpmean-priors.vpsd),x,'g','LineWidth',2); %input sd vp
title('P Wave Velocity Plot');
% Create ylabel
ylabel('Depth (m)');
% Create xlabel
xlabel('Vp (m/s)');
%set(gca,'Ydir','normal')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%VpVs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalising density plot
sizeCI = size(CI_density_VpVs);
len = sizeCI(1,1);
wid = sizeCI(1,2);
CI_density_VpVsN = zeros(len,wid);
for i = 1:len
    m = sum(CI_density_VpVs(i,:));
    for j = 1:wid
        CI_density_VpVsN(i,j) = CI_density_VpVs(i,j)/m;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find mode solution 
%calculating modal solution from CI_density plot
    mode_VpVs = zeros((dis-1),(dis-1));
    for i=1:(dis-1)
        [M,I]=max(CI_density_VpVs(i,:));
        mode_VpVs(i,I)=1;
    end
    VpVs_mode=zeros((dis-1),(dis-1));
    for i=1:(dis-1)
    VpVs_mode(i,:) = mode_VpVs(i,:).* VpVs_edge;
    end   
%removing the zeros in the solution
    VpVs_mode_n0 = nan((dis-1),1);
    for i=1:(dis-1)
    VpVs_mode_n0(i,1) = max(VpVs_mode(i,:));
    end
    VpVs_mode_n0(dis,1) = VpVs_mode_n0(dis-1,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % VpVs plot (Normalised) 
% figure
% imagesc(VpVs_edge,depth_edge,CI_density_VpVsN) %pdf plot
% hold on
% plot(VpVs_mode_n0,x,'k','LineWidth',2); % mode solution
% title('VpVs PDF');
% % Create ylabel
% ylabel('Depth (m)');
% % Create xlabel
% xlabel('VpVs modulus');
% colormap jet
% colorbar
% caxis([0 0.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%normalising density plot
sizeCI = size(CI_density_PR);
len = sizeCI(1,1);
wid = sizeCI(1,2);
CI_density_PRN = zeros(len,wid);
for i = 1:len
    m = sum(CI_density_PR(i,:));
    for j = 1:wid
        CI_density_PRN(i,j) = CI_density_PR(i,j)/m;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%find mode solution 
%calculating modal solution from CI_density plot
    mode_PR = zeros((dis-1),(dis-1));
    for i=1:(dis-1)
        [M,I]=max(CI_density_PR(i,:));
        mode_PR(i,I)=1;
    end
    PR_mode=zeros((dis-1),(dis-1));
    for i=1:(dis-1)
    PR_mode(i,:) = mode_PR(i,:).* PR_edge;
    end   
%removing the zeros in the solution
    PR_mode_n0 = nan((dis-1),1);
    for i=1:(dis-1)
    PR_mode_n0(i,1) = max(PR_mode(i,:));
    end
    PR_mode_n0(dis,1) = PR_mode_n0(dis-1,1);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % PR plot (Normalised) 
% figure
% imagesc(PR_edge,depth_edge,CI_density_PRN) %pdf plot
% hold on
% plot(PR_mode_n0,x,'k','LineWidth',2); % mode solution
% title('PR PDF');
% % Create ylabel
% ylabel('Depth (m)');
% % Create xlabel
% xlabel('PR modulus');
% colormap jet
% colorbar
% caxis([0 0.4])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot best model modal dispersion curves and observed data 
if running_mode == 1
figure;
plot( freq, data, '*');
hold on;
plot( freq, forward_model_best(:,2:4), 'r');
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
save('output_results.mat') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%      THE END                              %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

