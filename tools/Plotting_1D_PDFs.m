close all

%%%%%%%%%%%%%%%%%%%%%%INPUT PLOTTING PARAMETERS%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('output_data.mat') % load results
load('cmap.mat') % load colourbar for plotting PDF's

multi = 3; % multi = 0 non constrained inversion where num_layers = 1, 
           % multi = 1 constrained MuLTI inversion, 
           % multi = 2 MuLTI II inversion non constrainted, 
           % multi = 2.5 MuLTI II inversion constrainted, 
           % multi = 3 MuLTI III inversion (constrained or non-constrained) Vp, Vs and density PDF's, 
           % multi = 6 MuLTI III inversion, plot all PDF's including Vp,
           % Vs, density, shear mod, bulk mod, Vp/Vs and PR
           %(all constrained cases are assuming num_layers = 3)

solution = 1; % 0 = average ; 1 = mode (choose solution to plot on Vp, Vs and density PDF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOT 1D RESULTS%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VS plot (Normalised) 
figure
subplot(2,3,[1 4])
contourf(vs_edge,depth_edge,CI_density_limitN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
hold on
    if solution == 1
    plot(vs_mode_n0,x,'k','LineWidth',2); % mode solution
    else
    plot(AV,x,'k','LineWidth',2); % average solution
    end
title('Shear Wave Velocity PDF');
% Create ylabel
ylabel('Depth (m)');
% Create xlabel
xlabel('Vs (m/s)');
set(gca,'Ydir','reverse')
set(gca,'fontsize',14);
colormap(cm)
colorbar
caxis([0 0.3])
xlim([500 2500])

if multi == 0 %multi non constrained 
    subplot(2,3,[2 5])
    plot((priors.vp.*ones(dis)),x,'r','LineWidth',2); %input sd vp
    title('P Wave Velocity PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Vp (m/s)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([1500 4000])
    ylim([0 50])
    
    subplot(2,3,[3 6])
    plot((priors.density.*ones(dis)),x,'r','LineWidth',2); %input sd density
    title('Density PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Density (g/cc)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([0.5 1.05])
    ylim([0 50])
    
elseif multi == 1 %multi constrained 
   % calculate fixed layer constraints for Vp and density 
    [val1, idx1] = min(abs(x - priors.layer_depths(1))); 
    [val2, idx2] = min(abs(x - priors.layer_depths(2))); 
    vp = zeros(dis,1);
    density = zeros(dis,1);
    vp(1:idx1,1) = priors.vp(1);
    vp(idx1+1:idx2,1) = priors.vp(2);
    vp(idx2+1:dis,1) = priors.vp(3);
    density(1:idx1,1) = priors.density(1);
    density(idx1+1:idx2,1) = priors.density(2);
    density(idx2+1:dis,1) = priors.density(3);
    
    subplot(2,3,[2 5])
    plot(vp.',x,'r','LineWidth',2); %input sd vp
    title('P Wave Velocity PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Vp (m/s)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([1500 4000])
    ylim([0 50])
    
    subplot(2,3,[3 6])
    plot(density.',x,'r','LineWidth',2); %input sd density
    title('Density PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Density (g/cc)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([0.5 1.05])
    ylim([0 50])
    
elseif multi == 2 % MuLTI II inversion non constrainted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vp plot (Normalised) 
    subplot(2,3,[2 5])
    contourf(vp_edge,depth_edge,CI_density_vp_limitN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
        if solution == 1
        plot(vp_mode_n0,x,'k','LineWidth',2); % mode solution
        else
        plot(AV_vp,x,'k','LineWidth',2); % average solution
        end
    hold on
    plot(priors.vpmean,x,'r','LineWidth',2); %input mean vp
    hold on 
    plot((priors.vpmean+priors.vpsd),x,'g','LineWidth',2); %input sd vp
    hold on 
    plot((priors.vpmean-priors.vpsd),x,'g','LineWidth',2); %input sd vp
    title('P Wave Velocity PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Vp (m/s)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    %colorbar
    caxis([0 0.3])
    xlim([1500 4000])
    ylim([0 50])
    
    subplot(2,3,[3 6])
    plot((priors.density.*ones(dis)),x,'r','LineWidth',2); %input sd density
    title('Density PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Density (g/cc)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([0.5 1.05])
    ylim([0 50])
    
elseif multi == 2.5 % MuLTI II inversion constrainted
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vp plot (Normalised) 
    subplot(2,3,[2 5])
    contourf(vp_edge,depth_edge,CI_density_vp_limitN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
        if solution == 1
        plot(vp_mode_n0,x,'k','LineWidth',2); % mode solution
        else
        plot(AV_vp,x,'k','LineWidth',2); % average solution
        end
    hold on
    plot(priors.vpmean,x,'r','LineWidth',2); %input mean vp
    hold on 
    plot((priors.vpmean+priors.vpsd),x,'g','LineWidth',2); %input sd vp
    hold on 
    plot((priors.vpmean-priors.vpsd),x,'g','LineWidth',2); %input sd vp
    title('P Wave Velocity PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Vp (m/s)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    %colorbar
    caxis([0 0.3])
    xlim([1500 4000])
    ylim([0 50])
    
    % calculate fixed layer constraints for Vp and density
    [val1, idx1] = min(abs(x - priors.layer_depths(1))); 
    [val2, idx2] = min(abs(x - priors.layer_depths(2))); 
    density = zeros(dis,1);
    density(1:idx1,1) = priors.density(1);
    density(idx1+1:idx2,1) = priors.density(2);
    density(idx2+1:dis,1) = priors.density(3);
    
    subplot(2,3,[3 6])
    plot(density.',x,'r','LineWidth',2); %input sd density
    title('Density PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Density (g/cc)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([0.5 1.05])
    ylim([0 50])
    
else %multi == 3 || multi == 6 multi III and joint elastic PDFs

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Vp plot (Normalised) 
    subplot(2,3,[2 5])
    contourf(vp_edge,depth_edge,CI_density_vp_limitN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
        if solution == 1
        plot(vp_mode_n0,x,'k','LineWidth',2); % mode solution
        else
        plot(AV_vp,x,'k','LineWidth',2); % average solution
        end
    hold on
    plot(priors.vpmean,x,'r','LineWidth',2); %input mean vp
    hold on 
    plot((priors.vpmean+priors.vpsd),x,'g','LineWidth',2); %input sd vp
    hold on 
    plot((priors.vpmean-priors.vpsd),x,'g','LineWidth',2); %input sd vp
    title('P Wave Velocity PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Vp (m/s)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    %colorbar
    caxis([0 0.3])
    xlim([1500 4000])
    ylim([0 50])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Density plot (Normalised) 
    subplot(2,3,[3 6])
    contourf(den_edge,depth_edge,CI_density_den_limitN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
        if solution == 1
        plot(den_mode_n0,x,'k','LineWidth',2); % mode solution
        else
        plot(AV_den,x,'k','LineWidth',2); % average solution
        end
    hold on
    plot(priors.denmean,x,'r','LineWidth',2); %input mean denisty
    hold on 
    plot((priors.denmean+priors.densd),x,'g','LineWidth',2); %input sd density
    hold on 
    plot((priors.denmean-priors.densd),x,'g','LineWidth',2); %input sd density
    title('Density PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Density (g/cc)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    %colorbar
    caxis([0 0.3])
    xlim([0.5 1.05])
    ylim([0 50])
end

if multi == 6 %joint elastic PDFs
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shear plot (Normalised) 
    figure
    contourf(shear_edge,depth_edge,CI_density_shearN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
    plot(shear_mode_n0,x,'k','LineWidth',2); % mode solution
    title('Shear PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Shear modulus (10^x)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    colorbar
    caxis([0 0.3])
    %set(gca,'XScale','log') 
    %set(gca,'XMinorTick','on') 
    set(gca,'Ydir','reverse')
    xlim([5 7])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % bulk plot (Normalised) 
    figure
    contourf(bulk_edge,depth_edge,CI_density_bulkN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
    plot(bulk_mode_n0,x,'k','LineWidth',2); % mode solution
    title('bulk PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('bulk modulus (10^x)');
    colormap(cm)
    colorbar
    caxis([0 0.3])
    %set(gca,'XScale','log') 
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    xlim([6 7.2])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % VpVs plot (Normalised) 
    figure
    contourf(VpVs_edge,depth_edge,CI_density_VpVsN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
    plot(VpVs_mode_n0,x,'k','LineWidth',2); % mode solution
    title('VpVs PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('VpVs');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    colorbar
    caxis([0 0.3])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PR plot (Normalised) 
    figure
    contourf(PR_edge,depth_edge,CI_density_PRN(1:99,1:99), 1000,'edgecolor','none') %pdf plot
    hold on
    plot(PR_mode_n0,x,'k','LineWidth',2); % mode solution
    title('Poissons ratio PDF');
    % Create ylabel
    ylabel('Depth (m)');
    % Create xlabel
    xlabel('Poissons ratio');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',14);
    colormap(cm)
    colorbar
    caxis([0 0.3])
end

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot best model modal dispersion curves and observed data 
if running_mode == 1
figure;
plot( freq, forward_model_best(:,2:4),'k','linewidth',4);
hold on;
plot( freq, data, '*','markersize',8);
% Create ylabel
ylabel('Phase velocity (m/s)');
% Create xlabel
xlabel('Frequency (Hz)');
set(gca,'fontsize',14);
end


%%%%%%%%%%%%%%% Plot Statistics of the chain %%%%%%%%%%%%%%%%
figure
subplot(2,1,1)
semilogy(cov);
hold on
%line([burn_in burn_in],[0 cov(1)],'LineWidth',4,'Color',[1 0 0]);
xlabel('iterations','FontSize',14)
set(gca,'fontsize',14);

title('Data Misfit','FontSize',16)
subplot(2,1,2);
line([burn_in burn_in],[0 priors.npt_max+num_layers],'LineWidth',4,'Color',[1 0 0]);
hold on
plot(nnuclei)
xlabel('iterations','FontSize',14)
title('Number of nuclei','FontSize',16)
set(gca,'fontsize',14);


%%%%%%%%%%%%%% Plot Marginal posteriors %%%%%%%%%%%%%%%%%%%

% Plot histogram on change points
figure
subplot(2,1,1)
hist(change_points(1:bb),500)
title('Probability of change points','FontSize',14)
set(gca,'fontsize',14);

% Plot histogram on number of nuclei
subplot(2,1,2);
bar(nnucleihist)
title('Posterior Distribution on number of nuclei ','FontSize',14)
set(gca,'fontsize',14);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
