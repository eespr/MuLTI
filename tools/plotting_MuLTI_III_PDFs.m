close all

%load('MuLTI_III_output_data.mat')

solution = 1; % 0 = average ; 1 = mode ; for Vp Vs and density plots only 
% the mode solution will be plotted for shear modulus, bulk modulus, Vp:Vs ratio and Poisson's ratio. 

load('cmap_b.mat') % load colourbar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
hold on
plot(Vs_true,x,'r','LineWidth',2); %average solution
%title('Shear Wave Velocity PDF');
% Create ylabel
ylabel('Depth (m)');
% Create xlabel
xlabel('Vs (m/s)');
set(gca,'Ydir','reverse')
set(gca,'fontsize',10);
colormap(cm)
%colorbar
caxis([0 0.1])
xlim([500 2500])

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
    %title('P Wave Velocity PDF');
    % Create ylabel
    %ylabel('Depth (m)');
    % Create xlabel
    xlabel('Vp (m/s)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',10);
    colormap(cm)
    %colorbar
    caxis([0 0.1])
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
    %title('Density PDF');
    % Create ylabel
    %ylabel('Depth (m)');
    % Create xlabel
    xlabel('Density (g/cc)');
    set(gca,'Ydir','reverse')
    set(gca,'fontsize',10);
    colormap(cm)
    %colorbar
    caxis([0 0.1])
    xlim([0.5 1.05])
    ylim([0 50])

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot best model modal dispersion curves and observed data 
if running_mode == 1
%PDF ALL modes
figure
contourf(freq_edge,vs_edge,CI_density_ALL,1000,'edgecolor','none') %pdf plot
hold on;
errorbar(freq, data, fitting_error,'ok','MarkerFaceColor','k','linewidth',2) % observed data
% Create ylabel
ylabel('Phase velocity (m/s)');
% Create xlabel
xlabel('Frequency (Hz)');
colormap(cm)
colorbar
set(gca,'YDir','normal')
caxis([0 50000])
%PDF misfit
figure
contourf(freq_edge,misfit_edge,CI_density_misfit,1000,'edgecolor','none') %pdf plot
hold on
plot(freq, data.*0.05,'k','linewidth',4);
% Create ylabel
ylabel('Misfit (m/s)');
% Create xlabel
xlabel('Frequency (Hz)');
colormap(cm)
colorbar
set(gca,'YDir','normal')
caxis([0 50000])
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
