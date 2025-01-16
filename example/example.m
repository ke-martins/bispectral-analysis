% Example of script to compute spectral and bispectral quantities on a continuous free surface elevation timeseries
% Original dataset from the Anglet (France) 2018 dataset:
% [1] Mouragues, A., Bonneton, P., Castelle, B., Marieu, V., Jak McCarroll, R., Rodriguez‐Padilla, I., ... & Sous, D. (2020). High‐energy surf zone 
%     currents and headland rips at a geologically constrained mesotidal beach. Journal of Geophysical Research: Oceans, 125(10), e2020JC016259.
% [2] Martins, K., Bonneton, P., Lannes, D., & Michallet, H. (2021). Relation between orbital velocities, pressure, and surface elevation in nonlinear 
%     nearshore water waves. Journal of Physical Oceanography, 51(11), 3539-3556.
clear all; close all; clc

% Libraries
addpath('../')

% Loading data
% data = load('AST_Anglet_13OCT2018_CASE_A.mat');
data = load('AST_Anglet_13OCT2018_CASE_B.mat');

% Computing PSDs
p.wind     = 'hann'; % Windowing applied to each block of data
p.overlap  = 75;     % Overlap in percentage
p.nfft     = 256;    % For FFT, in s
psd_zeta  = fun_compute_spectrum( data.zeta , data.sf , p.nfft*data.sf , p.overlap , p.wind );

% Computing bispectrum
bis_AST  = fun_compute_bispectrum_H1982( data.zeta , data.sf , p.nfft*data.sf , p.overlap , 'rectangular' );
% bis_AST  = fun_compute_bispectrum_H2001( data.zeta , data.sf , p.nfft*data.sf , p.overlap , 'rectangular' );

% Computing root-mean square wavenumber
bis_AST.k_rms      = fun_compute_krms( data.h0 , bis_AST.f , bis_AST.P , bis_AST.B );
bis_AST.k_rms_info = 'Boussinesq estimate of root-mean square wavenumber (Herbers et al., 2002)';
bis_AST.k_rms_4th  = fun_compute_krms( data.h0 , bis_AST.f , bis_AST.P , bis_AST.B , 'second' );
bis_AST.k_rms_4th_info = 'Boussinesq estimate of root-mean square wavenumber (Herbers et al., 2002), with fourth order frequency dispersion effects';

%% Plot - PSD and wave phase velocity spectrum
% Figure
scrsz = get(0,'ScreenSize'); fig1 = figure(1); 
set(fig1,'Position',[500 350 scrsz(3)*0.35 scrsz(4)*0.40],'PaperUnits', 'Centimeters', 'PaperPosition', [0 0 20 10],'color','w');
set(0,'defaultAxesFontSize',8)

% K_L -- Approximation by Guo (2002) of the linear wave dispersion
nmid = (length(bis_AST.f)-1)/2 + 1; % Middle frequency (f = 0)
kL   = (2*pi*bis_AST.f(nmid:end)).^2/9.81 .* (1-exp(-((2*pi*bis_AST.f(nmid:end))*sqrt(data.h0/9.81)).^(5/2))).^(-2/5);
 
% Wave phase speed spectra
h(1) = subplot(2,2,1); hl_1 = nan(1,3);
hl_1(1) = plot( bis_AST.f , 0*bis_AST.f + sqrt(9.81*data.h0) ,'k--','LineWidth',0.5); hold on, grid on, box on
hl_1(2) = plot( bis_AST.f(nmid:end) , 2*pi*bis_AST.f(nmid:end) ./ kL , 'r', 'LineWidth', 1);
hl_1(3) = plot( bis_AST.f , 2*pi*bis_AST.f ./ bis_AST.k_rms , 'ko', 'markersize', 2., 'LineWidth', 0.5 ); hold off
set(gca, 'xlim', [0 0.5]), set(gca, 'xtick', [0 0.05 0.1:0.05:1],'Fontsize',9)
set(gca, 'ylim', [0 11]), set(gca, 'ytick', [0:2:30],'Fontsize',9)
ylabel( '$c(f) = 2\pi f/\kappa$ \,[m/s]', 'Interpreter', 'Latex', 'Fontsize', 11); 
set(gca,'TickDir','out');
text(0.012,10.45,'(a)','Fontsize',9,'FontWeight','bold')
text(0.035,0.9+3*1,'$h = 9.5\,$m','Fontsize',8,'Interpreter','Latex')
text(0.035,0.9+2*1,'$H_{m0} = 3.3\,$m','Fontsize',8,'Interpreter','Latex')
text(0.035,0.9+1,'$\mu = 0.24$','Fontsize',8,'Interpreter','Latex')
text(0.035,0.9,'$U_r = 0.7$','Fontsize',8,'Interpreter','Latex')
xtickangle(0)
% Legend
leg = legend( hl_1 , '$\sqrt{gh_0}$', '$2\pi f/\kappa_L$','$2\pi f/\kappa_{rms}$','Location','South','Interpreter','Latex'); leg.ItemTokenSize = [16,16];
% leg = legend( hl_1 ,'$2\pi f/\kappa_L$','$2\pi f/\kappa_B$','$2\pi f/\kappa_{rms}$','Location','NorthEast'); leg.ItemTokenSize = [16,16,16];
leg_pos = get(leg,'Position'); leg_pos(1) = leg_pos(1)-0.02; leg_pos(2) = leg_pos(2)+0.04; set(leg,'Position',leg_pos);
set(leg,'Fontsize',10), legend boxoff

% PSD - A3
h(2) = subplot(2,2,3);
semilogy( psd_zeta.f, psd_zeta.E , 'k' , 'LineWidth' , 1 ); hold on, grid on, box on
set(gca, 'xlim', [0 0.5]), set(gca, 'xtick', [0 0.05 0.1:0.05:1],'Fontsize',9)
set(gca, 'ylim', [10^-3 10^2]), set(gca, 'ytick', [10^-7 10^-6 10^-5 10^-4 10^-3 10^-2 10^-1 10^0 10^1 10^2],'Fontsize',9)
xlabel( '$f$ [Hz]', 'Interpreter', 'Latex', 'Fontsize', 11)
ylabel( '$E(f)$ [m$^2$/Hz]','Interpreter','Latex','Fontsize',11)
set(gca,'TickDir','out');
text(0.0105,4.25*10^1,'(b)','Fontsize',9,'FontWeight','bold')

% Add confidence levels 
coef_ulimit = psd_zeta.CI(2); coef_llimit = psd_zeta.CI(1);
f_CI = 0.075; psd_CI = 0.007;
std_upper_limit = psd_CI*coef_ulimit;
std_lower_limit = psd_CI*coef_llimit;
errorbar( f_CI , psd_CI , abs(psd_CI-std_lower_limit) , abs(psd_CI-std_upper_limit) , 'ko' , 'markerfacecolor' , 'k' , 'capsize' , 4 , 'markersize' , 2 , 'linewidth', 0.5)
text( 0.05 , 0.015 , '95% C.I.','Fontsize',7), pause(0.5)
% Legend
leg2 = legend( '$\zeta_{\rm AST}$','Location','North','Interpreter','Latex');
leg_pos = get(leg2,'Position'); leg_pos(1) = leg_pos(1)+0.04; leg_pos(2) = leg_pos(2)+0.; set(leg2,'Position',leg_pos);
set(leg2,'Fontsize',10,'NumColumns',3), legend boxoff; leg2.ItemTokenSize = [16,16,16,16];
xtickangle(0)

% Bicoherence (showing only a fraction, due to bispectrum symmetrical properties)
h(3) = subplot(2,2,2);
fid = find(and(bis_AST.f>=0,bis_AST.f<=0.5));
f_plot = bis_AST.f(fid);
bic = bis_AST.Bic(fid,fid);
ifc = find(f_plot>0.25,1,'first');
for ff = 2:ifc
  bic(ff+1:end,ff) = NaN;
end
for ff = ifc:2*ifc
  if ff <= numel(fid)
    bic(2*ifc-ff+1:end,ff) = NaN;
  end
end
[~,hcon,~] = contourf( f_plot , f_plot , bic.^2 ); hcon.LineWidth = 0.25; hold on, grid on
set(gca, 'xlim', [0 0.5]), set(gca, 'xtick', [0 0.1:0.1:1],'Fontsize',9)
set(gca, 'ylim', [0 0.275001]), set(gca, 'ytick', [0 0.05:0.05:1],'Fontsize',9)
xlabel( '$f$ [Hz]', 'Interpreter', 'Latex', 'Fontsize', 11)
ylabel( '$f$ [Hz]', 'Interpreter', 'Latex', 'Fontsize',11)
set(gca,'TickDir','out'); %set(gca, 'YScale', 'log');set(gca, 'XScale', 'log');
text(0.015 , 0.262 ,'(c)','Fontsize',9,'FontWeight','bold')
x = [0.06 0.08];    % adjust length and location of arrow 
y = [0.13 0.09];      % adjust hieght and width of arrow
annotation('textarrow',[0.675 0.6825]-0.03,...
  [0.415 0.35]+0.02,'String',' ($f_p,f_p$) ','FontSize',10,'Linewidth',0.75,'interpreter','latex','HeadWidth',4,'HeadLength',4)
xtickangle(0)

% Color bar
cmap = hot(128); colormap(flipud(cmap)); gca_pos = get(gca,'Position');
caxis([0 0.6]); pause(0.5)
hc = colorbar('XTick',0:0.1:0.6,'Location','northoutside');  
clabel = ylabel(hc,'$b^2(f,f)$ [-]','interpreter','latex','fontsize',10);
set(hc,'LineWidth',0.75); cl_pos = get(hc,'Position');
cl_pos(1) = cl_pos(1)+0.09;
cl_pos(2) = cl_pos(2)+0.03;
cl_pos(3) = 0.80*cl_pos(3);
cl_pos(4) = 0.5*cl_pos(4);
set(hc,'Position',cl_pos), pause(0.5)
set(gca,'Position',gca_pos)

% Positions
set(h(1),'Position',[0.08 0.60 0.41 0.375])
set(h(2),'Position',[0.08 0.14 0.41 0.375])
set(h(3),'Position',[0.595 0.14 0.395 0.69])

% Saving
% print(fig1,'-depsc','-r300','PSD_krms_and_Bic_H2001_CASE_A')
% print(fig1,'-depsc','-r300','PSD_krms_and_Bic_H2001_CASE_B')
print(fig1,'-depsc','-r300','PSD_krms_and_Bic_H1982_CASE_B')
