%% Paramètres
clear all; close all; clc;

tfin = 2000.0;     % temps final de la simulation
hbar = 1.0;	       % constante de planck normalise
mass = 1.0;        % mass de la particule
xL = -200.0;       % position bord de gauche
xR = 200.0;        % position bord de droite
omega = 0.004;     % parametre pour definir le potentiel
Delta = 0.0;       % parametre pour definir le potentiel
x0 = 0.0;          % centre du paquet d'onde initial
n = 16;            % nombre d onde du paquet d'onde initial
sigma_norm = 0.05; % ecart-type du paquet d'onde initial en unites de (xR-xL)

% Detecteur:
t_detect = -1.0;  % temps de detection (un temps negatif implique pas de detection)
xda = 130.0;      % limite de droite du senseur du detecteur
xdb = 230.0;      % limite de gauche du senseur du detecteur

% Numerique :
dt = 1.0;      % pas de temps
Ninters = 256; % nombre d'intervalles du maillage

% Simulations
    filename2  = ['Dt_',num2str(dt)];
    dt_loc = dt; fname2 = filename2;
    Delta_loc = Delta;
    t_d_loc = t_detect;
    x0_loc = x0;
    writeConfig;
    disp('Exercice8_Kervyn_LeMeur configuration.in');   
    system('Exercice8_Kervyn_LeMeur configuration.in'); 
    disp('Done.')

% chargement des résultats
    data = load([filename2,'_observables.out']);
    t = data(:,1);
    Pgauche = data(:,2);
    Pdroite = data(:,3);
    Ptot = data(:,4);
    E = data(:,5);
    xmoy = data(:,6);
    x2moy = data(:,7);
    pmoy = data(:,8);
    p2moy = data(:,9);
    DeltaxDeltap = data(:,10); %(Delta x)(Delta p)
    Deltax = data(:,11);
    Deltap = data(:,12);
    data = load([filename2,'_potential.out']);
    x = data(:,1);
    V = data(:,2);
    tpsi2 = load([filename2,'_psi2.out']);
    [nt,nx1]=size(tpsi2);
    psi2 = tpsi2(2:nt,2:nx1);

%% Figures du prof
ViewFormat;

figure
    plot(t,Pgauche, t,Pdroite,t,Ptot,'linewidth',lw)
    set(gca,'fontsize',fs)
    grid on
    xlabel('$t$')
    ylabel('$P$')
    legend('$P_{x<x_a}$','$P_{x>x_b}$','$P_{tot}$','Location','Best')
    %figure('Name',['Analyse de ' fichier], 'unit', 'norm', 'pos', [1,.1,.8,.3])

%subplot(2,2,1)
figure
    plot(x,V,'linewidth',lw)
    set(gca,'fontsize',fs)
    hold on
    plot(x([1,end]),E(1)*ones(1,2),'--','linewidth',lw)
    grid on
    xlabel('$x$')
    ylabel('$V$')
    legend('$V(x)$','$E$','Location','Best')
SaveIMG("PotentielHarmonique")

figure
    [X,T] = meshgrid(x,t);
    pcolor(X,T,psi2); % voir aussi contour, contourf
    set(gca,'fontsize',fs)
    shading interp
    colormap jet % il existe d'autres colormaps (help colormap)
    c = colorbar;
    xlabel('$x$')
    ylabel('$t$')
    ylabel(c,'$|\psi|^2$','Interpreter','Latex','FontSize',fs)
SaveIMG("densitéProbaHarmonique")

figure
    plot(t,E,'linewidth',lw)
    set(gca,'fontsize',fs)
    grid on
    xlabel('$t$')
    ylabel('$E$')

%% Etude de convergence en dt
clear all; close all; clc;

tfin = 2000.0;       % temps final de la simulation
hbar = 1.0;	      % constante de planck normalise
mass = 1.0;        % mass de la particule
xL = -200.0;       % position bord de gauche
xR = 200.0;        % position bord de droite
omega = 0.004;      % parametre pour definir le potentiel
Delta = 0.0;       % parametre pour definir le potentiel
x0 = 0.0;          % centre du paquet d'onde initial
n = 16;            % nombre d onde du paquet d'onde initial
sigma_norm = 0.05; % ecart-type du paquet d'onde initial en unites de (xR-xL)

% Detecteur:
t_detect = -1.0;  % temps de detection (un temps negatif implique pas de detection)
xda = 130.0;      % limite de droite du senseur du detecteur
xdb = 230.0;      % limite de gauche du senseur du detecteur

% Numerique :
nsimul = 10;
Nsteps = round(linspace(2000,20000,nsimul));
dt = tfin./Nsteps; % pas de temps
Ninters = 256; % nombre d'intervalles du maillage

% Simulations
esp_pos    = zeros(1,nsimul); %espérances
esp_qdm    = zeros(1,nsimul);
esp_deltax = zeros(1,nsimul);
esp_deltap = zeros(1,nsimul);
tf         = zeros(1,nsimul);
for i=1:nsimul
    filename2  = ['Dt_',num2str(dt(i))];
    dt_loc = dt(i); fname2 = filename2;
    x0_loc = x0;
    Delta_loc = Delta;
    t_d_loc = t_detect;
    writeConfig;
    disp('Exercice8_Kervyn_LeMeur configuration.in');   
    system('Exercice8_Kervyn_LeMeur configuration.in'); 
    disp('Done.')
    
    % chargement des résultats
    data = load([filename2,'_observables.out']);
    t = data(:,1);
    Pgauche = data(:,2);
    Pdroite = data(:,3);
    Ptot = data(:,4);
    E = data(:,5);
    xmoy = data(:,6);
    x2moy = data(:,7);
    pmoy = data(:,8);
    p2moy = data(:,9);
    DeltaxDeltap = data(:,10); %(Delta x)(Delta p)
    Deltax = data(:,11);
    Deltap = data(:,12);
    data = load([filename2,'_potential.out']);
    x = data(:,1);
    V = data(:,2);
    tpsi2 = load([filename2,'_psi2.out']);
    [nt,nx1]=size(tpsi2);
    psi2 = tpsi2(2:nt,2:nx1);
    
    tf(1,i)      = t(end);
    
    % espérance sur la position 
    esp_pos(1,i) = (xmoy(end)-xmoy(end-1))/(t(end)-t(end-1))*tfin+(xmoy(end)-(xmoy(end)-xmoy(end-1))/(t(end)-t(end-1))*t(end));
    % espérance sur la qdm
    esp_qdm(1,i) = (pmoy(end)-pmoy(end-1))/(t(end)-t(end-1))*tfin+(pmoy(end)-(pmoy(end)-pmoy(end-1))/(t(end)-t(end-1))*t(end));
    % espérance sur Delta x
    esp_deltax(1,i) = (Deltax(end)-Deltax(end-1))/(t(end)-t(end-1))*tfin+(Deltax(end)-(Deltax(end)-Deltax(end-1))/(t(end)-t(end-1))*t(end));
    % espérance sur Delta p
    esp_deltap(1,i) = (Deltap(end)-Deltap(end-1))/(t(end)-t(end-1))*tfin+(Deltap(end)-(Deltap(end)-Deltap(end-1))/(t(end)-t(end-1))*t(end));
end

%% Position
ViewFormat;
figure('Name',"Convergence espérance X")
    plot(dt.^2,esp_pos,'+','LineWidth',lw);
    
    % fit
    hold on
    P = polyfit(dt.^2,esp_pos,1); Pos_AS = P(2) 
    z = polyval(P, dt.^2);
    plot(dt.^2,z,'--','Linewidth',lw,'HandleVisibility','off');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$(\Delta t)^2$'); ylabel('$\langle x \rangle_f$');
SaveIMG("ConvergenceEsperanceX");

figure('Name',"Convergence LOGLOG espérance X")
    loglog(dt,abs((esp_pos-Pos_AS)/Pos_AS)*100,'+','Linewidth',lw,'HandleVisibility','off');
    hold on
        P = polyfit(log(dt),log(abs((esp_pos-Pos_AS)/Pos_AS)*100),1); 
        z  = polyval(P, log(dt));
        loglog(dt,exp(z),'--','Linewidth',1);
        legendStrings = string(P(1));
        leg = legend(legendStrings,'Location','northwest');
        title(leg, 'linear fit: slope')
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\Delta t$'); ylabel('$\varepsilon_{\langle x \rangle}$ [$\%$]');
SaveIMG("ConvergenceLOGLOGX");

%% Qdm
figure('Name',"Convergence espérance P")
    plot(dt.^2,esp_qdm,'+','LineWidth',lw);
    
    % fit
    hold on
    P = polyfit(dt.^2,esp_qdm,1); Qdm_AS = P(2)
    z = polyval(P, dt.^2);
    plot(dt.^2,z,'--','Linewidth',lw,'HandleVisibility','off');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$(\Delta t)^2$'); ylabel('$\langle p \rangle_f$');
SaveIMG("ConvergenceEsperanceP");

figure('Name',"Convergence LOGLOG espérance P")
    loglog(dt,abs((esp_qdm-Qdm_AS)/Qdm_AS)*100,'+','Linewidth',lw,'HandleVisibility','off');
    hold on
        P = polyfit(log(dt),log(abs((esp_qdm-Qdm_AS)/Qdm_AS)*100),1); 
        z  = polyval(P, log(dt));
        loglog(dt,exp(z),'--','Linewidth',1);
        legendStrings = string(P(1));
        leg = legend(legendStrings,'Location','northwest');
        title(leg, 'linear fit: slope')
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\Delta t$'); ylabel('$\varepsilon _{\langle p \rangle}$ [$\%$]');
SaveIMG("ConvergenceLOGLOGP");

%% Delta x
ViewFormat;
figure('Name',"Convergence espérance Delta X")
    plot(dt.^2,esp_deltax,'+','LineWidth',lw);
    
    % fit
    hold on
    P = polyfit(dt.^2,esp_deltax,1); DX_AS = P(2)
    z = polyval(P, dt.^2);
    plot(dt.^2,z,'--','Linewidth',lw,'HandleVisibility','off');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$(\Delta t)^2$'); ylabel('$\sigma_{x,f}$');
SaveIMG("ConvergenceEsperanceDX");

figure('Name',"Convergence LOGLOG espérance Delta X")
    loglog(dt,abs((esp_deltax-DX_AS)/DX_AS)*100,'+','Linewidth',lw,'HandleVisibility','off');
    hold on
        P = polyfit(log(dt),log(abs((esp_deltax-DX_AS)/DX_AS)*100),1); 
        z  = polyval(P, log(dt));
        loglog(dt,exp(z),'--','Linewidth',1);
        legendStrings = string(P(1));
        leg = legend(legendStrings,'Location','northwest');
        title(leg, 'linear fit: slope')
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\Delta t$'); ylabel('$\varepsilon _{\sigma_x}$ [$\%$]');
SaveIMG("ConvergenceLOGLOGDX");

%% Delta p
ViewFormat;
figure('Name',"Convergence espérance Delta P")
    plot(dt.^2,esp_deltap,'+','LineWidth',lw);
    % fit
    hold on
    P = polyfit(dt.^2,esp_deltap,1); DP_AS = P(2)
    z = polyval(P, dt.^2);
    plot(dt.^2,z,'--','Linewidth',lw,'HandleVisibility','off');
    grid minor; set(gca,'fontsize',fs);
    xlabel('$(\Delta t)^2$'); ylabel('$\sigma_{p,f}$');
SaveIMG("ConvergenceEsperanceDP");

figure('Name',"Convergence LOGLOG espérance Delta P")
    loglog(dt,abs((esp_deltap-DP_AS)/DP_AS)*100,'+','Linewidth',lw,'HandleVisibility','off');
    hold on
        P = polyfit(log(dt),log(abs((esp_deltap-DP_AS)/DP_AS)*100),1); 
        z  = polyval(P, log(dt));
        loglog(dt,exp(z),'--','Linewidth',1);
        legendStrings = string(P(1));
        leg = legend(legendStrings,'Location','northwest');
        title(leg, 'linear fit: slope')
    grid minor; set(gca,'fontsize',fs);
    xlabel('$\Delta t$'); ylabel('$\varepsilon _{\sigma_p}$ [$\%$]');
SaveIMG("ConvergenceLOGLOGDP");

%% Probabilité totale
figure('Name','Probabilité totale')
    plot(t,Pgauche, t,Pdroite,t,Ptot,'linewidth',lw)
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$P$');
    legend('$P_{x<x_a}$','$P_{x>x_b}$','$P_{tot}$','Location','Best')
    ylim([-0.1 1.1]);
SaveIMG("ProbabiliteTotale");

figure('Name','Vérif Probabilité totale')
    plot(t,abs(Pgauche+Pdroite-1)*100,'linewidth',lw)
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\varepsilon_P$ [\%]');
SaveIMG("CheckProbabiliteTotale");

%% Energie moyenne
figure('Name','Energie')
    plot(t,abs((E-E(1))/E(1))*100,'linewidth',lw)
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\varepsilon_E$ [\%]');
SaveIMG("Energie");

%% Incertitude Heisenberg
ViewFormat;
figure('Name','Heisenberg')
    plot(t,DeltaxDeltap,'linewidth',lw)
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\Delta x \Delta p$');
    hold on
    yline(0.5,'r--',{'$\hbar / 2$'},'Interpreter','Latex','LabelOrientation','horizontal','Fontsize',fs,'LineWidth',lw,'LabelVerticalAlignment','bottom');
    ylim([0.4,max(DeltaxDeltap)+0.1]);
SaveIMG("Heisenberg");

%% Cas classique 
ViewFormat;
k0      = 2*pi*n/(xR-xL);
Energie = hbar^2*k0^2/(2*mass);
classiX = 1/omega * sqrt(2*Energie/mass) * sin(omega*t);
classiP = sqrt(2*Energie*mass)*cos(omega*t);

h1 = figure('Name','cas classique X');
    plot(t,classiX,'--','Linewidth',lw);
    hold on
    plot(t,xmoy,'-','Linewidth',lw);
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle x \rangle$');
    xlim([t(1) t(end)]);
    ylim([min(classiX) max(classiX)]);
    MagInset(h1, -1, [1030 1350 min(classiX) -50], [850 1550 14 50], {'NW','SW';'NE','SE'});
    set(gca,'fontsize',fs-4);
SaveIMG("ClassiX");

figure('Name','Diff cas classique X')
    plot(t,abs(classiX-xmoy)/(1/omega*sqrt(2*Energie/mass))*100,'-','Linewidth',lw);
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\varepsilon_{\langle x \rangle,C}$ [\%]');
SaveIMG("DiffClassiX");

h1 = figure('Name','cas classique P');
    plot(t,classiP,'--','Linewidth',lw);
    hold on
    plot(t,pmoy,'-','Linewidth',lw);
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle p \rangle$');
    xlim([t(1) t(end)]);
    ylim([min(classiP) max(classiP)]);
    MagInset(h1, -1, [630 955 min(classiP) -0.2], [450 1170 0.05 0.2], {'NW','SW';'NE','SE'});
    set(gca,'fontsize',fs-4);
SaveIMG("ClassiP");

figure('Name','Diff cas classique P')
    plot(t,abs(classiP-pmoy)/sqrt(2*Energie*mass)*100,'-','Linewidth',lw);
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\varepsilon_{\langle p \rangle,C}$ [\%]');
SaveIMG("DiffClassiP");

%% Effet tunnel
clear all; close all; clc;

tfin = 2000.0;       % temps final de la simulation
hbar = 1.0;	         % constante de planck normalise
mass = 1.0;          % mass de la particule
xL = -200.0;         % position bord de gauche
xR = 200.0;          % position bord de droite
omega = 0.004;       % parametre pour definir le potentiel
Delta = [40 65 70 75 80];  % parametre pour definir le potentiel
x0 = -Delta;         % centre du paquet d'onde initial
n = 16;              % nombre d onde du paquet d'onde initial
sigma_norm = 0.05;   % ecart-type du paquet d'onde initial en unites de (xR-xL)

% Detecteur:
t_detect = -1.0;     % temps de detection (un temps negatif implique pas de detection)
xda = 130.0;         % limite de droite du senseur du detecteur
xdb = 230.0;         % limite de gauche du senseur du detecteur

% Numerique :
Nsteps = 5000;
dt = tfin./Nsteps; % pas de temps
Ninters = 256; % nombre d'intervalles du maillage

k0      = 2*pi*n/(xR-xL);
Energie = hbar^2*k0^2/(2*mass);

ViewFormat;
% Simulations
for i=1:length(Delta)
    filename2  = ['Delta_',num2str(Delta(i))];
    dt_loc = dt; 
    Delta_loc = Delta(i);
    x0_loc = x0(i);
    t_d_loc = t_detect;
    fname2 = filename2;
    writeConfig;
    disp('Exercice8_Kervyn_LeMeur configuration.in');   
    system('Exercice8_Kervyn_LeMeur configuration.in'); 
    disp('Done.')
    
    % chargement des résultats
    data = load([filename2,'_observables.out']);
    t = data(:,1);
    Pgauche = data(:,2);
    Pdroite = data(:,3);
    Ptot = data(:,4);
    E = data(:,5);
    xmoy = data(:,6);
    x2moy = data(:,7);
    pmoy = data(:,8);
    p2moy = data(:,9);
    DeltaxDeltap = data(:,10); %(Delta x)(Delta p)
    Deltax = data(:,11);
    Deltap = data(:,12);
    data = load([filename2,'_potential.out']);
    x = data(:,1);
    V = data(:,2);
    tpsi2 = load([filename2,'_psi2.out']);
    [nt,nx1]=size(tpsi2);
    psi2 = tpsi2(2:nt,2:nx1);
    
    figure('Name',"Potentiel Delta"+num2str(Delta(i)))
        plot(x,V,'linewidth',lw)
        set(gca,'fontsize',fs)
        hold on
        plot(x([1,end]),E(1)*ones(1,2),'--','linewidth',lw)
     %   plot(x,Energie*ones(1,length(x)),'--','Linewidth',lw)
        grid on
        xlabel('$x$')
        ylabel('$V$')
        legend('$V(x)$','$E$','Location','Best')
    SaveIMG("Potentiel Delta"+num2str(Delta(i)));
        
    figure('Name',"Probabilité Delta"+num2str(Delta(i)))
        plot(t,Pgauche, t,Pdroite,t,Ptot,'linewidth',lw)
        xlim([0 tfin]);
        set(gca,'fontsize',fs)
        grid on
        xlabel('$t$')
        ylabel('$P$')
        legend('$P_{x<x_a}$','$P_{x>x_b}$','$P_{tot}$','Location','Best')
    SaveIMG("Probabilité Delta"+num2str(Delta(i)));
    
    figure('Name',"Densité Probabilité Delta"+num2str(Delta(i)))
        [X,T] = meshgrid(x,t);
        pcolor(X,T,psi2); % voir aussi contour, contourf
        set(gca,'fontsize',fs)
        shading interp
        colormap jet % il existe d'autres colormaps (help colormap)
        c = colorbar;
        xlabel('$x$')
        ylabel('$t$')
        ylabel(c,'$|\psi|^2$','Interpreter','Latex')
    SaveIMG("Densité Probabilité Delta"+num2str(Delta(i)));
end

%% Détection de la particule
clear all; close all; clc;

tfin = 2500.0;       % temps final de la simulation
hbar = 1.0;	         % constante de planck normalise
mass = 1.0;          % mass de la particule
xL = -200.0;         % position bord de gauche
xR = 200.0;          % position bord de droite
omega = 0.004;       % parametre pour definir le potentiel
Delta = 70;          % parametre pour definir le potentiel
x0 = -Delta;         % centre du paquet d'onde initial
n = 16;              % nombre d onde du paquet d'onde initial
sigma_norm = 0.05;   % ecart-type du paquet d'onde initial en unites de (xR-xL)

% Detecteur:
t_detect = [-1.0 1000];     % temps de detection (un temps negatif implique pas de detection)
xda = 50;%130.0;         % limite de droite du senseur du detecteur
xdb = 150;%230.0;         % limite de gauche du senseur du detecteur

% Numerique :
Nsteps = 5000;
dt = tfin./Nsteps; % pas de temps
Ninters = 256; % nombre d'intervalles du maillage

k0      = 2*pi*n/(xR-xL);
Energie = hbar^2*k0^2/(2*mass);

% Simulations
mat_x    = zeros(length(t_detect),Nsteps+1);
mat_p    = zeros(length(t_detect),Nsteps+1);
mat_Dx   = zeros(length(t_detect),Nsteps+1);
mat_Dp   = zeros(length(t_detect),Nsteps+1);
mat_E    = zeros(length(t_detect),Nsteps+1);
mat_DxDp = zeros(length(t_detect),Nsteps+1);
mat_Pg   = zeros(length(t_detect),Nsteps+1);
mat_Pd   = zeros(length(t_detect),Nsteps+1);
for i=1:length(t_detect)
    filename2  = ['Delta_',num2str(Delta)];
    dt_loc = dt; 
    Delta_loc = Delta;
    x0_loc = x0;
    fname2 = filename2;
    t_d_loc = t_detect(i);
    writeConfig;
    disp('Exercice8_Kervyn_LeMeur configuration.in');   
    system('Exercice8_Kervyn_LeMeur configuration.in'); 
    disp('Done.')
    
% chargement des résultats
    data = load([filename2,'_observables.out']);
    t = data(:,1);
    Pgauche = data(:,2);
    Pdroite = data(:,3);
    E = data(:,5);
    xmoy = data(:,6);
    pmoy = data(:,8);
    Deltax = data(:,11);
    Deltap = data(:,12);
    DeltaxDeltap = data(:,10); %(Delta x)(Delta p)
    data1 = load([filename2,'_potential.out']);
    x = data1(:,1);
    
% affectation aux matrices pour plot superposé
    mat_x(i,:)    = xmoy';
    mat_p(i,:)    = pmoy';
    mat_Dx(i,:)   = Deltax';
    mat_Dp(i,:)   = Deltap';
    mat_E(i,:)    = E';
    mat_DxDp(i,:) = DeltaxDeltap;
    mat_Pg(i,:)   = Pgauche;
    mat_Pd(i,:)   = Pdroite;
end
% plot de f(x)
ViewFormat;
figure('Name','f(x)')
    f = zeros(1,length(x));
    for i=1:length(f)
        if((x(i)>=0)&&(x(i)<xda))
           f(i) = sin(pi*x(i)/(2*xda))^2; 
        end
        if((x(i)>=xda)&&(x(i)<xdb))
           f(i) = 1;
        end
        if((x(i)>=xdb)&&(x(i)<xR))
           f(i) = cos(pi*(x(i)-xdb)/(2*(xR-xdb)))^2;
        end
    end
    plot(x,f,'-','Linewidth',lw);
    set(gca,'fontsize',fs); grid on;
    xlabel('$x$'); ylabel('$f(x)$')
SaveIMG("Detection");

legendStrings = ["No detection","Det. at $t=$ "+num2str(t_detect(2))];

%position moy
figure('Name',"position moyenne")
for i=1:length(t_detect)
    plot(t,mat_x(i,:)','-','Linewidth',lw);
    hold on
end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle x \rangle$')
    legend(legendStrings,'Location','best');
SaveIMG("position moyenne");

%qdm moy
figure('Name',"qdm moyenne")
    for i=1:length(t_detect)
        plot(t,mat_p(i,:)','-','Linewidth',lw);
        hold on
    end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle p \rangle$')
    legend(legendStrings,'Location','best');
SaveIMG("qdm moyenne");

%energie moy
figure('Name',"E moyenne")
    for i=1:length(t_detect)
        plot(t,mat_E(i,:)','-','Linewidth',lw);
        hold on
    end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle E \rangle$')
    legend(legendStrings,'Location','east');
SaveIMG("E moyenne");

%incertitude sur la position
figure('Name',"Deltax moyen")
    for i=1:length(t_detect)
        plot(t,mat_Dx(i,:)','-','Linewidth',lw);
        hold on
    end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle \Delta x \rangle$')
    legend(legendStrings,'Location','northeast');
SaveIMG("Delta x moyen");

%incertitude sur la qdm
figure('Name',"Deltap moyen")
    for i=1:length(t_detect)
        plot(t,mat_Dp(i,:)','-','Linewidth',lw);
        hold on
    end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\langle \Delta p \rangle$')
    legend(legendStrings,'Location','best');
SaveIMG("Delta p moyen");

%Heisenberg
figure('Name',"DeltaxDeltap")
    for i=1:length(t_detect)
        plot(t,mat_DxDp(i,:)','-','Linewidth',lw);
        hold on
    end
    yline(0.5,'r--',{'$\hbar / 2$'},'Interpreter','Latex','LabelOrientation','horizontal','Fontsize',fs,'LineWidth',lw,'LabelVerticalAlignment','top');
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$\Delta x \Delta p$')
    legend(legendStrings,'Location','best');
SaveIMG("DeltaxDeltap");

%Probabilité gauche
figure('Name',"Pgauche")
    for i=1:length(t_detect)
        plot(t,mat_Pg(i,:)','-','Linewidth',lw);
        hold on
    end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$P_{x<x_a}$')
    legend(legendStrings,'Location','southwest');
SaveIMG("Pgauche");

%Probabilité droite
figure('Name',"Pdroite")
    for i=1:length(t_detect)
        plot(t,mat_Pd(i,:)','-','Linewidth',lw);
        hold on
    end
    set(gca,'fontsize',fs); grid on;
    xlabel('$t$'); ylabel('$P_{x>x_b}$')
    legend(legendStrings,'Location','northeast');
SaveIMG("Pdroite");