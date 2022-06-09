% 26 - 01 - 2022
close all;
clear all;
clc;

%% Flags
saveFigs = true;

% Log session
logNum = 89;
% logNum = 58; % no cw
logNum2 = 57; % with cw

% utility data
fontsize = 25;

% Dataset path
datasetPath = '../dvrkDynModelLib/x64/Release/';
fullPath = strcat(datasetPath,'LogSession_',int2str(logNum),'/');
fullPath2 = strcat(datasetPath,'LogSession_',int2str(logNum2),'/');

% Load data
tauMeas = load(strcat(fullPath,'tauMsr.txt'));
tauMod = load(strcat(fullPath,'tauModel.txt'));
tauCmd = load(strcat(fullPath,'tauCmd.txt'));
pdes = load(strcat(fullPath,'pdes.txt'));
Rdes = load(strcat(fullPath,'Rdes.txt'));
pee = load(strcat(fullPath,'pee.txt'));
Ree = load(strcat(fullPath,'Ree.txt'));
g = load(strcat(fullPath,'g.txt'));

% For comparison
pdes2 = load(strcat(fullPath2,'pdes.txt'));
pee2 = load(strcat(fullPath2,'pee.txt'));
g2 = load(strcat(fullPath2,'g.txt'));
tauMeas2 = load(strcat(fullPath2,'tauMsr.txt'));
tauMod2 = load(strcat(fullPath2,'tauModel.txt'));

logTitleFile = fopen(strcat(fullPath,'LogSessionInfo.txt'),'r');
logTitle = textscan(4,'%s','delimiter','\n');
logTitle = logTitle{1};
fclose(logTitleFile);

% Compute Cartesian position error
posErr = zeros(length(pdes),4);
posErr(:,1) = pdes(:,1);
posErr(:,2:4) = pdes(:,2:4) - pee(:,2:4);

posErr2 = zeros(length(pdes2),4);
posErr2(:,1) = pdes2(:,1);
posErr2(:,2:4) = pdes2(:,2:4) - pee2(:,2:4);

% Compute Cartesian orientation error
oriErr = zeros(length(Rdes),4);
abg_ee = zeros(length(Rdes),4);
abg_des = zeros(length(Rdes),4);
oriErr(:,1) = Rdes(:,1);
abg_ee(:,1) = Rdes(:,1);
abg_des(:,1) = Rdes(:,1);
for i = 1 : length(Rdes)

    Ri_des = reshape(Rdes(i,2:10),3,3)';
    Ri = reshape(Ree(i,2:10),3,3)';
    abg_des(i,2:4) = rotm2eul(Ri_des,'XYZ');
    abg_ee(i,2:4) = rotm2eul(Ri,'XYZ');
    oriErr(i,2:4) = abg_des(i,2:4) - abg_ee(i,2:4);
end

dynamicFlag = false;
if(contains(logTitle{1},'dynamic'))
    dynamicFlag = true;
end

trajTypeStr = '';
if(contains(logTitle{1},'rectilinear'))
    ymax = [2 2 5];
    trajTypeStr = 'rectTraj';
elseif(contains(logTitle{1},'spiral'))
    ymax = [5 5 5];
    trajTypeStr = 'spiralTraj';
elseif(contains(logTitle{1},'regulation'))
    ymax = [2 2 2];
    trajTypeStr = 'regulation';
else
    ymax = [2 2 2];
end


%% Plots

if ~dynamicFlag
    % % % 1. Torque comparison - no cw
    trqCompFig = figure;
    strMeasNames = {'$\tau_{1,sim}$';'$\tau_{2,sim}$';'$f_{3,sim}$'};
    strModNames = {'$\tau_{1,mod}$';'$\tau_{2,mod}$';'$f_{3,mod}$'};
    % strModNames = {'$\hat{\tau}_{1}$';'$\hat{\tau}_{2}$';'$\hat{f}_{3}$'};
    yStrNames = {'Torque [Nm]';'Torque [Nm]';'Force [N]'};
    tstart = tauMeas(1,1);
    tend = tauMeas(end,1);

    for i = 1 : 3
        subplot(3,1,i);
    %     if i == 1, title(logTitle,'interpreter','latex');hold on; end
        legStr = {};
        grid on; hold on;
        box on; hold on;
        plot(tauMeas(:,1),tauMeas(:,i+1),'LineWidth',2,'Color',[1 0 0]);hold on;legStr = [legStr;strMeasNames(i);];
        plot(tauMod(:,1),-tauMod(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;strModNames(i);];
        axis([tstart tend -ymax(i) ymax(i)]);
        xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
        ylabel(yStrNames(i),'interpreter','latex','FontSize',fontsize);
        set(gca,'FontSize',fontsize);
        legend(legStr,'interpreter','latex','FontSize',fontsize*2);
    end
    trqCompFig.WindowState = 'fullscreen';
    trqCompFig.PaperPositionMode='auto';
    trqCompFig.PaperOrientation='landscape';
    if saveFigs
        sigNameString = strcat(trajTypeStr,'-nocw-trq-comparison');
        filename_eps = strcat('./',sigNameString,'.eps');
        filename_pdf = strcat('./',sigNameString,'.pdf');
        exportgraphics(trqCompFig,filename_eps,'Resolution',300);
        exportgraphics(trqCompFig,filename_pdf,'Resolution',300);
    end

    % % % 1b. Torque comparison - with cw
    trqCompFig2 = figure;
    strMeasNames = {'$\tau_{1,sim}$';'$\tau_{2,sim}$';'$f_{3,sim}$'};
    strModNames = {'$\tau_{1,mod}$';'$\tau_{2,mod}$';'$f_{3,mod}$'};
    % strModNames = {'$\hat{\tau}_{1}$';'$\hat{\tau}_{2}$';'$\hat{f}_{3}$'};
    yStrNames = {'Torque [Nm]';'Torque [Nm]';'Force [N]'};
    tstart = tauMeas(1,1);
    tend = tauMeas(end,1);
    % ymax = [2 2 5];
    ymax = [5 5 5];
    for i = 1 : 3
        subplot(3,1,i);
    %     if i == 1, title(logTitle,'interpreter','latex');hold on; end
        legStr = {};
        grid on; hold on;
        box on; hold on;
        plot(tauMeas2(:,1),tauMeas2(:,i+1),'LineWidth',2,'Color',[1 0 0]);hold on;legStr = [legStr;strMeasNames(i);];
        plot(tauMod2(:,1),-tauMod2(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;strModNames(i);];
        axis([tstart tend -ymax(i) ymax(i)]);
        xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
        ylabel(yStrNames(i),'interpreter','latex','FontSize',fontsize);
        set(gca,'FontSize',fontsize);
        legend(legStr,'interpreter','latex','FontSize',fontsize*2);
    end
    trqCompFig2.WindowState = 'fullscreen';
    trqCompFig2.PaperPositionMode='auto';
    trqCompFig2.PaperOrientation='landscape';
    if saveFigs
        sigNameString = strcat(trajTypeStr,'-withcw-trq-comparison');
        filename_eps = strcat('./',sigNameString,'.eps');
        filename_pdf = strcat('./',sigNameString,'.pdf');
        exportgraphics(trqCompFig2,filename_eps,'Resolution',300);
        exportgraphics(trqCompFig2,filename_pdf,'Resolution',300);
    end
end


% 2. Cartesian position comparison (des vs actual)
CartPositionFig = figure;
pdes_str = {'$x_d$';'$y_d$';'$z_d$'};
p_str = {'$x$';'$y$';'$z$'};
perr_str = {'$e_x$';'$e_y$';'$e_z$'};
tstart = posErr(1,1);
tend = min(posErr(end,1),posErr2(end,1));
for i = 1 : 3
    subplot(3,1,i);
    grid on;hold on;
    box on; hold on;
    plot(pee(:,1),pee(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;
    plot(pdes(:,1),pdes(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;
    axis([0 tend min(pee(:,i+1))-0.02,pdes(1,i+1)+0.06]);
%     plot(posErr(:,1),posErr(:,i+1),'LineWidth',1,'Color',[0 0 0],'LineStyle','--');hold on;
    legend({p_str{i};pdes_str{i}},'interpreter','latex','FontSize',fontsize*2);
    xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
    ylabel('Position [m]','interpreter','latex','FontSize',fontsize);
    set(gca,'FontSize',fontsize);
end
CartPositionFig.WindowState = 'fullscreen';
CartPositionFig.PaperPositionMode='auto';
CartPositionFig.PaperOrientation='landscape';
if saveFigs
    sigNameString = strcat(trajTypeStr,'-cart-position');
    if dynamicFlag
        sigNameString = strcat('dynamic-',sigNameString);
    end
    filename_eps = strcat('./',sigNameString,'.eps');
    filename_pdf = strcat('./',sigNameString,'.pdf');
    exportgraphics(CartPositionFig,filename_eps,'Resolution',300);
    exportgraphics(CartPositionFig,filename_pdf,'Resolution',300);
end

% Cartesian position error
CartPositionErrFig = figure;
perr_str = {'$e_x$';'$e_y$';'$e_z$'};
scale = 0.8;
err_color = [scale, 0, 0;0, scale, 0;0, 0, scale];
tstart = posErr(1,1);
tend = min(posErr(end,1),posErr2(end,1));
grid on;hold on;
box on; hold on;
for i = 1 : 3
%     plot(posErr(:,1),abs(posErr(:,i+1)),'LineWidth',2,'Color',err_color(i,:));hold on;
    plot(posErr(:,1),posErr(:,i+1),'LineWidth',2,'Color',err_color(i,:));hold on;
end
% axis([0 tend -0.00 0.015])
axis([0 tend min([min(posErr(:,2)),min(posErr(:,3)),min(posErr(:,4))]) max([max(posErr(:,2)),max(posErr(:,3)),max(posErr(:,4))])])
xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
ylabel('Position error [m]','interpreter','latex','FontSize',fontsize);
set(gca,'FontSize',fontsize);
legend(perr_str,'interpreter','latex','FontSize',fontsize*2);
CartPositionErrFig.WindowState = 'fullscreen';
CartPositionErrFig.PaperPositionMode='auto';
CartPositionErrFig.PaperOrientation='landscape';
if saveFigs
    sigNameString = strcat(trajTypeStr,'-cart-position-error');
    if dynamicFlag
        sigNameString = strcat('dynamic-',sigNameString);
    end
    filename_eps = strcat('./',sigNameString,'.eps');
    filename_pdf = strcat('./',sigNameString,'.pdf');
    exportgraphics(CartPositionErrFig,filename_eps,'Resolution',300);
    exportgraphics(CartPositionErrFig,filename_pdf,'Resolution',300);
end

if ~dynamicFlag
    % 2. Cartesian orientation error
    CartOrientationFig = figure;
    pdes_str = {'$\alpha_d$';'$\beta_d$';'$\gamma_d$'};
    p_str = {'$\alpha$';'$\beta$';'$\gamma$'};
    perr_str = {'$e_x$';'$e_y$';'$e_z$'};
    tstart = posErr(1,1);
    tend = min(posErr(end,1),posErr2(end,1));
    ymax = [0.25 0.25 0.25];
    for i = 1 : 3
        subplot(3,1,i);
        grid on;hold on;
        box on; hold on;
        plot(abg_ee(:,1),abg_ee(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;
        plot(abg_des(:,1),abg_des(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;
    %     plot(posErr(:,1),posErr(:,i+1),'LineWidth',1,'Color',[0 0 0],'LineStyle','--');hold on;
        axis([tstart tend -ymax(i) ymax(i)]);
        legend({p_str{i};pdes_str{i}},'interpreter','latex','FontSize',fontsize*2);
        xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
        ylabel('Orientation [rad]','interpreter','latex','FontSize',fontsize);
        set(gca,'FontSize',fontsize);
    end
    CartOrientationFig.WindowState = 'fullscreen';
    CartOrientationFig.PaperPositionMode='auto';
    CartOrientationFig.PaperOrientation='landscape';
    if saveFigs
        sigNameString = strcat(trajTypeStr,'-cart-orientation');
        filename_eps = strcat('./',sigNameString,'.eps');
        filename_pdf = strcat('./',sigNameString,'.pdf');
        exportgraphics(CartOrientationFig,filename_eps,'Resolution',300);
        exportgraphics(CartOrientationFig,filename_pdf,'Resolution',300);
    end
end


if dynamicFlag
    % 3. Commanded torques
    tauCmdFig = figure;
    grid on; hold on;
    box on; hold on;
    tau_color = [scale, 0, 0;0, scale, 0;0, 0, scale];
    strTauNames = {'$u_{1}$';'$u_{2}$';'$u_{3}$'};
    yStrNames = {'Torque [Nm]';'Torque [Nm]';'Force [N]'};
    tstart = tauCmd(1,1);
    tend = tauCmd(end,1);
%     ymax = [50 50 50];
    legStr = {};
    for i = 1 : 3
        subplot(3,1,i);
        grid on;hold on;
        box on;hold on;
        legStr = {};
%         plot(tauCmd(:,1),tauCmd(:,i+1),'LineWidth',2,'Color',tau_color(i,:));hold on;legStr = [legStr;strTauNames(i);];
        plot(tauMod(:,1),tauMod(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;strcat('$\tau_{mod,',int2str(i),'}$');];
        plot(tauMeas(:,1),-tauMeas(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;legStr = [legStr;strcat('$\tau_{sim,',int2str(i),'}$');];
        plot(tauCmd(:,1),tauCmd(:,i+1),'LineWidth',2,'Color',[0 1 0],'LineStyle','--');hold on;legStr = [legStr;strcat('$\tau_{cmd,',int2str(i),'}$');];
        xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
        ylabel(yStrNames(i),'interpreter','latex','FontSize',fontsize);
        set(gca,'FontSize',fontsize);
        legend(legStr,'interpreter','latex','FontSize',fontsize);
    end
%     axis([tstart tend -ymax(i) ymax(i)]);
%     xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
%     ylabel(yStrNames(i),'interpreter','latex','FontSize',fontsize);
%     set(gca,'FontSize',fontsize);
%     legend(legStr,'interpreter','latex','FontSize',fontsize*2);
    tauCmdFig.WindowState = 'fullscreen';
    tauCmdFig.PaperPositionMode='auto';
    tauCmdFig.PaperOrientation='landscape';
    if saveFigs
        sigNameString = strcat(trajTypeStr,'-tauCmd');
        if dynamicFlag
            sigNameString = strcat('dynamic-',sigNameString);
        end
        filename_eps = strcat('./',sigNameString,'.eps');
        filename_pdf = strcat('./',sigNameString,'.pdf');
        exportgraphics(tauCmdFig,filename_eps,'Resolution',300);
        exportgraphics(tauCmdFig,filename_pdf,'Resolution',300);
    end
end

% 2. Cartesian position error for gravity comparison
% figure;
% posErrNames = {'$x \mbox{(no grav.)}$';'$y \mbox{(no grav.)}$';'$z \mbox{(no grav.)}$'};
% posErrNames2 = {'$x \mbox{(with grav.)}$';'$y \mbox{(with grav.)}$';'$z \mbox{(with grav.)}$'};
% tstart = posErr(1,1);
% tend = min(posErr(end,1),posErr2(end,1));
% ymax = [0.05 0.1 0.05];
% logTitle = 'Cartesian position error';
% for i = 1 : 3
%     subplot(3,1,i);
%     if i == 1, title(logTitle,'interpreter','latex');hold on; end
%     legStr = {};
%     grid on; hold on;
%     box on; hold on;
%     plot(posErr(:,1),abs(posErr(:,i+1)),'LineWidth',2,'Color',[1 0 0]);hold on;legStr = [legStr;posErrNames(i);];
%     plot(posErr2(:,1),abs(posErr2(:,i+1)),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;posErrNames2(i);];
%     axis([tstart tend -ymax(i) ymax(i)]);
%     xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
%     ylabel('Position [m]','interpreter','latex','FontSize',fontsize);
%     set(gca,'FontSize',fontsize);
%     legend(legStr,'interpreter','latex','FontSize',fontsize);
% end

% Gravity vector
% figure;
% gNames = {'$g_1 \mbox{(no cw.)}$';'$g_2 \mbox{(no cw.)}$';'$g_3 \mbox{(no cw.)}$'};
% gNames2 = {'$g_1 \mbox{(with cw.)}$';'$g_2 \mbox{(with cw.)}$';'$g_3 \mbox{(with cw.)}$'};
% yGStr = {'Torque [Nm]','Torque [Nm]','Force[N]'};
% tstart = g(1,1);
% tend = min(g(end,1),g2(end,1));
% logTitle = 'Cartesian position error';
% for i = 1 : 3
%     subplot(3,1,i);
%     if i == 1, title(logTitle,'interpreter','latex');hold on; end
%     legStr = {};
%     grid on; hold on;
%     box on; hold on;
%     plot(g(:,1),g(:,i+1),'LineWidth',2,'Color',[1 0 0]);hold on;legStr = [legStr;gNames(i);];
%     plot(g2(:,1),g2(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;gNames2(i);];
%     gmin = min(min(g(:,i+1)),min(g2(:,i+1)));
%     gmax = max(max(g(:,i+1)),max(g2(:,i+1)));
%     axis([tstart tend gmin gmax]);
%     xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
%     ylabel(yGStr{i},'interpreter','latex','FontSize',fontsize);
%     set(gca,'FontSize',fontsize);
%     legend(legStr,'interpreter','latex','FontSize',fontsize);
% end
% 

