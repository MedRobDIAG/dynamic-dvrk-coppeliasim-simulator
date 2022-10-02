% 26 - 01 - 2022
close all;
clear all;
clc;

%% Flags
saveFigs = false;
dynamicFlag = false; % keep it false, it will be changed automatically later
jointExcitationFlag = false; % keep it false, it will be changed automatically later

% Log legend
% 1 : kinematic control - 10ms rectilinear trajectory - no cw
% 2 : kinematic control - 10ms rectilinear trajectory - w cw 
% 3 : kinematic control - 10ms spiral trajectory - no cw
% 4 : kinematic control - 10ms spiral trajectory - w cw 
% 5 : kinematic control - 5ms joint exciting trajectories - w cw
% 6 : kinematic control - 5ms joint exciting trajectories - no cw

% 7 : dynamic control - regulation - w cw
% 8 : dynamic control - rectilinear trajectory - w cw
% 9 : dynamic  control - spiral trajectory - w cw

% Log session
logNum = 1; % no cw
logNum_cw = 2; % with cw

% utility data
fontsize = 25;

% Dataset path
datasetPath = '../dvrkDynModelLib/x64/Release/';
fullPath = strcat(datasetPath,'LogSession_',int2str(logNum),'/');
fullPath_cw = strcat(datasetPath,'LogSession_',int2str(logNum_cw),'/');

% Load data
logTitleFile = fopen(strcat(fullPath,'LogSessionInfo.txt'),'r');
logTitle = textscan(4,'%s','delimiter','\n');
logTitle = logTitle{1};
fclose(logTitleFile);
if(contains(logTitle{1},'dynamic'))
    dynamicFlag = true;
end

if(contains(logTitle{1},'excitation'))
    jointExcitationFlag = true;
end



tauMeas = load(strcat(fullPath,'tauMsr.txt'));
tauMod = load(strcat(fullPath,'tauModel.txt'));
pdes = load(strcat(fullPath,'pdes.txt'));
Rdes = load(strcat(fullPath,'Rdes.txt'));
pee = load(strcat(fullPath,'pee.txt'));
Ree = load(strcat(fullPath,'Ree.txt'));
g = load(strcat(fullPath,'g.txt'));


% For comparison
pdes_cw = load(strcat(fullPath_cw,'pdes.txt'));
pee_cw = load(strcat(fullPath_cw,'pee.txt'));
g_cw = load(strcat(fullPath_cw,'g.txt'));
tauMeas_cw = load(strcat(fullPath_cw,'tauMsr.txt'));
tauMod_cw = load(strcat(fullPath_cw,'tauModel.txt'));

if dynamicFlag || (~dynamicFlag && jointExcitationFlag)
    tauCmd_ncw = load(strcat(fullPath,'tauCmd.txt'));
    tauCmd_cw = load(strcat(fullPath_cw,'tauCmd.txt'));
end

% Compute Cartesian position error
posErr_ncw = zeros(length(pdes),4);
posErr_ncw(:,1) = pdes(:,1);
posErr_ncw(:,2:4) = pdes(:,2:4) - pee(:,2:4);

posErr_cw = zeros(length(pdes_cw),4);
posErr_cw(:,1) = pdes_cw(:,1);
posErr_cw(:,2:4) = pdes_cw(:,2:4) - pee_cw(:,2:4);

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
elseif(contains(logTitle{1},'excitation'))
    trajTypeStr = 'joint-excitationTraj';
end

%% RMSE tracking accuracy
start_t = 0.5;
start_idx = find(abs(posErr_ncw(:,1)-start_t)<1e-2);start_idx = start_idx(1);
rmsTrackPos_ncw = rms(posErr_ncw(:,2:4));
rmsTrackPos_tot_ncw = sqrt(rmsTrackPos_ncw(1)^2+rmsTrackPos_ncw(2)^2+rmsTrackPos_ncw(3)^2);

start_idx = find(abs(posErr_cw(:,1)-start_t)<1e-2);start_idx = start_idx(1);
rmsTrackPos_cw = rms(posErr_cw(:,2:4));
rmsTrackPos_tot_cw = sqrt(rmsTrackPos_cw(1)^2+rmsTrackPos_cw(2)^2+rmsTrackPos_cw(3)^2);

start_idx = find(abs(oriErr(:,1)-start_t)<1e-2);start_idx = start_idx(1);
rmsTrackOri_cw = rms(oriErr(:,2:4));
rmsTrackOri_tot_cw = sqrt(rmsTrackOri_cw(1)^2+rmsTrackOri_cw(2)^2+rmsTrackOri_cw(3)^2);

%% RMSE Torques
start_t = 0.5;
start_idx = find(abs(tauMeas(:,1)-start_t)<1e-2);start_idx = start_idx(1);
if ~jointExcitationFlag
    tauErr_ncw = tauMeas(start_idx:end,2:4)-(-tauMod(start_idx:end,2:4));
else
    tauErr_ncw = tauMeas(start_idx:end,2:4)-(-tauCmd_ncw(start_idx:end,2:4));
end

relRMS_ncw = norm(tauErr_ncw)/norm(tauMod(start_idx:end,2:4));
RMS_ncw = sqrt(sum(tauErr_ncw.^2)/size(tauErr_ncw,1))
RMS_12_ncw = sqrt(RMS_ncw(1)^2+RMS_ncw(2)^2);
RMS_3_ncw = RMS_ncw(3);
RMS_tot_ncw = RMS_12_ncw + RMS_3_ncw;

start_idx = find(abs(tauMeas_cw(:,1)-start_t)<1e-2);start_idx = start_idx(1);
if ~jointExcitationFlag 
    tauErr_cw = tauMeas_cw(start_idx:end,2:4)-(-tauMod_cw(start_idx:end,2:4));
else
    tauErr_cw = tauMeas_cw(start_idx:end,2:4)-(-tauCmd_cw(start_idx:end,2:4));
end

relRMS_cw = norm(tauErr_cw)/norm(tauMod_cw(start_idx:end,2:4));
RMS_cw = sqrt(sum(tauErr_cw.^2)/size(tauErr_cw,1))
RMS_12_cw = sqrt(RMS_cw(1)^2+RMS_cw(2)^2);
RMS_3_cw = RMS_cw(3);
RMS_tot_cw = RMS_12_cw + RMS_3_cw;


%% Plots

if ~dynamicFlag
    % % % 1. Torque comparison - no cw
    trqCompFig = figure;
    strMeasNames = {'$\tau_{1,sim}$';'$\tau_{2,sim}$';'$f_{3,sim}$'};
    strModNames = {'$\tau_{1,mod}$';'$\tau_{2,mod}$';'$f_{3,mod}$'};
    strWangNames = {'$\tau_{1,W}$';'$\tau_{2,W}$';'$f_{3,W}$'};
    yStrNames = {'Torque [Nm]';'Torque [Nm]';'Force [N]'};
    tstart = tauMeas(1,1);
    tend = tauMeas(end,1);   
    ymax = [5 5 5];

    for i = 1 : 3
        subplot(3,1,i);
        legStr = {};
        grid on; hold on;
        box on; hold on;
        plot(tauMod(:,1),-tauMod(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;strModNames(i);];
        plot(tauMeas(:,1),tauMeas(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;legStr = [legStr;strMeasNames(i);];
        if(jointExcitationFlag)
            plot(tauCmd_ncw(:,1),-tauCmd_cw(:,i+1),'LineWidth',2,'Color',[0 1 0]);hold on;legStr = [legStr;strWangNames(i);];
        end
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
        filename_eps = strcat(fullPath,sigNameString,'.eps');
        filename_pdf = strcat(fullPath,sigNameString,'.pdf');
        exportgraphics(trqCompFig,filename_eps,'Resolution',300);
        exportgraphics(trqCompFig,filename_pdf,'Resolution',300);
    end
    
    % % % 1b. Torque comparison - with cw
    trqCompFig2 = figure;
    strMeasNames = {'$\tau_{1,sim}$';'$\tau_{2,sim}$';'$f_{3,sim}$'};
    strModNames = {'$\tau_{1,mod}$';'$\tau_{2,mod}$';'$f_{3,mod}$'};
    strWangNames = {'$\tau_{1,W}$';'$\tau_{2,W}$';'$f_{3,W}$'};
    yStrNames = {'Torque [Nm]';'Torque [Nm]';'Force [N]'};
    tstart = tauMeas_cw(1,1);
    tend = tauMeas_cw(end,1);
    ymax = [5 5 5];
    
    
    for i = 1 : 3
        subplot(3,1,i);
        legStr = {};
        grid on; hold on;
        box on; hold on;
        plot(tauMod_cw(:,1),-tauMod_cw(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;strModNames(i);];
        plot(tauMeas_cw(:,1),tauMeas_cw(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;legStr = [legStr;strMeasNames(i);];
        if(jointExcitationFlag)
            plot(tauCmd_cw(:,1),-tauCmd_cw(:,i+1),'LineWidth',2,'Color',[0 1 0]);hold on;legStr = [legStr;strWangNames(i);];
        end
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
        filename_eps = strcat(fullPath_cw,sigNameString,'.eps');
        filename_pdf = strcat(fullPath_cw,sigNameString,'.pdf');
        exportgraphics(trqCompFig2,filename_eps,'Resolution',300);
        exportgraphics(trqCompFig2,filename_pdf,'Resolution',300);
    end
end

if ~jointExcitationFlag
    % 2. Cartesian position comparison (des vs actual)
    CartPositionFig = figure;
    pdes_str = {'$x_d$';'$y_d$';'$z_d$'};
    p_str = {'$x$';'$y$';'$z$'};
    perr_str = {'$e_x$';'$e_y$';'$e_z$'};
    tstart = posErr_ncw(1,1);
    tend = min(posErr_ncw(end,1),posErr_cw(end,1));
    for i = 1 : 3
        subplot(3,1,i);
        grid on;hold on;
        box on; hold on;
        plot(pee(:,1),pee(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;
        plot(pdes(:,1),pdes(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;
        axis([0 tend min(pee(:,i+1))-0.02,pdes(1,i+1)+0.06]);
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
        filename_eps = strcat(fullPath,sigNameString,'.eps');
        filename_pdf = strcat(fullPath,sigNameString,'.pdf');
        exportgraphics(CartPositionFig,filename_eps,'Resolution',300);
        exportgraphics(CartPositionFig,filename_pdf,'Resolution',300);
    end

    % Cartesian position error
    CartPositionErrFig = figure;
    perr_str = {'$e_x$';'$e_y$';'$e_z$'};
    scale = 0.8;
    err_color = [scale, 0, 0;0, scale, 0;0, 0, scale];
    tstart = posErr_ncw(1,1);
    tend = min(posErr_ncw(end,1),posErr_cw(end,1));
    grid on;hold on;
    box on; hold on;
    for i = 1 : 3
    %     plot(posErr(:,1),abs(posErr(:,i+1)),'LineWidth',2,'Color',err_color(i,:));hold on;
        plot(posErr_ncw(:,1),posErr_ncw(:,i+1),'LineWidth',2,'Color',err_color(i,:));hold on;
    end
    % axis([0 tend -0.00 0.015])
    axis([0 tend min([min(posErr_ncw(:,2)),min(posErr_ncw(:,3)),min(posErr_ncw(:,4))]) max([max(posErr_ncw(:,2)),max(posErr_ncw(:,3)),max(posErr_ncw(:,4))])])
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
        filename_eps = strcat(fullPath,sigNameString,'.eps');
        filename_pdf = strcat(fullPath,sigNameString,'.pdf');
        exportgraphics(CartPositionErrFig,filename_eps,'Resolution',300);
        exportgraphics(CartPositionErrFig,filename_pdf,'Resolution',300);
    end

    if ~dynamicFlag
        % 2. Cartesian orientation error
        CartOrientationFig = figure;
        pdes_str = {'$\alpha_d$';'$\beta_d$';'$\gamma_d$'};
        p_str = {'$\alpha$';'$\beta$';'$\gamma$'};
        perr_str = {'$e_x$';'$e_y$';'$e_z$'};
        tstart = posErr_ncw(1,1);
        tend = min(posErr_ncw(end,1),posErr_cw(end,1));
        ymax = [0.25 0.25 0.25];
        for i = 1 : 3
            subplot(3,1,i);
            grid on;hold on;
            box on; hold on;
            plot(abg_ee(:,1),abg_ee(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;
            plot(abg_des(:,1),abg_des(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;
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
            filename_eps = strcat(fullPath,sigNameString,'.eps');
            filename_pdf = strcat(fullPath,sigNameString,'.pdf');
            exportgraphics(CartOrientationFig,filename_eps,'Resolution',300);
            exportgraphics(CartOrientationFig,filename_pdf,'Resolution',300);
        end
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
    legStr = {};
    for i = 1 : 3
        subplot(3,1,i);
        grid on;hold on;
        box on;hold on;
        legStr = {};
        plot(tauMod(:,1),tauMod(:,i+1),'LineWidth',2,'Color',[0 0 1]);hold on;legStr = [legStr;strcat('$\tau_{mod,',int2str(i),'}$');];
        plot(tauMeas(:,1),-tauMeas(:,i+1),'LineWidth',2,'Color',[1 0 0],'LineStyle','--');hold on;legStr = [legStr;strcat('$\tau_{sim,',int2str(i),'}$');];
        xlabel('Time [s]','interpreter','latex','FontSize',fontsize);
        ylabel(yStrNames(i),'interpreter','latex','FontSize',fontsize);
        set(gca,'FontSize',fontsize);
        legend(legStr,'interpreter','latex','FontSize',fontsize);
    end
    tauCmdFig.WindowState = 'fullscreen';
    tauCmdFig.PaperPositionMode='auto';
    tauCmdFig.PaperOrientation='landscape';
    if saveFigs
        sigNameString = strcat(trajTypeStr,'-tauCmd');
        if dynamicFlag
            sigNameString = strcat('dynamic-',sigNameString);
        end
        filename_eps = strcat(fullPath,sigNameString,'.eps');
        filename_pdf = strcat(fullPath,sigNameString,'.pdf');
        exportgraphics(tauCmdFig,filename_eps,'Resolution',300);
        exportgraphics(tauCmdFig,filename_pdf,'Resolution',300);
    end
end

