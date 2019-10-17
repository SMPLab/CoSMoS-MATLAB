%%
%%
% Developed by Salma Hobbi, Chandra Rupa Rajulapati, and Simon Michael Papalexiou
% 
% Global Institue for Water Security (GIWS) and Global Water Future (GWF)
% University of Saskatchewan 
% 
% Package:		‘CoSMoS-MATLAB’
% Title:		Complete Stochastic Modelling Solution (CoSMoS)
% Version:		v0.9 (beta)
% License:		GPL-3
% Depends:		MATLAB (tested on MATLAB R2018)
% Coded by:		Salma Hobbi
% Conceptual design by:	Salma Hobbi, Chandra Rupa Rajulapati, Simon Michael Papalexiou 
% Contact:		Salma Hobbi (salma.hobbi@usask.ca)
% Repository:		SMPLab 
% Date/Publication:	2019-10-17 16:30:00 CST
% Funding: 		The package was partly funded by the Global institute for Water Security (GIWS; https://www.usask.ca/water/) and the Global Water Futures (GWF; https://gwf.usask.ca/) program.
% 
% Description:  A single framework, unifying, extending, and improving a general-purpose modelling strategy, based on the assumption that any process can emerge by transforming a specific 'parent' Gaussian process (Papalexiou, 2018).
% 
% Redistribution and use the source code, with or –without modification, are permitted if the following conditions are met:
% •	Redistribution of the source code must retain the developers' and affiliated institution names. 
% •	Redistribution of source code must retain all the information in this section along with the following disclaimer
% This software is designed for educational and research purposes only, and the commercial use of this package is strictly prohibited. It is provided by developers 'as is' and any express or implied warranties of merchantability for a particular purpose are disclaimed. Authors, developers and the affiliate institutions accept neither responsibility nor liability for direct/indirect errors, mistakes, and consequential damages; loss of use, data, or profits; or business interruption. User is responsible for any such errors or damages as the package works based on his/her discretion.
% Neither the University of Saskatchewan nor the authors and contributors may be used to endorse or promote products derived from this software without specific prior written permission.

%%
cd([pwd,'\Results'])
load SimulationInfo.mat
load  GeneratedTS.mat
load Distribution.mat
load SimulationInfo.mat
load ACTF.mat
load ACS.mat
Inch_SS = get(0,'screensize');
Inch_SS( Inch_SS == 0 ) = [];
Inch_SS = min( Inch_SS );

figure('Name','GUI plots')
%% Timeseries plot
%% set No of TS
NoTs=1;   %%%%%%% Specify the time series number to visualize
subplot(2,2,1)
plot(GeneratedTS(:,NoTs),'k-')
hold on
ylabel('Synthetic TS','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])
xlabel('Time','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])

set(gca, 'ycolor',[208/255,208/255,208/255]);
        set(gca, 'xcolor',[208/255,208/255,208/255]);
set(gca,'color',[208/255,208/255,208/255]);
      %% exceendance probability plot
subplot(2,2,2)
semilogy(Distribution{1,NoTs}(:,3),Distribution{1,NoTs}(:,1),'k.','MarkerSize',10)
hold on
semilogy(Distribution{1,NoTs}(:,3),Distribution{1,NoTs}(:,2),'k-')
legend('Empirical',['F_{' num2str(SimulationInfo(1).Val{2,2}) '}(X)'],'Location','best','color','none','FontSize',8,'FontWeight','bold')
ylabel('Exceedance Pr','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])
xlabel('Variable','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])
     set(gca,'color',[208/255,208/255,208/255]);
set(gca, 'ycolor',[208/255,208/255,208/255]);
        set(gca, 'xcolor',[208/255,208/255,208/255])

     %% synthetic acs plot  
subplot(2,2,3)
plot(SimulationInfo(6).Val(:,1),ACS{1,NoTs}(:,2),'b-','LineWidth',2)
    hold on
   
plot(SimulationInfo(6).Val(:,1),ACS{1,NoTs}(:,3),'k--')
 hold on
plot(SimulationInfo(6).Val(:,1),ACS{1,NoTs}(:,4),'k.','MarkerSize',10)
xlabel('Lag (\tau)','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255]) ;
ylabel('Auto correlation','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255]) ;
legend({'Target \rho_{x}(\tau)','Gaussian \rho_{z}(\tau)','Empirical'},'color','none','Location','best','FontSize',9)


set(gca,'fontweight','bold','fontsize',8);
     set(gca,'color',[208/255,208/255,208/255]);
set(gca, 'ycolor',[208/255,208/255,208/255]);
        set(gca, 'xcolor',[208/255,208/255,208/255])

        %% autocorrelation values plot

subplot(2,2,4)
plot(ACTF(2,:),ACTF(1,:),'k.-','MarkerSize',12)
hold on
plot(ACTF(1,:),ACTF(1,:),'k--','MarkerSize',12)

xlim([0 max(ACTF(2,:))]);
    set(gca, 'XTick',0:0.1:1,'fontsize',8,'color',[208/255,208/255,208/255])
    set(gca, 'ycolor',[208/255,208/255,208/255])
        set(gca, 'xcolor',[208/255,208/255,208/255])
 xlabel('Target autocorrelation \rho_{x}','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255]) ;
 ylabel('Gaussian \rho_{z}','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255]) ;
        set(gca,'fontweight','bold','fontsize',8)
set(gca,'color',[208/255,208/255,208/255]);
set(gcf,'color',[0.35,0.35,0.35]);

