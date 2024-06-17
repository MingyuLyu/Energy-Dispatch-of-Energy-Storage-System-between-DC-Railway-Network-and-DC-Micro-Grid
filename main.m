function main()
clc; 
close all; 
clear;

%% Constants 
time_step=1/12;    %time step
hours=24;        %Total test hours
T=1/time_step*hours; % 60 minutes per hour
Price_grid=17.54; % electriciy price
Price_ev=28;     % EV charging price
Price_sell=11;   % Energy selling price

einv=0.95;   

Qc=500;        % ESS capacity (kWh)
% Qev=60;         % EV battery capacity (kWh)
einv=0.95;
ech =0.95;      %充電効率
edisch =0.95;   %放電効率
Price_deg=5.06;    % battery degredation cost
% Price_deg=0;    % battery degredation cost
%Pgmax=100;       % Maximum grid power 
% Pevmax=50;         % Maximum EV charging power
InitialSOE = 0.4;  % Initial SOE of ESS
minSOE=0.2;         % Minimum SOE of ESS
maxSOE=0.8;         % Maximun SOE of ESS

%%%%%%%%%%%%%%%% RealTime Electricity Price %%%%%%%%%%%%%%%%%%%%%%%%%%

Price_gridr=[15.96 15.96 15.96 15.57 15.57 15.36 15.55 15.78 15.97 18.50 22 22.9 ...
    19.38 21.29 21.12 21.21 17.42 17.22 21.38 26.71 26.93 26.94 26.93...
    26.94 27.25 27.21 26.94 26.98 26.94 26.94 26.94 26.93 25.76 25.13...
    23.27 19.28 18.96 17.22 15.96 15.57 15.97 15.68 15.68 15.36 15.34 ...
    15.29 15.34 15.57] ;    %real time electricity price
Price_gridreal=interp1(1:6:T,Price_gridr,1:1:T,'spline');                                                      % Power of PV



%%%%%%%%%%%%%%%% EV power %%%%%%%%%%%%%%%%%%%%%%%%%%

load('EVpower.mat','EVpower2');
Pev=EVpower2'/10; % Power of EV charging   % good result
% Pev=EVpower2'*0;

%%%%%%%%%%%%%%%% PV power %%%%%%%%%%%%%%%%%%%%%%%%%%
ppv=csvread('PVpower_1.csv')/3;
% ppv=csvread('PVpower_1.csv')*0;
Ppv=zeros(1,T);

for i=1:48
    Ppv(1,T/48*i)=ppv(i);
end
% ppva(1:40)=ppv(9:48); % time axis start from 0:00 to 4:00
% ppva(1:40)=ppv(9:48)*0; % time axis start from 0:00 to 4:00   % good result
Ppv=interp1(1:6:T,ppv,1:1:T,'spline');                                                      % Power of PV
Ppv(1,1:12)=0;
Ppv(Ppv<0)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   Single Trian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pr2ess_%
% t=[4 299 320 350 368 388 397 406 413 416 423 428 434 438 442 445 448 ...
%     451 453 459 462 465 468 472 475 478 482 485 489 492 496 499 502 504 509 ...
%     513 521 525 532 537 542 552 562 576 590 597 608 616 625 641 655 666 677 685 ...
%     702 715 726 736 745 762 776 786 797 806 819 835 846 853 864 872 883 895 ...
%     906 913 925 936 943 955 966 970 976 985 992 999 1002 1007 1110 1118 1028 ...
%     1036 1043 1051 1056 1064 1069 1079 1084 1093 1096 1108 1115 1119 1128 ...
%     1140 1146 1150 1156 1162 1169 1176 1180 1185 1191 1202 1205 1209 1222 1227 1237 ...
%     1246 1257 1264 1269 1274 1280 1291 1295 1300 1308 1313 1321 1331 1338 1356 ...
%     1367 1377 1389 1397 1405 1421];
% pra=[2 14 25 327 344 365 383 397 418 440 447 456 461 470 480 483 489 491 495 502 ...
%     505 509 514 517 520 524 532 536 539 543 547 551 555 562 567 573 577 582 588 ...
%     593 599 608 613 622 632 639 652 660 670 682 689 700 712 726 742 750 760 772 ...
%     780 790 802 810 820 832 846 882 870 880 912 905 910 922 930 940 952 962 969 ...
%     982 995 1002 1012 1020 1030 1042 1050 1054 1060 1065 1072 1080 1084 1089 1094 1102 ...
%     1110 1115 1120 1133 1140 1143 1149 1156 1163 1170 1175 1192 1198 1205 1214 1216 ...
%     1224 1233 1239 1246 1252 1255 1261 1267 1272 1278 1284 1290 1296 1301 1308 1315 ...
%     1323 1334 1336 1343 1351 1360 1365 1373 1382 1392 1401 1410 1419 1429];
% t=pra-pr2ess_t;
% B=abs(t)>5;
% 
% Pr2ess=1.11*[zeros(1,74),B,zeros(1,74)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%   Multiple Trian   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
Rep2=xlsread('2-squeez.xlsx','B2:EHM63');   %Read regenerative energy informatio (headway 2min)
Rep10=xlsread('10-squeez.xlsx','B2:EHM13'); %Read regenerative energy informatio (headway 10min)
Rep5=xlsread('5-squeeze.xlsx','B2:EHM25');   %Read regenerative energy informatio (headway 5min)

Rep2=sum(Rep2);                         
Rep10=sum(Rep10);
Rep5=sum(Rep5);

Rep2=reshape(Rep2,[3600*time_step,1/time_step]);
Rep10=reshape(Rep10,[3600*time_step,1/time_step]);
Rep5=reshape(Rep5,[3600*time_step,1/time_step]);

Rep2=-1*sum(Rep2)/(1000*60/time_step);                         
Rep10=-1*sum(Rep10)/(1000*60/time_step);
Rep5=-1*sum(Rep5)/(1000*60/time_step);


Rep2_5=xlsread('squeez-2-5.xlsx','B2:EHM43');   %Read regenerative energy informatio (headway 2-5min)
Rep5_2=xlsread('squeez-5-2.xlsx','B2:EHM43'); %Read regenerative energy informatio (headway 5-2min)
Rep5_10=xlsread('squeez-5-10.xlsx','B2:EHM25');   %Read regenerative energy informatio (headway 5-10min)
Rep10_5=xlsread('squeez-10-5.xlsx','B2:EHM21');   %Read regenerative energy informatio (headway 10-5min)

Rep2_5=sum(Rep2_5);                         
Rep5_2=sum(Rep5_2);
Rep5_10=sum(Rep5_10);
Rep10_5=sum(Rep10_5);  

Rep2_5=reshape(Rep2_5,[3600*time_step,1/time_step]);
Rep5_2=reshape(Rep5_2,[3600*time_step,1/time_step]);
Rep5_10=reshape(Rep5_10,[3600*time_step,1/time_step]);
Rep10_5=reshape(Rep10_5,[3600*time_step,1/time_step]);

Rep2_5=-1*sum(Rep2_5)/(1000*60/time_step);                         
Rep5_2=-1*sum(Rep5_2)/(1000*60/time_step);
Rep5_10=-1*sum(Rep5_10)/(1000*60/time_step);
Rep10_5=-1*sum(Rep10_5)/(1000*60/time_step); 

Rep2_=Rep2;
Rep2_(1:1/(2*time_step))=Rep2(1+1/(2*time_step):1/time_step);
Rep10_1=Rep10;
Rep10_1(1:1/(2*time_step))=Rep10(1+1/(2*time_step):1/time_step);
Rep5_=Rep5;
Rep5_(1:1/(2*time_step))=Rep5(1+1/(2*time_step):1/time_step);

Rep2_5_=zeros(1,1/time_step);
Rep2_5_(1:1/(2*time_step))=Rep2_5(1+1/(2*time_step):1/time_step);
Rep2_5_(1+1/(2*time_step):1/time_step)=Rep5(1+1/(2*time_step):1/time_step);

Rep5_2_=zeros(1,1/time_step);
Rep5_2_(1:1/(2*time_step))=Rep5(1+1/(2*time_step):1/time_step);
Rep5_2_(1+1/(2*time_step):1/time_step)=Rep5_2(1+1/(2*time_step):1/time_step);

Rep5_10_=zeros(1,1/time_step);
Rep5_10_(1:1/(2*time_step))=Rep5(1+1/(2*time_step):1/time_step);
Rep5_10_(1+1/(2*time_step):1/time_step)=Rep5_10(1+1/(2*time_step):1/time_step);

Rep10_5_=zeros(1,1/time_step);
Rep10_5_(1:1/(2*time_step))=Rep10(1+1/(2*time_step):1/time_step);
Rep10_5_(1+1/(2*time_step):1/time_step)=Rep10_5(1+1/(2*time_step):1/time_step);


% Pr2ess=[Rep10,Rep10_1,Rep10_1,Rep10_5_,repmat(Rep5_,1,2),Rep5_10_,repmat(Rep10_1,1,5),Rep10_5_,repmat(Rep5_,1,2),Rep5_10,repmat(Rep10_1,1,3),Rep10_1,fliplr(Rep10),zeros(1,3/time_step)]/3;
Pr2ess=[Rep10,Rep10_5_,Rep5_,Rep5_2_,repmat(Rep2_,1,2),Rep2_5_,repmat(Rep5_,1,5),Rep5_2_,repmat(Rep2_,1,2),Rep2_5_,repmat(Rep5_,1,3),Rep5_10_,fliplr(Rep10),zeros(1,3/time_step)]*5;
% Pr2ess=[zeros(1,3/time_step),Rep2_5,zeros(1,3/time_step),Rep5_2,zeros(1,3/time_step),Rep5_10,zeros(1,3/time_step),Rep10_5,zeros(1,8/time_step)];
%Pr2ess=[repmat(Rep10,1,2),repmat(Rep5,1,2),repmat(Rep2,1,2),repmat(Rep5,1,7),repmat(Rep2,1,2),repmat(Rep5,1,5),repmat(Rep10,1,1),zeros(1,36)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nvars=6*T;     % Total variables


%% Constraints 

% Construction of x:
% 1T: PV2Grid
% 2T: PV2EV
% 3T: PV2ESS
% 4T: Grid2EV
% 5T: ESS2EV
% 6T: ESS2Grid
%%%%%%%%%%%%%%%%%%%%%%%%%%% Linear inequality constraints %%%%%%%%%%%%%%%%%%%%%





% energy going in ess 
pv2ess=tril(ones(T-1))*ech*time_step/Qc;  
pv2ess=[zeros(1,T-1);pv2ess];                  % The first time step is depending on the initial SOE
pv2ess=[pv2ess,zeros(T,1)];
% energy going out of ess
ess2ev=-1*tril(ones(T-1))*time_step/Qc;
ess2ev=[zeros(1,T-1);ess2ev];                  % The first time step is depending on the initial SOE
ess2ev=[ess2ev,zeros(T,1)];

ess2g=-1*tril(ones(T-1))*time_step/Qc;
ess2g=[zeros(1,T-1);ess2g];                     % The first time step is depending on the initial SOE
ess2g=[ess2g,zeros(T,1)];


SOEup=[zeros(T),zeros(T),pv2ess,zeros(T),ess2ev,ess2g];
SOElow=-1*SOEup;

% % Pessout=Pess2ev+Pess2grid
% pessout=zeros(T,9*T);
% for i=1:T
%     pessout(i,T+i)=1;
%     pessout(i,2*T+i)=1;
%     pessout(i,6*T+i)=1;
% end
% 
% % Pessin=Ppv2ess+Pgrid2ess+Pr2ess
% pessin=zeros(T,9*T);
% for i=1:T
%     pessin(i,4*T+i)=-ech;
%     pessin(i,7*T+i)=-ech;
% end

A=[SOEup;SOElow];

% regenerate energy
pr2ess=zeros(T,1);
for i=1:T
    pr2ess(i)=sum(Pr2ess(1,1:i));
end
pr2ess=pr2ess*ech*time_step/Qc;


bup=maxSOE*ones(T,1)-InitialSOE-pr2ess;
blow=-minSOE*ones(T,1)+InitialSOE+pr2ess;

% bpessout=1000*ones(T,1);
% bpessin=-750*ones(T,1)+ech*Pr2ess;

b=[bup;blow];

%%%%%%%%%%%%%%%%%%%%% Linear equality constraints %%%%%%%%%%%%%%%%%%%%%%%%

ev=zeros(T,6*T);    % Pev=Pgrid2ev+Pess2ev+Ppv2ev
for i=1:T
    ev(i,T+i)=1;
    ev(i,3*T+i)=edisch;
    ev(i,4*T+i)=einv;
end

pv=zeros(T,6*T);    % Ppv=Ppv2grid+Ppv2ev+Ppv2ess
for i=1:T
    pv(i,i)=1;
    pv(i,T+i)=edisch;
    pv(i,2*T+i)=einv;
end


%%% Final SOC return to 0.5 %%%
finpv2ess=ones(1,T)*ech*time_step/Qc;

finess2ev=-1*ones(1,T)*time_step/Qc;

finess2g=-1*ones(1,T)*time_step/Qc;

finsoe=[zeros(1,T),zeros(1,T),finpv2ess,zeros(1,T),finess2ev,finess2g];


% finpv2ess_1=ones(1,T-1)*ech*time_step/Qc;
% finpv2ess_1=[finpv2ess_1,0];
% finess2ev_1=-1*ones(1,T-1)*time_step/Qc;
% finess2ev_1=[finess2ev_1,0];
% finess2g_1=-1*ones(1,T-1)*time_step/Qc;
% finess2g_1=[finess2g_1,0];
% finsoe_1=[zeros(1,T),zeros(1,T),finpv2ess_1,zeros(1,T),finess2ev_1,finess2g_1];


% finpv2ess_2=ones(1,T-2)*ech*time_step/Qc;
% finpv2ess_2=[finpv2ess_2,zeros(1,5)];
% finess2ev_2=-1*ones(1,T-2)*time_step/Qc;
% finess2ev_2=[finess2ev_2,zeros(1,5)];
% finess2g_2=-1*ones(1,T-2)*time_step/Qc;
% finess2g_2=[finess2g_2,zeros(1,5)];
% finsoe_2=[zeros(1,T),zeros(1,T),finpv2ess_2,zeros(1,T),finess2ev_2,finess2g_2];


% finpv2ess_5=ones(1,T-5)*ech*time_step/Qc;
% finpv2ess_5=[finpv2ess_5,zeros(1,5)];
% finess2ev_5=-1*ones(1,T-5)*time_step/Qc;
% finess2ev_5=[finess2ev_5,zeros(1,5)];
% finess2g_5=-1*ones(1,T-5)*time_step/Qc;
% finess2g_5=[finess2g_5,zeros(1,5)];
% finsoe_5=[zeros(1,T),zeros(1,T),finpv2ess_5,zeros(1,T),finess2ev_5,finess2g_5];


FinSOE=0.4-InitialSOE-sum(Pr2ess)*ech*time_step/Qc;
% 
Aeq=[ev;pv;finsoe];
beq=[Pev;Ppv';FinSOE];

% Aeq=[ev;pv];
% beq=[Pev;Ppv];

% upper and lower bounds
lb=zeros(6*T,1);

ub=[Inf(4*T,1);0.25*Qc*ones(2*T,1)];

%% Function

Cost=zeros(1,6*T);
% 1T: PV2Grid
for i=1:T
    Cost(1,i)=-time_step*einv*Price_sell;
    Cost(1,T+i)=-time_step*einv*Price_ev;
%     Cost(1,3*T+i)=-time_step*(Price_ev-Price_grid);
    Cost(1,3*T+i)=-time_step*(Price_ev-Price_gridreal(1,i));
    Cost(1,4*T+i)=-time_step*(edisch*Price_ev-Price_deg);
%    Cost(1,5*T+i)=-1*time_step*(Price_sell-Price_deg);
Cost(1,5*T+i)=-1*time_step*(-Price_deg);
end

f=Cost;
% Construction of x:
% 1T: PV2Grid
% 2T: PV2EV
% 3T: PV2ESS
% 4T: Grid2EV
% 5T: ESS2EV
% 6T: ESS2Grid
% objfun=@Fun;
% F
% options=gaoptimset('paretoFraction',0.3,'populationsize',200,'generations',50,'stallGenLimit',100,'TolFun',1e-10,'PlotFcns',@gaplotpareto);
%%
% options=optimoptions('paretoFraction',0.3,'populationsize',200,'generations',300,'stallGenLimit',200,'TolFun',1e-10,'PlotFcns',@gaplotpareto);

% [x,fval,exitflag,output,population,scores] =gamultiobj(objfun,nvars,A,b,Aeq,beq,lb,ub,options);
x=intlinprog(f,nvars,A,b,Aeq,beq,lb,ub);

%% Results
% 因为gamultiobj是以目标函数分量取极小值为目标，
% 因此在y=Fun(x)里取相反数的目标函数再取相反数画出原始情况
% figure(1)
% plot(fval(:,1),fval(:,2),'pr')
% xlabel('f_1(x)')
% ylabel('f_2(x)')
% title('Pareto front')
% grid on

%% Plot
% Construction of x:
% 1T: PV2Grid
% 2T: PV2EV
% 3T: PV2ESS
% 4T: Grid2EV
% 5T: ESS2EV

xplot=zeros(6,T);
for i=1:6
    xplot(i,:)=x((i-1)*T+1:i*T,1);
end 

t=1:1:T;



%% Figure
figure(1)
yyaxis left
plot(t,xplot,'linewidth',2);
grid on
colororder('default')
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
set(gca,'FontSize',20);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',25);

yyaxis right
plot(t,Price_gridreal,'--','linewidth',0.5);
colororder('m');
ylabel('Grid Electricity Price/yen','Fontname', 'Times New Roman','FontSize',25);

% title('Power/ kW','Fontname', 'Times New Roman','FontSize',20);

legend('PV2Grid','PV2EV','PV2ESS','Grid2EV','ESS2EV','ESS2Grid','Price');
% figure(3)
% plot(t,x(opt,T+1:2*T));
% grid on
% title('ESS2EV/ kW','Fontname', 'Times New Roman','FontSize',20);
% 
% figure(4)
% plot(t,x(opt,2*T+1:3*T));
% grid on
% title('ESS2Grid/ kW','Fontname', 'Times New Roman','FontSize',20);
% 
% figure(5)
% plot(t,x(opt,3*T+1:4*T));
% grid on
% title('PV2Grid/ kW','Fontname', 'Times New Roman','FontSize',20);
% 
% figure(6)
% plot(t,x(opt,4*T+1:5*T));
% grid on
% title('Grid2ESS/ kW','Fontname', 'Times New Roman','FontSize',20);
% 
% figure(7)
% plot(t,x(opt,5*T+1:6*T));
% grid on
% title('Grid2R/ kW','Fontname', 'Times New Roman','FontSize',20);
% 
% figure(8)
% plot(t,x(opt,6*T+1:7*T));
% grid on
% title('ESS2R/ kW','Fontname', 'Times New Roman','FontSize',20);
% 
% figure(9)
% plot(t,x(opt,7*T+1:8*T));
% grid on
% title('PV2ESS/ kW','Fontname', 'Times New Roman','FontSize',20);

%% SOE

SOE=zeros(1,T);
for i=1:T
    SOE(1,i)=InitialSOE+time_step*(sum(xplot(3,1:i))-sum(xplot(5,1:i))-sum(xplot(6,1:i))+sum(Pr2ess(1:i)))/Qc;
end
figure(2)
plot(t,SOE,'linewidth',2);
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'YTick',0:0.2:1);
set(gca,'YTickLabel',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'XLim',[0 T]);
set(gca,'YLim',[0 1]);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('SOE','Fontname', 'Times New Roman','FontSize',25);
set(gca,'FontSize',20);
grid on
title('ESS SOE','Fontname', 'Times New Roman','FontSize',20);

%% Toatal EV charging power
figure(3)
t=1:1:T;
bar(t,Pev);
title('EV charging power demand profile','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',25);
grid on
% set(AX,'xTick',1:time_step:T,'linewidth',2,'fontsize',20,'fontname','Times');
% set(gca,'XtickLabel',0:1:T);

%% Toatal Regenerative power
figure(4)
t=1:1:T;
bar(t,Pr2ess,0.5);
title('Railway regenerative power profile','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',25);
grid on
%% Toatal PV power
figure(5)
t=1:1:T; 
bar(t,Ppv,0.5);
title('PV power profile','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',25);
grid on
%% Total Figure
figure(6)
t=1:1:T; 
subplot(3,1,1)
bar(t,Pr2ess,0.5);
title('Railway regenerative power profile','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',10);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
set(gca,'FontSize',20);
xlabel('time','Fontname', 'Times New Roman','FontSize',20);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',20);
grid on

subplot(3,1,2)
bar(t,Ppv,0.5);
title('PV power profile','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',10);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
set(gca,'FontSize',20);
xlabel('time','Fontname', 'Times New Roman','FontSize',20);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',20);
grid on


subplot(3,1,3)
bar(t,Pev);
title(' EV charging power demand profile','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',10);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
set(gca,'FontSize',20);
xlabel('time','Fontname', 'Times New Roman','FontSize',20);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',20);
grid on


%% Electricity Price
figure(7)
t=1:1:T; 
% plot(t,Price_gridreal,'linewidth',3);
plot(1:6:T,Price_gridr,'linewidth',3);
hold on
plot(1:6:T,Price_gridr,'.','MarkerSize',35,'MarkerFaceColor',[0 114 189]/225 ,'MarkerEdgeColor',[1 114 189]/225);
title('Real-time electricity price','Fontname', 'Times New Roman','FontSize',20);
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('yen/kWh','Fontname', 'Times New Roman','FontSize',25);
grid on

%% Figure without PV EV
figure(8)
yyaxis left
plot(t,xplot,'linewidth',2);
grid on
colororder('default')
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'XLim',[0 T]);
set(gca,'FontSize',20);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('Power/kW','Fontname', 'Times New Roman','FontSize',25);

yyaxis right 
plot(t,SOE,'--','linewidth',1);
set(gca,'linewidth',2,'fontsize',15);
set(gca,'XTick',1:1/time_step:T);
set(gca,'XTickLabel',[4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 1 2 3 4]);
set(gca,'YTick',0:0.2:1);
set(gca,'YTickLabel',[0 0.2 0.4 0.6 0.8 1]);
set(gca,'XLim',[0 T]);
set(gca,'YLim',[0 1]);
xlabel('time','Fontname', 'Times New Roman','FontSize',25);
ylabel('SOE','Fontname', 'Times New Roman','FontSize',25);
set(gca,'FontSize',20);
grid on
legend('PV2Grid','PV2EV','PV2ESS','Grid2EV','ESS2EV','ESS2Grid','SOE');
% title('Power/ kW','Fontname', 'Times New Roman','FontSize',20);
