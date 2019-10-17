function varargout = CoSMoS(varargin)
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
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CoSMoS_OpeningFcn, ...
                   'gui_OutputFcn',  @CoSMoS_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before CoSMoS is made visible.
function CoSMoS_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CoSMoS (see VARARGIN)

% Choose default command line output for CoSMoS
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes CoSMoS wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CoSMoS_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% Start of the CoSMoS 

% --- Executes on selection change in popupmenu_dist.
function popupmenu_dist_Callback(hObject, eventdata, handles)
%% choose distribution 
delete('*.mat');

contents_dist = get(handles.popupmenu_dist,'String');
dist1 = contents_dist{get(handles.popupmenu_dist,'Value')};

if(strcmp(dist1, 'weibull' )||strcmp(dist1, 'Pareto Type II'))
set(handles.Shape2_dist, 'enable', 'off')
set(handles.loc_dist, 'enable', 'off')
set(handles.Shape1_dist, 'enable', 'on')

elseif strcmp(dist1, 'Burr Type III' )||strcmp(dist1, 'Burr Type XII')||strcmp(dist1, 'Generalized Gamma')
    set(handles.Shape1_dist, 'enable', 'on')
set(handles.Shape2_dist, 'enable', 'on')
set(handles.loc_dist, 'enable', 'off')

elseif(strcmp(dist1, 'Normal' ))
set(handles.loc_dist, 'enable', 'on')
set(handles.Shape2_dist, 'enable', 'off')
set(handles.Shape1_dist, 'enable', 'off')

end

% --- Executes during object creation, after setting all properties.
function popupmenu_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu_acs.
function popupmenu_acs_Callback(hObject, eventdata, handles)
%% choose autocorrelation structure

contents = get(handles.popupmenu_acs,'String');
ACS_Dist = contents{get(handles.popupmenu_acs,'Value')}; 

if(strcmp(ACS_Dist, 'G-Logarithmic' )||strcmp(ACS_Dist, 'Pareto Type II' )||strcmp(ACS_Dist, 'weibull' ))
set(handles.Shape2_acs, 'enable', 'off')
set(handles.Shape1_acs, 'enable', 'on')
set(handles.Scale_acs, 'enable', 'on')

elseif strcmp(ACS_Dist, 'Fractional Gaussian Noise' )
        set(handles.Shape2_acs, 'enable', 'off')
        set(handles.Shape1_acs, 'enable', 'on')
        set(handles.Scale_acs, 'enable', 'off')

           R = msgbox({'Enter Hurst parameter in shape text box'});
ACS_Dist='FGN';
pause(2);
delete(R);
elseif strcmp(ACS_Dist, 'Burr Type XII' )
    set(handles.Shape2_acs, 'enable', 'on')
    set(handles.Shape1_acs, 'enable', 'on')
    set(handles.Scale_acs, 'enable', 'on')

elseif strcmp(ACS_Dist, 'Markovian' )
        set(handles.Shape2_acs, 'enable', 'off')
        set(handles.Shape1_acs, 'enable', 'off')
        set(handles.Scale_acs, 'enable', 'on')
                   R = msgbox({'Enter one-lag parameter in scale text box'});
pause(2);
delete(R);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_acs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in pushbutton_generate.
function pushbutton_generate_Callback(hObject, eventdata, handles)
%% assign distribution parameters
tic
contents_dist = get(handles.popupmenu_dist,'String');
dist1 = contents_dist{get(handles.popupmenu_dist,'Value')};
x=str2double(get(handles.scale_dist,'string'));
y=str2double(get(handles.Shape1_dist,'string'));
z=str2double(get(handles.Shape2_dist,'string'));
location=str2double(get(handles.loc_dist,'string'));
 if strcmp( dist1, 'Normal' )
            param=[location  x];
        elseif strcmp(dist1, 'Pareto Type II' )||strcmp(dist1, 'weibull' )
            param=[y x];
        elseif strcmp(dist1, 'Burr Type III' )||strcmp(dist1, 'Burr Type XII' )||strcmp(dist1, 'Generalized Gamma' )
            param=[y z x];
         
 end
%% assign intermittency and compute zp
P0_int=str2double(get(handles.P0_int,'string'));
if isnan(P0_int)
            P0_int=0;
 z_p = icdf('Normal',P0_int,0,1);
        else
        P0_int=P0_int;
        z_p = icdf('Normal',P0_int,0,1);
end
%% assign acs parameters
contents = get(handles.popupmenu_acs,'String');
ACS_Dist = contents{get(handles.popupmenu_acs,'Value')};
if strcmp(ACS_Dist, 'Fractional Gaussian Noise' )
    ACS_Dist='FGN';
end
x_acs=str2double(get(handles.Scale_acs,'string'));
y_acs=str2double(get(handles.Shape1_acs,'string'));
z_acs=str2double(get(handles.Shape2_acs,'string'));
if strcmp( ACS_Dist, 'FGN' )
            Input_acs_param=[y_acs];
        elseif strcmp( ACS_Dist, 'Markovian')
            Input_acs_param=[x_acs];
elseif strcmp( ACS_Dist, 'G-Logarithmic'  )||strcmp(ACS_Dist, 'Pareto Type II' )||strcmp(ACS_Dist, 'weibull' )
    Input_acs_param=[x_acs y_acs];
elseif strcmp(ACS_Dist, 'Burr Type XII' )
Input_acs_param=[x_acs y_acs z_acs];
end
ACS_Dist_param=Input_acs_param;
%% assign sample size and lag
Xlag=str2double(get(handles.Xlag,'string'));
LTS=str2double(get(handles.LTS_tag,'string'));
TS_number=str2double(get(handles.TS_num,'string'));
if isnan(Xlag)
    Xlag=10;
end
if isnan(LTS)
    LTS=10000;
end
if isnan(TS_number)
    TS_number=1;
end
     Xlagged_S=(0:1:Xlag)';
%% compute actf variables using Eqs.8 and 9 explained in the reference
        ACTF = actpnts_t(dist1,param, P0_int);
        bc=fitactf(ACTF);
%      ACS_Dist_param=Input_acs_param;
  %% compute fitted acs values according to selected acs 
  ACS_fit=fitAcs(ACS_Dist,ACS_Dist_param,Xlagged_S);
 %% calculate gussian acs  using Eq. 27    
    ACS_Zlag=ACS_Z_lag(ACS_fit,bc);

    %% simulate synthetic timeseries using Eq. 31
    GeneratedTS=AR_p_LTS(ACS_Zlag,dist1,LTS, param,P0_int,TS_number);
    %% compute empirical and synthetic acs
    for m=1:TS_number
        ACS_Syn{1,m}=acs_syn(GeneratedTS(:,m),size(ACS_fit,1));
        [Distribution{1,m},cumf{1,m},tem{1,m}]=Exceed_pr(GeneratedTS(:,m),dist1,param);
Distribution{1,m}=[Distribution{1,m},cumf{1,m},tem{1,m}];
ACS{1,m}=[Xlagged_S(:,1),ACS_fit,ACS_Zlag, ACS_Syn{1,m}(:,2)];
    end
   
[ACS_syn_show,Emperic_exceed_show,dist_show,tem_show,cumf_show]=displ(dist1,ACS_Syn,Distribution);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   Ra = msgbox({'CoSMoS is ploting required graphs'});
cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes3)
cla(handles.axes4)
 axes(handles.axes1)
 %% Timeseries plot
plot(GeneratedTS(:,1),'k-')
hold on
ylabel('Synthetic TS','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])
xlabel('Time','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])

     set(gca,'color',[208/255,208/255,208/255]);
set(gca, 'ycolor',[208/255,208/255,208/255]);
        set(gca, 'xcolor',[208/255,208/255,208/255]);

     %% synthetic acs plot
      axes(handles.axes2)

semilogy(tem_show,Emperic_exceed_show,'k.','MarkerSize',10)
hold on
semilogy(tem_show,cumf_show,'k-')
legend('Empirical',['F_{' num2str(dist_show) '}(X)'],'Location','best','color','none','FontSize',8,'FontWeight','bold')

ylabel('Exceedance Pr','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])
xlabel('Variable','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255])
     set(gca,'color',[208/255,208/255,208/255]);
set(gca,'color',[208/255,208/255,208/255]);
set(gca, 'ycolor',[208/255,208/255,208/255]);
        set(gca, 'xcolor',[208/255,208/255,208/255])
        %% exceendance probability plot
      axes(handles.axes3)
plot(Xlagged_S(:,1),ACS_fit,'b-','LineWidth',2)
    hold on
   
plot(Xlagged_S(:,1),ACS_Zlag,'k--')
 hold on
plot(Xlagged_S(:,1),ACS_syn_show(:,2),'k.','MarkerSize',10)
xlabel('Lag (\tau)','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255]) ;
ylabel('Auto correlation','FontSize',8,'FontWeight','bold','color',[208/255,208/255,208/255]) ;
legend({'Target \rho_{x}(\tau)','Gaussian \rho_{z}(\tau)','Empirical'},'color','none','Location','best','FontSize',9)


set(gca,'fontweight','bold','fontsize',8);
     set(gca,'color',[208/255,208/255,208/255]);
set(gca,'color',[208/255,208/255,208/255]);
set(gca, 'ycolor',[208/255,208/255,208/255]);
        set(gca, 'xcolor',[208/255,208/255,208/255])

        %% autocorrelation values plot

 axes(handles.axes4)
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

hold off
PWD_main=pwd;
if ispc
    mkdir([pwd,'\Results'])
    cd([pwd,'\Results'])
else
     mkdir([pwd,'/Results'])
    cd([pwd,'/Results'])
end
delete('*.txt');
delete('*.mat');
       %% save calculated values to Reposrt text file
       Rd = msgbox({'CoSMoS is exporting results into MAT files'});
     
       SimulationInfo=savetostruct(dist1,dist_show,ACS_Dist,param,ACS_Dist_param,Xlagged_S,LTS,P0_int,TS_number);
          save GeneratedTS.mat GeneratedTS
          save Distribution.mat Distribution
          save ACS.mat ACS
          save SimulationInfo.mat SimulationInfo
          save ACTF.mat ACTF
          pause(3)
     delete(Ra);
     delete(Rd);
  cd(PWD_main);
 RunTime = toc 
  
  %% Functions called in gui

function res = acs_syn(data,Init_nLags)


Var=data;
    for i=1:Init_nLags
    Var1=Var(i:end);
    Var2=corrcoef(Var(1:(end-i+1)), Var1, 'rows','complete');
    %Auto_Corr_Var(i,1)=Var2(1,2);
    
    
    Emperic_ACS(i,1)=Var2(1,2);
    end
    
    
Xlagged_S=(0:(size( Emperic_ACS,1)-1))';

res=[Xlagged_S,Emperic_ACS];

function C = cintegral_t(dist1,RhoZ,param,P0_int)
z_p = icdf('Normal',P0_int,0,1);

lowlim=max(-8,z_p);
if strcmp( dist1, 'weibull' )
    for i=1:size(RhoZ,2)
        Fun=@(w,r)((Quantile_weibull(w,P0_int,param)).*(Quantile_weibull(r,P0_int,param)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2)))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
       
    end
elseif strcmp( dist1, 'Normal' )
    mu=param(1);
    sigma=param(2);
    for i=1:size(RhoZ,2)
        Fun=@(w,r)((norminv(Fi(w,P0_int),mu,sigma)).*(norminv(Fi(r,P0_int),mu,sigma)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2)))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
    end
elseif strcmp( dist1, 'Generalized Gamma' )
    
    for i=1:size(RhoZ,2)
        Fun=@(w,r)((Quantile_GG((((0.5*(erfc(-w/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*(Quantile_GG((((0.5*(erfc(-r/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2)))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
    end
elseif strcmp( dist1, 'Generalized Pareto' )
    k=param(1);
    sigma=param(2);
    loc=param(3);
    for i=1:size(RhoZ,2)
        Fun=@(w,r)((loc+(sigma/k).*(1-(1-(((0.5*(erfc(-w/sqrt(2)))) - P0_int)/(1-P0_int))).^k)).*(loc+(sigma/k).*(1-(1-(((0.5*(erfc(-r/sqrt(2)))) - P0_int)/(1-P0_int))).^k)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2)))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
    end
elseif strcmp( dist1, 'Pareto Type II' )
    for i=1:size(RhoZ,2)
        Fun=@(w,r)((Quantile_PII((((0.5*(erfc(-w/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*(Quantile_PII((((0.5*(erfc(-r/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2)))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
    end
elseif strcmp( dist1, 'Burr Type III' )
    
    for i=1:size(RhoZ,2)
        Fun=@(w,r)((Quantile_BIII((((0.5*(erfc(-w/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*(Quantile_BIII((((0.5*(erfc(-r/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2)))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
    end
elseif strcmp( dist1, 'Burr Type XII' )
    for i=1:size(RhoZ,2)
        Fun=@(w,r)(Quantile_BXII((((0.5*(erfc(-w/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*(Quantile_BXII((((0.5*(erfc(-r/sqrt(2)))) - P0_int)/(1-P0_int)),param)).*((1./(2*pi*sqrt(1-RhoZ(i)^2))).*exp(-(w.^2-2*RhoZ(i)*w.*r+r.^2)./(2*(1-RhoZ(i)^2))));
        C(1,i)=integral2(Fun,lowlim,8,lowlim,8);
    end
end
function bc = fitactf (res)  %% bc(2)=b; bc(1)=c; Equation 27

RhoZ=res(1,:);
RhoX=res(2,:);

x0=[10,10];
ub=[];
lb=[0.001,0];
options = optimoptions('lsqnonlin','Display','off');
fun = @(b)((((1 + b(1).*RhoX).^(1 - b(2))) - 1)/(((1 + b(1)).^(1 - b(2))) - 1))-RhoZ;
bc = lsqnonlin(fun,x0,lb,ub,options);

function res = ACS_Z_lag(ACS_Weibull_lag,b)  


res = ((((1 + b(1).*ACS_Weibull_lag).^(1 - b(2))) - 1)/(((1 + b(1)).^(1 - b(2))) - 1));


function [Empirical_Exceedance,cumf,tem]=Exceed_pr(P_Z,dist1,param)
tem=sort(nonzeros(P_Z));
N=size(tem,1);
tem1=(1:N)';
Empirical_Exceedance=1-(tem1/(N+1));
if(strcmp(dist1, 'weibull' ))
    b=param(2);
    c=param(1);
    cumf=1-wblcdf(tem,b,c);
elseif strcmp(dist1, 'Pareto Type II')
     b=param(2);
    c=param(1);
    cumf=(1+(c*tem/b)).^(-1/c);
elseif strcmp(dist1, 'Burr Type III' )
    gam1=param(1);%shape1
gam2=param(2);%shape2
a=param(3);%scale;;
cumf=1-(1+((1/gam1)*((tem/a).^(-1/gam2)))).^(-gam1*gam2);
    elseif strcmp(dist1, 'Burr Type XII')
        a=param(3); % scale parameter
gam1=param(1); % first shape parameter
gam2=param(2); % second shape parameter
cumf=(gam2*(tem/a).^(gam1)+1).^(-1/(gam1.*gam2));
elseif strcmp(dist1, 'Generalized Gamma')
a=param(3);
gam1=param(1);
gam2=param(2);
% PDF_gg=(p/(a.^d)).*(data.^(d-1)).*(exp(-(data/a)).^p)/(gamma(d/p));
cumf=gammainc(gam1/gam2,(tem/a).^gam2);
elseif strcmp(dist1, 'Normal')
    mu=param(1);% mean
sigma=param(2); % scale
% cumf=exp(-((tem/b).^c));
cumf=1-normcdf(tem,mu,sigma);
end

function res = actpnts_t(dist1,param, P0_int)


RhoZ=[0:0.1:0.9 0.95];

C=cintegral_t(dist1,RhoZ,param,P0_int);

if strcmp( dist1, 'weibull' )
    pd=makedist(dist1,'a',param(2),'b',param(1));
    Mean=mean(pd);
    Var=var(pd);
    [Var_R,Mean_R]=rev_m(Mean,Var,P0_int);
    elseif strcmp( dist1, 'Normal' )
         Mean=param(1);
    Var=param(2).^2;
        [Var_R,Mean_R]=rev_m(Mean,Var,P0_int);
elseif strcmp( dist1, 'Generalized Gamma' )
    Mean=moment_GG(1,param);
    M_x=moment_GG(2,param);
      Var=M_x-Mean.^2;
      [Var_R,Mean_R]=rev_m(Mean,Var,P0_int);
elseif strcmp( dist1, 'Generalized Pareto' )
    k=param(1);
    sigma=param(2);
    theta=param(3);
    pd=makedist('Generalized Pareto','k',param(1),'sigma',param(2),'theta',param(3));
    Mu_x_2=mean(pd).*(1-P0_int);
    Var_x_2=var(pd).*(1-P0_int)+P0_int*(1-P0_int).*(mean(pd).^2);
elseif strcmp( dist1,'Pareto Type II')
    Mean=moment_PII(1,param);
    M_x=moment_PII(2,param);
    Var=M_x-Mean.^2;
    [Var_R,Mean_R]=rev_m(Mean,Var,P0_int);
elseif strcmp( dist1, 'Burr Type III' )
    Mean=moment_BIII(1,param);
    M_x=moment_BIII(2,param);
    Var=M_x-Mean.^2;
    [Var_R,Mean_R]=rev_m(Mean,Var,P0_int);
elseif strcmp( dist1, 'Burr Type XII' )
    Mean=moment_BXII(1,param);
    M_x=moment_BXII(2,param);
     Var=M_x-Mean.^2;
   [Var_R,Mean_R]=rev_m(Mean,Var,P0_int);
end
RhoX=(C-(Mean_R.^2))/(Var_R);
res=[RhoZ;RhoX];
res(:,end) = [1 1];
% % new AR_LTS
function P_Z = AR_p_LTS(ACS_Zlag, dist,LTS,param,P0_int,TS_number)


ACS_Z_lag_from1=ACS_Zlag(2:end);
P=ones(size(ACS_Z_lag_from1,1));
 for g=1:size(ACS_Z_lag_from1,1)
     for f=1:size(ACS_Z_lag_from1,1)
             if g==f
                 P(g,f)==1;
             else
                 y=abs(g-f);
             P(g,f)=ACS_Z_lag_from1(y,1);
     end
     end
 end
 
a=inv(P)*ACS_Z_lag_from1;
%%%% White noise

var_epsilon=sqrt(1-sum(a.*ACS_Z_lag_from1));

p=size(ACS_Z_lag_from1,1);
for m=1:TS_number
Z_t1T(:,m)=random('Normal',0,var_epsilon,[p,1]);
Z_tT=Z_t1T;
for i=p+1:(LTS+p)
    for j=1:p
        mu1(j,m)=a(j).*Z_t1T(p-(j-1),m);
    end

   Z_tT(i,m)=sum(mu1(:,m))+random('Normal',0,var_epsilon);
   Z_t1T(:,m)=Z_tT((end-(p-1)):end,m);
end
Z_tTo(:,m)=Z_tT(p+1:end,m);
end
Z_tT=Z_tTo;
pd = makedist('Normal','mu',0,'sigma',1);
CDF_Z=cdf(pd, Z_tT);
P0=P0_int;
for m=1:size(Z_tT,2)
for t=1:size(Z_tT,1)
if CDF_Z(t,m)-P0<0
    P_Z(t,m)=0;
elseif strcmp(dist, 'weibull' ) 
    P_Z1(t,m)=(CDF_Z(t,m)-P0)/(1-P0);
    P_Z(t,m)=wblinv(P_Z1(t,m),param(2),param(1));
    elseif strcmp(dist, 'Normal' ) 
    P_Z1(t,m)=(CDF_Z(t,m)-P0)/(1-P0);
    P_Z(t,m)=norminv(P_Z1(t,m),param(1),param(2));
elseif strcmp(dist, 'Pareto Type II' )
    P_Z1(t,m)=(CDF_Z(t,m)-P0)/(1-P0);
    P_Z(t,m)=Quantile_PII(P_Z1(t,m),param);
     elseif strcmp(dist, 'Burr Type XII' )
    P_Z1(t,m)=(CDF_Z(t,m)-P0)/(1-P0);
    P_Z(t,m)=Quantile_BXII(P_Z1(t,m),param);
    elseif strcmp(dist, 'Burr Type III' )
    P_Z1(t,m)=(CDF_Z(t,m)-P0)/(1-P0);
    P_Z(t,m)=Quantile_BIII(P_Z1(t,m),param);
    elseif strcmp(dist, 'Generalized Gamma' )
    P_Z1(t,m)=(CDF_Z(t,m)-P0)/(1-P0);
    P_Z(t,m)=Quantile_GG(P_Z1(t,m),param);
end
end
end
function InvPII=Quantile_PII(data,param)
InvPII=((((1-data).^(-param(1)))-1)./param(1))*param(2);

function M_PII=moment_PII(r,param)
M_PII=(((param(2)./param(1)).^r * gamma((1/param(1)) - r) .* gamma(r+ 1))/gamma(1/param(1)));
function [Var_R,Mean_R]=rev_m(Mean,Var,P0_int)
Var_R=Var.*(1-P0_int)+P0_int*(1-P0_int).*(Mean.^2);
Mean_R=Mean.*(1-P0_int);
function InvBIII=Quantile_BIII(data,param)
InvBIII=param(1,3)*(param(1,1)*(((data).^(-1/(param(1)*param(1,2))))-1)).^(-param(1,2));
function InvBXII=Quantile_BXII(data,param)
InvBXII=param(3)*(-(1/param(2))*(1-(1-data).^(-param(1)*param(2)))).^(1/param(1));
function InvGG=Quantile_GG(data,param)
InvGG= param(3)*(gammaincinv(data,(param(1)/(param(2))))).^(1/param(2)); 
function M_BurrXII=moment_BXII(r,param)    
M_BurrXII=(param(3)^r * param(2).^(-1 - r/param(1)) * beta((r + param(1))/param(1), (1 - r * param(2))/(param(1) * param(2)))/param(1));
function M_BurrGG=moment_GG(r,param)   
M_BurrGG=param(3).^r.*(gamma((r/param(2)) + (param(1)/param(2)))/gamma(param(1)/param(2)));

function Inv_weibull=Quantile_weibull(w,P0_int,param)
    gam=param(1);
    lam=param(2);
Inv_weibull=((lam)*(-log(1-Fi(w,P0_int))).^(1/gam));

function M_BurrIII=moment_BIII(r,param)
M_BurrIII=((param(3)^r * gamma((r + param(1)) * param(2)) * gamma(1 - r * param(2)))/(param(1)^(r * param(2)) * gamma(param(1) * param(2))));
function u=Fi(r,P0_int)
u=(((0.5*(erfc(-r/sqrt(2)))) - P0_int)/(1-P0_int));
function res = acs(data)
% 
N=length(data);
Init_nLags=length(data)/3;
Init_nLags=round(Init_nLags);
Var=data;
    for i=1:Init_nLags
    Var1=Var(i:end);
    Var2=corrcoef(Var(1:(end-i+1)), Var1, 'rows','complete');
    %Auto_Corr_Var(i,1)=Var2(1,2);
    if Var2(1,2)<0.05
        break
    end
    
    Emperic_ACS(i,1)=Var2(1,2);
    end
    
    
Xlagged_S=(0:(size( Emperic_ACS,1)-1))';

res=[Xlagged_S,Emperic_ACS];
function  [ACS_syn_show,Emperic_exceed_show,dist_show,tem_show,cumf_show]=displ(dist1,Auto_Corr_Var_synthetic,Empirical_Exceedance)
ACS_syn_show=Auto_Corr_Var_synthetic{1,1};
Emperic_exceed_show=Empirical_Exceedance{1,1}(:,1);
tem_show=Empirical_Exceedance{1,1}(:,3);
cumf_show=Empirical_Exceedance{1,1}(:,2);
if strcmp(dist1, 'Pareto Type II' )
      dist_show='PII';
      elseif strcmp(dist1, 'Burr Type XII' )
           dist_show='BXII';
            elseif strcmp(dist1, 'Burr Type III' )
           dist_show='BIII';
            elseif strcmp(dist1, 'Generalized Gamma' )
           dist_show='GG';
           elseif strcmp(dist1, 'Normal' )
           dist_show='Norm';
           elseif strcmp(dist1, 'weibull' )
           dist_show='Weib';
end
%%ghabli
% function  [ACS_syn_show,Emperic_exceed_show,dist_show,tem_show,cumf_show]=displ(dist1,Auto_Corr_Var_synthetic,Empirical_Exceedance,tem,cumf)
% ACS_syn_show=Auto_Corr_Var_synthetic{1,1};
% Emperic_exceed_show=Empirical_Exceedance{1,1};
% tem_show=tem{1,1};
% cumf_show=cumf{1,1};
% if strcmp(dist1, 'Pareto Type II' )
%       dist_show='PII';
%       elseif strcmp(dist1, 'Burr Type XII' )
%            dist_show='BXII';
%             elseif strcmp(dist1, 'Burr Type III' )
%            dist_show='BIII';
%             elseif strcmp(dist1, 'Generalized Gamma' )
%            dist_show='GG';
%            elseif strcmp(dist1, 'Normal' )
%            dist_show='Norm';
%            elseif strcmp(dist1, 'weibull' )
%            dist_show='Weib';
% end
function  ACS_fit=fitAcs(ACS_Dist,ACS_Dist_param, Xlagged_S)
if strcmp(ACS_Dist, 'weibull' )
        ACS_Weibull=@(ACS_Dist_param) exp(-(Xlagged_S./ACS_Dist_param(1)).^(ACS_Dist_param(2)));
        ACS_fit=ACS_Weibull(ACS_Dist_param);
    elseif strcmp(ACS_Dist, 'Pareto Type II' )
        ACS_ParetoII=@(ACS_Dist_param) ((1+ACS_Dist_param(2).*(Xlagged_S./ACS_Dist_param(1))).^(-1/ACS_Dist_param(2)));
        ACS_fit=ACS_ParetoII(ACS_Dist_param);
    elseif strcmp(ACS_Dist, 'Burr Type XII' )
        ACS_BurrXII=@(ACS_Dist_param)((1+ACS_Dist_param(2).*(Xlagged_S./ACS_Dist_param(1)).^(ACS_Dist_param(3))).^(-1/(ACS_Dist_param(3).*ACS_Dist_param(1))));
        ACS_fit=ACS_BurrXII(ACS_Dist_param);
    elseif strcmp(ACS_Dist, 'G-Logarithmic' )
        ACS_Gen_Logarithmic=@(ACS_Dist_param)((1+log(1+ACS_Dist_param(2).*(Xlagged_S./ACS_Dist_param(1)))).^(-1/ACS_Dist_param(2)));
        ACS_fit= ACS_Gen_Logarithmic(ACS_Dist_param);
        elseif strcmp(ACS_Dist, 'Markovian' )
        ACS_Marcov=@(ACS_Dist_param)(ACS_Dist_param(1).^(Xlagged_S));
        ACS_fit= ACS_Marcov(ACS_Dist_param);
            elseif strcmp(ACS_Dist, 'FGN' )
FGN=@(ACS_Dist_param) (0.5*(abs(Xlagged_S+1).^(2*ACS_Dist_param(1)) + abs(Xlagged_S-1).^(2*ACS_Dist_param(1)) - 2*abs(Xlagged_S).^(2*ACS_Dist_param(1))));
ACS_fit= FGN(ACS_Dist_param);
end

  function P=savetostruct(dist1,dist_show,ACS_Dist,param,ACS_Dist_param,Xlagged_S,LTS,P0_int,TS_number)
P(1,1).InputParameters='Marginal Distribution';
P(2,1).InputParameters='Autocorrelation Structure';
P(3,1).InputParameters='Intermittency';
P(4,1).InputParameters='Sample Size';
P(5,1).InputParameters='No of TS';
P(6,1).InputParameters='Max Lag';
if (strcmp(dist1, 'weibull' )||strcmp(dist1,  'Pareto Type II'))&&(strcmp(ACS_Dist, 'G-Logarithmic' )||strcmp(ACS_Dist, 'Pareto Type II' )||strcmp(ACS_Dist, 'weibull' ))
P(1,1).Val={'Name','Short Form','Scale','Shape';dist1,dist_show,param(2),param(1)};
P(2,1).Val={'Name','Scale','Shape';ACS_Dist,ACS_Dist_param(1),ACS_Dist_param(2)};
elseif (strcmp(dist1, 'Normal'))&&(strcmp(ACS_Dist, 'G-Logarithmic' )||strcmp(ACS_Dist, 'Pareto Type II' )||strcmp(ACS_Dist, 'weibull' ))
P(1,1).Val={'Name','Short Form','Location','Scale';dist1,dist_show,param(2),param(1)};
P(2,1).Val={'Name','Scale','Shape';ACS_Dist,ACS_Dist_param(1),ACS_Dist_param(2)};
elseif(strcmp(dist1, 'weibull' )||strcmp(dist1,'Pareto Type II' ))&&(strcmp(ACS_Dist, 'FGN' ))
P(1,1).Val={'Name','Short Form','Scale','Shape';dist1,dist_show,param(2),param(1)};
P(2,1).Val={'Name','Hurst';ACS_Dist,ACS_Dist_param(1)};
elseif(strcmp(dist1, 'Normal' ))&&(strcmp(ACS_Dist, 'FGN' ))
P(1,1).Val={'Name','Short Form','Location','Scale';dist1,dist_show,param(2),param(1)};
P(2,1).Val={'Name','Hurst';ACS_Dist,ACS_Dist_param(1)};
elseif  (strcmp(dist1, 'weibull' )||strcmp(dist1,'Pareto Type II' ))&&(strcmp(ACS_Dist, 'Markovian' ))
P(1,1).Val={'Name','Short Form','Scale','Shape';dist1,dist_show,param(1,2),param(1,1)};
P(2,1).Val={'Name','1-lag';ACS_Dist,ACS_Dist_param(1)};
elseif(strcmp(dist1, 'Normal' ))&&(strcmp(ACS_Dist, 'Markovian' ))
P(1,1).Val={'Name','Short Form','Location','Scale';dist1,dist_show,param(2),param(1)};
P(2,1).Val={'Name','1-lag';ACS_Dist,ACS_Dist_param(1)};
elseif  (strcmp(dist1, 'Burr Type XII' )||strcmp(dist1, 'Burr Type III' )||strcmp(dist1, 'Generalized Gamma' ))&&(strcmp(ACS_Dist, 'FGN' ))
P(1,1).Val={'Name','Short Form','Scale','Shape1','Shape2';dist1,dist_show,param(1,3),param(1,1),param(1,2)};
P(2,1).Val={'Name','Hurst';ACS_Dist,ACS_Dist_param(1)};
elseif  (strcmp(dist1, 'Burr Type XII' )||strcmp(dist1, 'Burr Type III' )||strcmp(dist1, 'Generalized Gamma' ))&&(strcmp(ACS_Dist, 'Markovian' ))
P(1,1).Val={'Name','Short Form','Scale','Shape1','Shape2';dist1,dist_show,param(1,3),param(1,1),param(1,2)};
P(2,1).Val={'Name','1-lag';ACS_Dist,ACS_Dist_param(1)};
elseif  (strcmp(dist1, 'Burr Type XII' )||strcmp(dist1, 'Burr Type III' )||strcmp(dist1, 'Generalized Gamma' ))&&(strcmp(ACS_Dist, 'Pareto Type II' )||strcmp(ACS_Dist, 'weibull' )||strcmp(ACS_Dist, 'G-Logarithmic' ))
P(1,1).Val={'Name','Short Form','Scale','Shape1','Shape2';dist1,dist_show,param(1,3),param(1,1),param(1,2)};
P(2,1).Val={'Name','Scale','Shape';ACS_Dist,ACS_Dist_param(1),ACS_Dist_param(2)};
elseif  (strcmp(dist1, 'Burr Type XII' )||strcmp(dist1, 'Burr Type III' )||strcmp(dist1, 'Generalized Gamma' ))&&(strcmp(ACS_Dist, 'Burr Type XII'))
P(1,1).Val={'Name','Short Form','Scale','Shape1','Shape2';dist1,dist_show,param(1,3),param(1,1),param(1,2)};
P(2,1).Val={'Name','Scale','Shape1','Shape2';ACS_Dist,ACS_Dist_param(1),ACS_Dist_param(2),ACS_Dist_param(3)};
elseif (strcmp(dist1, 'weibull' )||strcmp(dist1,'Pareto Type II'))&&(strcmp(ACS_Dist, 'Burr Type XII'))
P(1,1).Val={'Name','Short Form','Scale','Shape';dist1,dist_show,param(1,2),param(1,1)};
P(2,1).Val={'Name','Scale','Shape1','Shape2';ACS_Dist,ACS_Dist_param(1),ACS_Dist_param(2),ACS_Dist_param(3)};
elseif (strcmp(dist1, 'Normal'))&&(strcmp(ACS_Dist, 'Burr Type XII'))
P(1,1).Val={'Name','Short Form','Location','Scale';dist1,dist_show,param(1,2),param(1,1)};
P(2,1).Val={'Name','Scale','Shape1','Shape2';ACS_Dist,ACS_Dist_param(1),ACS_Dist_param(2),ACS_Dist_param(3)};
end
P(3,1).Val=P0_int;
P(4,1).Val=LTS;
P(5,1).Val=TS_number;
P(6,1).Val=Xlagged_S;
  %% Gui components

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function scale_dist_Callback(hObject, eventdata, handles)
% hObject    handle to scale_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of scale_dist as text
%        str2double(get(hObject,'String')) returns contents of scale_dist as a double


% --- Executes during object creation, after setting all properties.
function scale_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scale_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Shape1_dist_Callback(hObject, eventdata, handles)
% hObject    handle to Shape1_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shape1_dist as text
%        str2double(get(hObject,'String')) returns contents of Shape1_dist as a double


% --- Executes during object creation, after setting all properties.
function Shape1_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shape1_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Shape2_dist_Callback(hObject, eventdata, handles)
% hObject    handle to Shape2_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shape2_dist as text
%        str2double(get(hObject,'String')) returns contents of Shape2_dist as a double

% --- Executes during object creation, after setting all properties.
function Shape2_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shape2_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function LTS_tag_Callback(hObject, eventdata, handles)
% hObject    handle to LTS_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LTS_tag as text
%        str2double(get(hObject,'String')) returns contents of LTS_tag as a double


% --- Executes during object creation, after setting all properties.
function LTS_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LTS_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function P0_int_Callback(hObject, eventdata, handles)
% hObject    handle to P0_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P0_int as text
%        str2double(get(hObject,'String')) returns contents of P0_int as a double


% --- Executes during object creation, after setting all properties.
function P0_int_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P0_int (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Scale_acs_Callback(hObject, eventdata, handles)
% hObject    handle to Scale_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Scale_acs as text
%        str2double(get(hObject,'String')) returns contents of Scale_acs as a double


% --- Executes during object creation, after setting all properties.
function Scale_acs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Scale_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Shape1_acs_Callback(hObject, eventdata, handles)
% hObject    handle to Shape1_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shape1_acs as text
%        str2double(get(hObject,'String')) returns contents of Shape1_acs as a double


% --- Executes during object creation, after setting all properties.
function Shape1_acs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shape1_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Shape2_acs_Callback(hObject, eventdata, handles)
% hObject    handle to Shape2_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Shape2_acs as text
%        str2double(get(hObject,'String')) returns contents of Shape2_acs as a double


% --- Executes during object creation, after setting all properties.
function Shape2_acs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Shape2_acs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N_Callback(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N as text
%        str2double(get(hObject,'String')) returns contents of N as a double


% --- Executes during object creation, after setting all properties.
function N_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes during object creation, after setting all properties.
function axes3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes3


% --- Executes during object creation, after setting all properties.
function axes4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes4



function loc_dist_Callback(hObject, eventdata, handles)
% hObject    handle to loc_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loc_dist as text
%        str2double(get(hObject,'String')) returns contents of loc_dist as a double


% --- Executes during object creation, after setting all properties.
function loc_dist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loc_dist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1


% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function Xlag_Callback(hObject, eventdata, handles)
% hObject    handle to Xlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Xlag as text
%        str2double(get(hObject,'String')) returns contents of Xlag as a double


% --- Executes during object creation, after setting all properties.
function Xlag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Xlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function TS_num_Callback(hObject, eventdata, handles)
% hObject    handle to TS_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TS_num as text
%        str2double(get(hObject,'String')) returns contents of TS_num as a double


% --- Executes during object creation, after setting all properties.
function TS_num_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TS_num (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
