% Dati
% Data
function [epsilon_0,nu,sigma,Sigma_rod,V_cyto,theta,alpha,tol_fix,t_fin,n_step_t,t_step,u_tent,v_tent,v_dark,tol_stat,flag_Ca_clamp, ...
    n_step_R,lambda,mu,nu_RG,K1,K2,K3,K4,M,Rec_tot,RK_tot,I,I_hys,Phi,E_tot,G_tot,RGS9_tot,k_GE,k_cat,k_b,k_f,k_hyd,k_st,alpha_max,alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    B_ca,F,j_cg_max,f_ca,m_cg,K_cg_max,K_cg_min,K_caM,n_caM,j_ex_sat,K_ex]=data



% PARAMETRI GEOMETRICI
% GEOMETRICAL PARAMETERS

% raggio del bastoncello
% rod radius
R=0.7-0.015; % [\mu m]
% altezza del bastoncello
% rod height
H=23.6; % [\mu m]



% spessore dei dischi
% disc thickness
epsilon_0=14.5e-3; %  [\mu m]

% rapporto fra spessore degli interdischi e spessore dei dischi
% ratio between interdisc and disc thickness 
nu=14.5/14.5; % []

% rapporto fra spessore dell'outer shell e spessore dei dischi
% ratio between outer-shell and disc thickness
sigma=15/14.5; % []



% DATI SULLE INCISURE
% DATA ON INCISURES





% % I dati sulle incisure sono allocati nella matrice inc, con n_inc colonne
% % ogni colonna ha tre componenti, che sono la lunghezza dell'incisura; la posizione angolare; l'apertura
% % la lunghezza di ciascuna incisura deve essre minore del raggio
% % le anomalie a cui sono posizionate le incisure devono essere tutte positive 
% % e crescenti in senso antiorario al crescere dell'indice di incisura (anche oltre 2*pi se necessario)
% % se sono presenti solo due incisure, l'angolo dalla prima alla seconda deve essere minore del piatto
% % l'angolo di apertura (beanza) delle incisure, assunte essere settori
% % circolari, è misurato rispetto al verice dell'incisura
% % MOUSE
% n_inc=1;
% % width of a single incisure
% l_b=0.2593;
% % total area of incisures
% % A_tot=0.80;
% % length of a single incisure
% l_r=0.3111;
% % angle of a single incisure
% theta_inc=l_b/l_r;
% % generate the variable inc
% inc=[l_r*ones(1,n_inc); 0:2*pi/n_inc:2*pi*(1-1/n_inc); theta_inc*ones(1,n_inc)];




% calcolo dei parametri geometrici del ROD che entrano nel modello GWS
% a_inc=sum(inc(1,:)^2*inc(3,:)/2);
% a_pat=2*pi*R*sigma*epsilon_0+a_inc;
a_int=nu/(1+nu)*pi*R^2;
%a_tot=a_int+a_pat;

Sigma_rod=2*pi*R*H;
%V_cyto=a_tot*H;   
disp('ATTENZIONE ALLA DEFINIZIONE DI V_CYTO')
V_cyto=a_int*H;   





% parametro theta nell'integrazione temporale
% theta-parameter of the theta-method
theta=0.5; % []

% parametro di rilassamento nell'iterazione sul punto fisso: deve essere 0<alpha<=1 (alpha=1 per eliminare il rilassamento)
% relaxation parameter in the fixed-point iteration: it must be 0<alpha<=1   (alpha=1 means no relaxation)
alpha=1; % []

% errore relativo sul punto fisso nel prolema homogeneizzato
% relative error on the fixed-point iteration of the homogenized problem 
tol_fix=1e-5; % []


% durata della simulazione
% time horizon
% t_fin=1.5; % [s]
t_fin=10; % [s]

% numero di passi di integrazione
% number of integration steps
n_step_t=1000000;
% n_step_t=1000;

% lunghezza dello step temporale
t_step=t_fin/n_step_t;


% DATI INIZIALI DI TENTATIVO PER LO STEADY-STATE
% INITIAL VALUES TO BE USED IN SEARCHING FOR THE STEADY-STATE

% valore iniziale di tentativo del cGMP
% trial cGMP value
u_tent=3;  % [\mu M]

% valore iniziale di tentativo  del Ca2+ 
% trial Ca2+ value
v_tent=1; % [\mu M]

v_dark=0.344;

% tolleranza sui flussi da rendere nulli per calcolare lo stato stazionario
% maximum unbalancing of fluxes allowed in the steady-state 
tol_stat=1e-8; % [(\mu m) (\mu M)/s = 10^(-9) mole/(m^2 s)]



% flag_Ca_clamp=false: modello standard; flag_model=true modello con clamp del calcio
% flag_Ca_clamp=false: standard model; flag_model=true: Calcium clamp model
flag_Ca_clamp=false;



% DATI SULLA R*
% DATA ON R*


% numero di passi di shutoff della R*
% number of steps required for a complete shut off of R*
n_step_R=7;


% tassi di fosforilazione ed attività catalitiche della R
lambda_0=10.5;
%%%%%%%%%%%%%%%%%%%lambda_0=7;
mu_0=60;



lambda=lambda_0*(n_step_R-1:-1:0);
mu=mu_0*[0 0 0 1 1 1 1];


% Parametri per descrivere il decadimento dell'attività catalitica della
% rodopsina
kv=0.5;
nu_RG_max=330;
%%%%%%%%%%%%%%%%%%%%%%%%%%%nu_RG_max=500;
nu_RG=nu_RG_max*exp(-kv*(0:n_step_R-1));

% Feedback del calcio sulla rodopsina
K1=4.5; %muM
K2=230; %muM
K3=3.4; %muM
K4=K3;
M=6000;  %muM
Rec_tot=34;  %muM
RK_tot=7;   % muM




% Intensità della corrente pregressa, supposta costante
I=100*0/V_cyto;   % fotoisomerizzazioni per micron cubo per secondo

% Storia della corrente di background
I_hys=0*100/V_cyto*ones(1,n_step_t+1); % fotoisomerizzazioni per micron cubo per secondo
%I_hys(1:250001)=100/V_cyto*ones(1,250001);


% Intensità dei flash aggiuntivi. Per ogni flash si mette una
% coppia di numeri, il tempo del flash (in secondi) e la sua intensità (in
% fotoisomerizzazioni totali prodotte).
% Phi_flash=[0 200 0.2 200];
%Phi_flash=[0 1000  1.2  5000];
 Phi_flash=[0 120000];

 % Crea il vettore Phi per l'output, lungo n_step_t, che ha tutti zeri
 % tranne agli indici corrispondenti alla emissione del flash (indice k non
 % nullo implica flash emesso al tempo (k-1)*t_step
 Phi=zeros(1,n_step_t);
 time_flash=Phi_flash(1:2:end);
 val_flash=Phi_flash(2:2:end)/V_cyto;
 
 indici_flash=fix(time_flash/t_step)+1;
 Phi(indici_flash)=val_flash;
 
 

 
 


% DATI SULLA G*
% DATA ON G*

% tasso di produzione di E* per ogni G* per densità unitaria (#/\mu m^3) di E spenta
% production rate of E* per each G* per unit density (#/\mu m^3)of basal E
% k2 in shraiman paper
k_GE=1*1/2*nu*epsilon_0; % [\mu m^2/s]

% Densità superficiale di G totali
G_s=2500;

% DATI SULLA E*
% DATA ON E*

% Densità superficiale di PDE totali
% nota: la densità di subunit è il doppio di PDE_s
% superficial density of total PDE
% note: subunit density is 2*PDE_s
PDE_s=750;  % [molecole/(\mu m)^2]


% Calcolo delle E e G totali per unità di volume
    E_tot=2*PDE_s/(1/2*nu*epsilon_0);
    G_tot=G_s/(1/2*nu*epsilon_0);
    
    
    nu_RG=nu_RG/G_tot;

    
    RGS9_tot=185/(1/2*nu*epsilon_0);



% Tasso di reazione
k_cat=52.8/9;

k_b=13.8;

k_f=0.07*(1/2*nu*epsilon_0)*4;

%k_f=0.07*(1/2*nu*epsilon_0)*2;



%k_f=0.07*(1/2*nu*epsilon_0);









%k_E=4*k_E;
% Tasso di idrolisi del cGMP da parte della PDE attiva al buio (spenta)
% rate of hydrolysis of cGMP by dark-activated PDE
% k_hyd=7e-5;  % [(\mu m)^3/(molecole s)]
Beta_dark=2.9;
k_hyd=0.5*nu*epsilon_0*Beta_dark/PDE_s;  % [(\mu m)^3/(molecole s)]

% % attivazione della PDE spenta = k_hyd*[PDE]_s = flusso
% gamma_0= k_hyd*PDE_s;  %  [\mu m s^(-1)]

% Costante di idrolisi del cGMP da parte della PDE attivata
% rate of hydrolysis of cGMP by light-activated PDE
k_st=0.9; % [(\mu m)^3/molecole s]






% DATI SULLA CICLASI
% DATA ON CYCLASE
 
% Tasso massimo di sintesi di cGMP
% maximum rate of production of cGMP
alpha_max=50; % [\mu M s^{-1}]
% alpha_max=76.5; % [\mu M s^{-1}]

% Tasso minimo di sintesi di cGMP
% minumum rate of production of cGMP
 alpha_min=0.02*alpha_max; % [\mu M s^{-1}]
%alpha_min=5.5; % [\mu M s^{-1}]

% Dati sulla ciclasi con il modello di Makino
A1=0.4*0;
A2=0.6*0;
KK1=0.133;
KK2=0.047;
m1=2.1;
m2=1.9;

% L'attività miniuma della ciclase è legata solo alla presenza del GCAP1
%alpha_min=alpha_min*(A1>0);



% DATI COMUNI AI CANALI DEL CALCIO
% COMMON DATA ON CALCIUM CHANNELS

% Buffer del Ca2+, assunto costante
% buffer of Ca2+, here assumed to be constant
B_ca=20; % []

% Costante di Faraday
% Faraday constant
F=96500/1e21; % [C / (\mu M (\mu m)^3)]


% DATI SUL CANALE cGMP-OPERATO
% DATA ON THE cGMP-OPERATED CHANNEL

% Corrente di scambio massima
% maximum exchange current
% j_cg_max=7000e-12; % [A]
j_cg_max=3550e-12; % [A]

% Frazione di corrente attivata dal c_GMP costituita da ioni Ca2+
% fraction of cGMP-activated current given by Ca2+
% f_ca=0.17; % []
f_ca=0.06; % []

% Esponente di Hill
% Hill coefficient
% m_cg=2.5; % []
m_cg=3.5; % []

% Concentrazione di cGMP alla metà dell'apertura dei canali
% half-maximum-activation concentration, here assumed to be constant
%%%K_cg=20;%32; % [\mu M]

K_cg_max=32;
K_cg_min=13;
K_caM=60;
n_caM=2;


% K_cg_max=32;
% K_cg_min=32;
% K_caM=60;
% n_caM=2;




% DATI SULLO SCAMBIATORE
% DATA ON Ca2+ EXCHANGER

% Corrente di scambio alla saturazione
% exchange current at saturation
% j_ex_sat=17e-12; % [A]
j_ex_sat=1.8e-12; % [A]

% Concentrazione di calcio alla metà dell'apertura dei canali
% half-maximum-activation concentration
% K_ex=1.5; %[\mu M]
K_ex=1.6; %[\mu M]


return
