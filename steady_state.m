% soluzione numerica del sistema non lineare che descrive la soluzione stazionaria
function [u_ss,v_ss,E_ss,GE_ss,G_ss,G_s_ss,RGS9_ss,RGS9GE_ss,R_ss,RK_dark]=steady_state(u_tent,v_tent,v_dark,K1,K2,K3,K4,M,Rec_tot,RK_tot, ...
    n_step_R,mu,lambda,I,E_tot,G_tot,RGS9_tot,k_GE,nu_RG,k_b,k_f,k_cat,tol_fix,tol_stat,k_hyd,k_st,alpha_max,alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    B_ca,F,j_cg_max,f_ca,m_cg,K_cg_max,K_cg_min,K_caM,n_caM,j_ex_sat,K_ex,V_cyto)




% Il numero di rodopsine in ciascun stato di fosforilazione sono date dalla
% risoluzione di un sistema algebrico risolto per induzione.
% R_j=lambda_j-1R_j-1/(mu_j+lambda_j), essendo R_1=I/(mu_1+lambda_1). Note
% le R_j in funzione del calcio si calcola il livello di equilibrio di cGMP
% e Ca e reitera sul valore di Ca ottenuto fino a che è sufficientemente
% vicino a quello con cui si sono calcolate le R_j.






% Il primo valore di tentativo del calcio si assume pari al valore dark
% passato alla routine
u_ss=u_tent;
v_ss=v_tent;

% Vettore contenente le R* in ciascun stato di fosforilazione
R_ss=zeros(n_step_R,1);

% ritiene in lambd_max il valore a calcio dark delle lambda
lambda_max=lambda;


% Calcola RK allo stato dark

    Ca=v_dark;
    
    
    % I tassi ti fosforilazione lambda dipendono dal calcio, i tassi di legame con l'arrestina sono indipendenti dal calcio.    
    C1=(Ca/K1)^2*(1/K3+1/K4*M/K2)*Rec_tot;
    C2=1+(Ca/K1)^2*(1+M/K2);


    % Calcolo di Rec/Rec_tot
    aa=C1*C2;
    bb=(C1*(RK_tot/Rec_tot-1)+C2);
    Rec_Rectot=(-bb+sqrt(bb^2+4*aa))/(2*aa);

    RK_dark=RK_tot*(1+C1*Rec_Rectot)^(-1);


% Iterazione sul calcio

Ca=0;
while abs(v_ss-Ca)>tol_fix
  
    % Si assume come valore di tentativo del calcio l'ultimo valore di
    % calcio calcolato
    
     
    Ca=v_ss;
    
    
    % I tassi ti fosforilazione lambda dipendono dal calcio, i tassi di legame con l'arrestina sono indipendenti dal calcio.    
    C1=(Ca/K1)^2*(1/K3+1/K4*M/K2)*Rec_tot;
    C2=1+(Ca/K1)^2*(1+M/K2);


    % Calcolo di Rec/Rec_tot
    aa=C1*C2;
    bb=(C1*(RK_tot/Rec_tot-1)+C2);
    Rec_Rectot=(-bb+sqrt(bb^2+4*aa))/(2*aa);


    RK=RK_tot*(1+C1*Rec_Rectot)^(-1);
    
   
    
    lambda=lambda_max*RK/RK_dark;
 

   

    % Calcolo delle R*
    R_ss(1)=I/(mu(1)+lambda(1));
    for cont=2:n_step_R
        R_ss(cont)=lambda(cont-1)*R_ss(cont-1)/(mu(cont)+lambda(cont));
    end
   
 
    
    
    % Valori di primo tentativo per E_ss, GE_ss, G_ss, G_s_ss, RGS9GE_ss
    x0=[0,0,0];
    
    % Si passa a fsolve i logaritmi di x0
    x0=log(x0);
    % Calcolo della G_st e E_st
    % il numero di subunit è il doppio delle molecole di PDE
    x = fsolve(@fun_GE, x0, optimset('TolFun', tol_stat,'MaxFunEvals',1e9), RGS9_tot, G_tot, E_tot, nu_RG, k_GE, k_b, k_f, k_cat, R_ss);
    % Soluzione in termini di G e E
    GE_ss=exp(x(1));
    G_s_ss=exp(x(2));
    RGS9GE_ss=exp(x(3));
    
    G_ss=G_tot-GE_ss-G_s_ss-RGS9GE_ss;
    E_ss=E_tot-GE_ss-RGS9GE_ss;
    RGS9_ss=RGS9_tot-RGS9GE_ss;


  
    % Valori di primo tentativo per cGMP e Ca
    x0=[u_ss,v_ss];
    % si passa a fsolve i logaritmi di x0
    x0=log(x0);
    % Calcolo di cGMP e Ca
    x = fsolve(@fun_cGCa, x0, optimset('TolFun', tol_stat,'MaxFunEvals',1e9), k_hyd, E_tot, alpha_max, alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg_max,K_cg_min,K_caM,n_caM, j_ex_sat, K_ex, V_cyto, k_st, GE_ss); 

    % mette la soluzione in u_ss e v_ss
    u_ss=exp(x(1));
    v_ss=exp(x(2));
    
end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function y = fun_GE(x, RGS9_tot, G_tot, E_tot, nu_RG, k_GE, k_b, k_f, k_cat, R_ss)
    % estrae le incognite u,v da x scrivendole in forma esponenziale per
    % evitare numeri negativi
    GE_ss=exp(x(1));
    G_s_ss=exp(x(2));
    RGS9GE_ss=exp(x(3));
    
    G_ss=G_tot-GE_ss-G_s_ss-RGS9GE_ss;
    E_ss=E_tot-GE_ss-RGS9GE_ss;
    RGS9_ss=RGS9_tot-RGS9GE_ss;

    % output
    y=zeros(1,3);

    % prima equazione:
    % equazione per la G_st
    y(1)=G_ss*nu_RG*R_ss-k_GE*GE_ss*G_s_ss;
    y(2)=k_f*RGS9_ss*GE_ss-k_b*RGS9GE_ss-k_cat*RGS9GE_ss;
    y(3)=-k_GE*E_ss*G_s_ss+k_cat*RGS9GE_ss;
    
    
    
    
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = fun_cGCa(x, k_hyd, E_tot, alpha_max, alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    B_ca, F, j_cg_max, f_ca, m_cg, K_cg_max,K_cg_min,K_caM,n_caM, j_ex_sat, K_ex, V_cyto, k_st, GE_ss)
    % estrae le incognite u,v da x scrivendole in forma exponenziale per
    % evitare numeri negativi
    u=exp(x(1));
    v=exp(x(2));

    % output
    y=zeros(1,2);

    % prima equazione:
    % bilancio dei flussi del cGMP sulle facce dei dischi
    y(1)=-k_hyd*(E_tot-GE_ss)/2*u-k_st*GE_ss*u+...
        (alpha_min+(alpha_max-alpha_min)*(A1/(1+(v/KK1)^m1)+A2/(1+(v/KK2)^m2)));

    % seconda equazione:
    % bilancio dei flussi del Ca2+ sulla membrana plasmatica
    K_cg=K_cg_max + (K_cg_min-K_cg_max)/(1+(v/K_caM)^n_caM);

    y(2)=-j_ex_sat/(V_cyto*B_ca*F)*v/(K_ex+v)+...
        j_cg_max*f_ca/(V_cyto*B_ca*F*2)*u^m_cg/(K_cg^m_cg+u^m_cg);


return
