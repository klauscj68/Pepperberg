% Programma per la soluzione del problema GWS della diffusione dei secondi messaggeri 
% nel cytosol del bastoncello (problema della fototrasduzione) con
% adattamento



close all;



% legge i dati
[epsilon_0,nu,sigma,Sigma_rod,V_cyto,theta,alpha,tol_fix,t_fin,n_step_t,t_step,u_tent,v_tent,v_dark,tol_stat,flag_Ca_clamp, ...
    n_step_R,lambda,mu,nu_RG,K1,K2,K3,K4,M,Rec_tot,RK_tot,I,I_hys,Phi,E_tot,G_tot,RGS9_tot,k_GE,k_cat,k_b,k_f,k_hyd,k_st,alpha_max,alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    B_ca,F,j_cg_max,f_ca,m_cg,K_cg_max,K_cg_min,K_caM,n_caM,j_ex_sat,K_ex]=data;
    




% Calcolo della soluzione stazionaria. Le incognite sono Rj, G*, E*, cGMP e
% Ca (tutte grandezze volumiche). Si assume come dato una intensità di luce
% che agisce in maniera costante da tempo sufficiente affinché sia stato
% raggiunto uno stato stazionario
[u_ss,v_ss,E_ss,GE_ss,G_ss,G_s_ss,RGS9_ss,RGS9GE_ss,R_ss,RK_dark]=steady_state(u_tent,v_tent,v_dark,K1,K2,K3,K4,M,Rec_tot,RK_tot, ...
    n_step_R,mu,lambda,I,E_tot,G_tot,RGS9_tot,k_GE,nu_RG,k_b,k_f,k_cat,tol_fix,tol_stat,k_hyd,k_st,alpha_max,alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    B_ca,F,j_cg_max,f_ca,m_cg,K_cg_max,K_cg_min,K_caM,n_caM,j_ex_sat,K_ex,V_cyto);




    K_cg=K_cg_max + (K_cg_min-K_cg_max)/(1+(v_ss/K_caM)^n_caM);

% Calcolo della corrente dark
    dens_curr_cGMP_in = j_cg_max/Sigma_rod*u_ss.^m_cg./(K_cg^m_cg+u_ss.^m_cg);
    dens_curr_Ca_in   = j_ex_sat/Sigma_rod*v_ss./(K_ex+v_ss);
    curr_in=(dens_curr_cGMP_in+dens_curr_Ca_in)*Sigma_rod;

    


% Valori dei flash


    esp=(5:0.5:8);
   
    
    
 Phi_values=exp(esp)/V_cyto;



T_s=zeros(1,length(esp));

% vettore dei tempi
t=(0:t_step:t_fin);

drop=zeros(length(esp),length(t));



dati=zeros(length(esp),n_step_t+1);

n_flash=length(Phi_values);


E_flash=zeros(n_flash,n_step_t+1);
RGS9GE_flash=zeros(n_flash,n_step_t+1);
R_flash=zeros(n_flash,n_step_t+1);


% ciclo sui flash

for cc=1:length(Phi_values)
    
    disp(cc);

    Phi(1)=Phi_values(cc);



% Integrazione nel tempo





% vettore soluzione
R=zeros(n_step_R,n_step_t+1);
G=zeros(1,n_step_t+1);
G_s=zeros(1,n_step_t+1);
E=zeros(1,n_step_t+1);
GE=zeros(1,n_step_t+1);
RGS9=zeros(1,n_step_t+1);
RGS9GE=zeros(1,n_step_t+1);

u=zeros(1,n_step_t+1);
v=zeros(1,n_step_t+1);

% Condizioni iniziali
R(:,1)=R_ss;
G(1)=G_ss;
G_s(1)=G_s_ss;
E(1)=E_ss;
GE(1)=GE_ss;
RGS9(1)=RGS9_ss;
RGS9GE(1)=RGS9GE_ss;
u(1)=u_ss;
v(1)=v_ss;


% ciclo dei tempi
for cont=2:n_step_t+1
    
    % Iterazione sul valore di u e v
    % Valori di tentativo pari ai valori al passo precedente
    % Valori all'istante precedente
    R_p=R(:,cont-1);
    G_p=G(cont-1);
    G_s_p=G_s(cont-1);
    E_p=E(cont-1);
    GE_p=GE(cont-1);
    RGS9_p=RGS9(cont-1);
    RGS9GE_p=RGS9GE(cont-1);
    u_p=u(cont-1);
    v_p=v(cont-1);

    
    
    R_tent=R_p;
    G_tent=G_p;
    G_s_tent=G_s_p;
    E_tent=E_p;
    GE_tent=GE_p;
    RGS9_tent=RGS9_p;
    RGS9GE_tent=RGS9GE_p;
    u_tent=u_p;
    v_tent=v_p;
    
    % Valore theta dell'intensità di corrente nel tempo
    I_th=I_hys(cont-1)*(1-theta)+I_hys(cont)*theta;
    
    % inizializza l'errore
    err=1;
    
    
    
    % Si aggiunge l'eventuale effetto di un flash aggiuntivo
    R_p(1)=R_p(1)+Phi(cont-1);
    
    % Iterazioni 
    while err>tol_fix
       
    % calcola il valore theta delle incognite
        R_th=R_p*(1-theta)+R_tent*theta;
        G_th=G_p*(1-theta)+G_tent*theta;
        G_s_th=G_s_p*(1-theta)+G_s_tent*theta;
        E_th=E_p*(1-theta)+E_tent*theta;
        GE_th=GE_p*(1-theta)+GE_tent*theta;
        RGS9_th=RGS9_p*(1-theta)+RGS9_tent*theta;
        RGS9GE_th=RGS9GE_p*(1-theta)+RGS9GE_tent*theta;
        u_th=u_p*(1-theta)+u_tent*theta;
        v_th=v_p*(1-theta)+v_tent*theta;
        
        
 
        
        [R_t]=integra_R(R_p,R_th,I_th,v_th,K1,K2,K3,K4,M,Rec_tot,RK_tot,RK_dark,lambda,mu,t_step);
        [G_t,G_s_t,E_t,GE_t,RGS9_t,RGS9GE_t]=integra_GE(G_p,G_s_p,E_p,GE_p,RGS9_p,RGS9GE_p,R_th,G_th,G_s_th,E_th,GE_th,RGS9_th,RGS9GE_th,nu_RG,k_GE,k_cat,k_b,k_f,t_step);

         [u_t,v_t]=integra_uv(u_p,v_p,GE_th,u_th,v_th,k_hyd,E_tot,k_st,alpha_max,alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    j_ex_sat,V_cyto,B_ca,F,K_ex,j_cg_max,f_ca,m_cg,K_cg_max,K_cg_min,K_caM,n_caM,t_step);        
        err=norm([R_t'-R_tent' G_t-G_tent G_s_t-G_s_tent E_t-E_tent GE_t-GE_tent RGS9GE_t-RGS9GE_tent u_t-u_tent v_t-v_tent]); 
        
        
        
        % Aggiorna i valori di tentativo con i valori calcolati per le
        % incognite introducendo il rilassamento
        R_tent=R_t*alpha+R_p*(1-alpha);
        G_tent=G_t*alpha+G_p*(1-alpha);
        G_s_tent=G_s_t*alpha+G_s_p*(1-alpha);
        E_tent=E_t*alpha+E_p*(1-alpha);
        GE_tent=GE_t*alpha+GE_p*(1-alpha);
        RGS9_tent=RGS9_t*alpha+RGS9_p*(1-alpha);
        RGS9GE_tent=RGS9GE_t*alpha+RGS9GE_p*(1-alpha);
        u_tent=u_t*alpha+u_p*(1-alpha);
        v_tent=v_t*alpha+v_p*(1-alpha);
                        
    end
    

    
    % Assegna la soluzione
    R(:,cont)=R_t;
    G(cont)=G_t;
    G_s(cont)=G_s_t;
    E(cont)=E_t;
    GE(cont)=GE_t;
    RGS9(cont)=RGS9_t;
    RGS9GE(cont)=RGS9GE_t;
    u(cont)=u_t;
    v(cont)=v_t;
    
 
end

E_flash(cc,:)=E;
RGS9GE_flash(cc,:)=RGS9GE;
R_flash(cc,:)=R(1,:);

    K_cg=K_cg_max + (K_cg_min-K_cg_max)./(1+(v/K_caM).^n_caM);


% Calcola la corrente
 % la corrente è calcolate in modo esatto
    dens_curr_cGMP = j_cg_max/Sigma_rod*u.^m_cg./(K_cg.^m_cg+u.^m_cg);
    dens_curr_Ca   = j_ex_sat/Sigma_rod*v./(K_ex+v);
    

    curr_tot=(dens_curr_cGMP+dens_curr_Ca)*Sigma_rod;
    drop(cc,:)=1-curr_tot/curr_in;

    
 
% Calcolo di Ts definito tagliando il drop a liv

liv=0.8;
[indici]=find(drop(cc,:)>liv);


ind_inf=min(indici);
ind_sup=max(indici);





T_s(cc)=t(ind_sup)-t(ind_inf);








    %dati(cc,:)=drop;



end

Phi_values=esp;

save('pepp_GCAP','Phi_values','T_s','drop')
