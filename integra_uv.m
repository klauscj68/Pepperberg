% integrazione delle equazioni per R*
function [u_t,v_t]=integra_uv(u_p,v_p,E_th,u_th,v_th,k_hyd,E_tot,k_st,alpha_max,alpha_min, A1, A2, KK1, KK2, m1, m2, ...
    j_ex_sat,V_cyto,B_ca,F,K_ex,j_cg_max,f_ca,m_cg,K_cg_max,K_cg_min,K_caM,n_caM,t_step)








    % Calcolo di u e v col metodo theta
    
    
    %u_t=u_p+(-k_hyd*(E_tot-E_th)/2*u_th-k_st*E_th*u_th +(alpha_min +(alpha_max-alpha_min)*k_cyc^m_cyc/(k_cyc^m_cyc+v_th^m_cyc)))*t_step;

    u_t=u_p+(-k_hyd*(E_tot-E_th)/2*u_th-k_st*E_th*u_th +(alpha_min+(alpha_max-alpha_min)*(A1/(1+(v_th/KK1)^m1)+A2/(1+(v_th/KK2)^m2))))*t_step;

    K_cg=K_cg_max + (K_cg_min-K_cg_max)/(1+(v_th/K_caM)^n_caM);

    v_t=v_p+(-j_ex_sat/(V_cyto*B_ca*F)*v_th/(K_ex+v_th)+j_cg_max*f_ca/(V_cyto*B_ca*F*2)*u_th^m_cg/(K_cg^m_cg+u_th^m_cg))*t_step;
    

 
return
