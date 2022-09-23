% integrazione delle equazioni per R*
function   [G_t,G_s_t,E_t,GE_t,RGS9_t,RGS9GE_t]=integra_GE(G_p,G_s_p,E_p,GE_p,RGS9_p,RGS9GE_p,R_th,G_th,G_s_th,E_th,GE_th,RGS9_th,RGS9GE_th,nu_RG,k_GE,k_cat,k_b,k_f,t_step)








    % Calcolo di G* e E* col metodo theta
    

    G_t=G_p+t_step*(-G_th*nu_RG*R_th+k_cat*RGS9GE_th);

    
    G_s_t=G_s_p+t_step*(G_th*nu_RG*R_th-k_GE*E_th*G_s_th);

        
    E_t=E_p+t_step*(-k_GE*E_th*G_s_th+k_cat*RGS9GE_th);

 
    GE_t=GE_p+t_step*(k_GE*E_th*G_s_th-k_f*RGS9_th*GE_th+k_b*RGS9GE_th);

    
    RGS9_t=RGS9_p+t_step*(k_cat*RGS9GE_th+k_b*RGS9GE_th-k_f*RGS9_th*GE_th);


    RGS9GE_t=RGS9GE_p+t_step*(k_f*RGS9_th*GE_th-k_b*RGS9GE_th-k_cat*RGS9GE_th);

    
    
    
 
    
    
 
return
