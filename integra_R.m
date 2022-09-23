% integrazione delle equazioni per R*
function [R_t]=integra_R(R_p,R_th,I_th,v_th,K1,K2,K3,K4,M,Rec_tot,RK_tot,RK_dark,lambda,mu,t_step)




% Il numero di rodopsine in ciascun stato di fosforilazione sono date dalla
% risoluzione di un sistema di equazioni differenziali ordinarie. I
% coefficienti di questo sistema sono funzioni del calcio.








% ritiene in lambd_max il valore a calcio dark delle lambda
lambda_max=lambda;
Ca=v_th;
        


% I tassi ti fosforilazione lambda dipendono dal calcio, i tassi di legame con l'arrestina sono indipendenti dal calcio.    
C1=(Ca/K1)^2*(1/K3+1/K4*M/K2)*Rec_tot;
C2=1+(Ca/K1)^2*(1+M/K2);


% Calcolo di Rec/Rec_tot
aa=C1*C2;
bb=(C1*(RK_tot/Rec_tot-1)+C2);
Rec_Rectot=(-bb+sqrt(bb^2+4*aa))/(2*aa);
RK=RK_tot*(1+C1*Rec_Rectot)^(-1);
lambda=lambda_max*RK/RK_dark;


% % Si esclude il feedback del calcio sulla rodopsina
% lambda=lambda_max;
    % Calcolo delle R*
    % Matrice dei coefficienti del sistema ode
    M=diag(-(lambda+mu),0)+diag(lambda(1:end-1),-1);
    
    % Calcolo del valore di 
    
    R_t=R_p+M*R_th*t_step;
    
    % Aggiunge il pezzo dovuto all'intensità di corrente
    R_t(1)=R_t(1)+t_step*I_th;
    
    
    
 
return
