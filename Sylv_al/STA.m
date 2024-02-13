function [X1, X2] = STA(A, B, C1, C2, SU, SV, k, maxit, tol, p)
    % Initialize variables
    r = size(C1,2);
    [U1,~]=qr(C1,0);
    [V1,~]=qr(C2,0);
    [QU1,beta1] = qr(SU*C1,0);
    [QV1,beta2] = qr(SU*C2,0);
    
    Tu1 = beta1;
    Tv1 = beta2;
    
    U_1 = U1;
    V_1 = V1;
    H1 = U1'*A*U1;
    G1 = vV1'*B'*V1;
    Hhat_0 = Tu1*H1*inv(Tu1);
    Ghat_0 = Tv1*G1*inv(Tv1);
    
    
    for d= 1:maxit
        Utilde = A*Ud ;
        Vtilde = B'*Vd;
        max = max{1,d-k+1};
        Ed = [zeros((d-1)*r,r),eye(r)];
        

        for i = max:d 
            Utilde = Utilde - Ui*Ui'*Utilde;
            Vtilde = Vtilde - Vi*Vi'*Vtilde;
        end
        
        j = d+1;
        h_(d,d-1) = hjd;
        [Uj, hjd]= qr(Utilde,0);
        [Vj,gjd] = qr(Vtilde,0);
        
        
        %Update Qu,d+1 Tu,d+1, which are Sx(d+1)r , and (d+1)r x (d+1)r dimension
        U_j = [U_d, Uj];
        V_d = [V_d, Vj];
        [Quj, Tuj] = qr(SU*U_j,0);
        [Qvj,Tvj]= qr(SV*V_j,0);
        
        
        %Update Hd ,H
        tow_d = Tud((d-1)*r:d*r, (d-1)*r:d*r);
        theta_d = Tvd((d-1)*r:dr, (d-1)*r:dr);

        Hhat = Tuj(1:dr,dr:jr)*hjd*inv(tow_d);
        Ghat = Tvj(1:dr,dr:jr)*gjd*inv(theta_d);
        
        Hh = [];
        Hg=[];
        for k = 1:d
            Hh = [Hh, Ui'*Utilde];
            Hg = [Hg, Vi'*Vtilde];
        end
        
        %update Hhat_d+1

        topleft = Hhat_(d-1) + Tud(1:(d-1)*r,(d-1)*r:dr)*h_(d,d-1)*inv(tow_(d-1))*E(d-1)';
        downleft = tow_(d)*h_(d,d-1)*inv(tow_(d-1))*E(d-1)';
        
        Hhat_new = topleft*Tud(1:(d-1)*r,(d-1)*r:dr)*inv(tow_d)
        +Tud(1:(d-1)*r,1:(d-1)*r)*Hh*inv(tow_d)+ Tud(1:(d-1)*r,(d-1)*r:dr)*inv(tow_(d-1));
        
        
        
        
        
        
    end
        
        
        
        
       
            
        
    
    
    




end
