function [X1,X2]= STArnoldi_AL2(A, B, C1, C2, SU, SV, k, maxit, tol, p)
%Truncated Arnoldi method for Sylvester equations
    %Initialization
    %Compute skinny QRs:
    n1 =size(C1,1);
    n2 =size(C2,1);
    r = size(C1,2);
    [U1,~] = qr(C1,0);
    [V1,~] = qr(C2,0);
    [Qu1,beta1] = qr(SU*C1,0);
    [Qv1, beta2] = qr(SV*C2,0);
    Tu1 = beta1;
    Tv1 = beta2;
    U=[U1];
    V=[V1];
    Hd = [];
    Gd= [];
    %intialize Hd hat and Gd hat
   %Hdhat = Tu1*U1'*A*U1/Tu1;
   %Gdhat = Tv1*V1'*B'*V1/Tv1;
    for d=1:maxit
        Ud = U(1:n1, (d-1)*r+1:d*r);
        Vd = V(1:n2, (d-1)*r+1: d*r);
        
        Utilde = A*Ud;
        Vtilde = B'*Vd;
        if 1<d-k+1
            max = d-k+1;
        else
            max = 1;
        end
        
        for i=max:d
            Ui = U(:,(i-1)*r+1:i*r);
            Vi = V(:, (i-1)*r+1:i*r);
            hid = Ui'*Utilde;
            gid = Vi'*Vtilde;
           
            Utilde = Utilde - Ui*hid;
            Vtilde = Vtilde - Vi*gid;
        end
        
        %Compute skinny QRs of Utilde and Vtilde
        [Udp1,hdp1d] = qr(Utilde,0);
        [Vdp1,gdp1d] = qr(Vtilde,0);
        U = [U,Udp1];
        V = [V,Vdp1];
        


        
        %Update Hd and Gd
        Hd = HGd_AL2(A,U,d,r,n1);
        Gd = HGd_AL2(B',V,d,r,n2);
        
         
         Ed = [zeros((d-1)*r,r);eye(r)];
         E1 = [eye(r);zeros((d-1)*r,r)];
         
         
         if mod(d,p) == 0 
             %solve lyap equation lyap(A,B,C) solves AX+XB+C = 0
             C = -E1*beta1*beta2'*E1';
             Y = lyap(Hd,Gd',C);
             
             rho = sqrt(norm(hdp1d*Ed' * Y, 'fro')^2 + norm(Y *Ed*gdp1d', 'fro')^2);
             disp(rho);
             if rho < tol
                 break;
             end
         end

         
    end
    
   Rank = rank(Y);
   [Y1, Y2] = LowRankApproximation(Y, Rank);


   X1 = (U(:,1:d*r)/Tud) * Y1;
   X2 = (V(:,1:d*r)/Tvd) * Y2;

end