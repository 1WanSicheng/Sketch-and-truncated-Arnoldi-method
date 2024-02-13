function Hdhat = HGhat_computation(A,U,d,r,n1,Tud)
H= zeros(d*r);
Uisetd = zeros(d*r,n1);

for i = 1:d 
    Ui = U(:,(i-1)*r+1:i*r);
    Uisetd((i-1)*r+1:i*r,:)= Ui';
end

for j =1:d
    Uj = U(:,(j-1)*r+1:j*r);
    Utilde_j = A*Uj;
    H(:,(j-1)*r+1:j*r) = Uisetd*Utilde_j;
end

Hdhat  =Tud*H/Tud;
end

