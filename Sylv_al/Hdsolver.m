function L= Hdsolver(A,B,C,D,d,r)
   % A is Hhat_d-1; B is Tud, C is Utilde , and D is hd,d-1 here.
    
    Edmin1 = [zeros((d-2)*r,r),eye(r)];
    H = [];
    
    for i = 1:d-1
        hid = C(:,r*(i-1)+1:r*i)'*C;
        H = [H,hid'];
    end
    hdd = C(:,r*(d-1)+1:r*d)'*C;
    
    % get Tu,d-1 here
    Tudmin1 = B(1:(d-1)*r,1:(d-1)*r);
    THd = B(1:(d-1)*r, (d-1)*r+1:d*r);
    towd = B((d-1)*r+1:d*r,(d-1)*r+1:d*r);
    towdmin1 = Tudmin1((d-2)*r+1:(d-1)*r, (d-2)*r+1:(d-1)*r);
    
    lefttop = A + THd*D*towdmin1*Edmin1';
    righttop = lefttop*THd*inv(towd) + Tudmin1*H'*inv(towd)+ THd*hdd*inv(towdmin1);
    leftdown = towd*D*inv(towdmin1)*Edmin1';
    rightdown = leftdown*THd*inv(towd)+hdd*inv(towd);
    
    
    L = [lefttop,-righttop;
        leftdown,-rightdown];
    
  
end
