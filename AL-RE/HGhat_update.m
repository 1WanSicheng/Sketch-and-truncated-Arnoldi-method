function Hdp1hat = HGhat_update(A,U,hdp1d,Hdhat,THdp1,Tudp1,tow_d,d,r)
    Ed = [zeros((d-1)*r,r);eye(r)];
    
    
    %get hid 
    Udp1 = U(:,d*r+1:(d+1)*r);
    Utildenew = A*Udp1;
    H=[];
    for i = 1:d
        Ui =U(:,(i-1)*r+1:i*r);
        hidp1 = Ui'*Utildenew;
        H=[H,hidp1'];
    end
        
    
    tow_dp1 = Tudp1(d*r+1:(d+1)*r, d*r+1:(d+1)*r);
    hdp1dp1 = Udp1'*Utildenew;
    Tud = Tudp1(1:d*r,1:d*r);
    
    %LT means the left top enties in Hd+1 hat
    %LD,RD are the left down and right down, respectively
    LT = Hdhat + THdp1*hdp1d/tow_d*Ed';
   
   
    Hhatnew = -LT*THdp1/tow_dp1+ Tud*H'/tow_dp1+THdp1*hdp1dp1/tow_dp1;
    LD = tow_dp1*hdp1d/tow_d*Ed';
    RD = tow_dp1*(-hdp1d/tow_d*Ed'*THdp1+hdp1dp1)/tow_dp1;
    
    Hdp1hat = [LT,Hhatnew;
        LD,RD];
    
    
    
end