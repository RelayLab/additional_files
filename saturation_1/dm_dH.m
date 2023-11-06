function [dm]=dm_dH(H,H_array,dmdH_array)

% k1 = 0.00001;
% k2 = 1.9;
% k3 = 19000;
% k4 = 1000;
% k5 = 10000;

%% 1
%       if abs(H)<=37
%          m1 = 0;
%       else
%          m1 = - 1.85 * sign(H) / (H^2) / m0;
%       end
%% 2
    %p=[1.66510056949964e-12,7.17355216449709e-25,-1.03526515276425e-07,-3.26948558744223e-20,0.00176226272167064,2.78664895508172e-16,-12.4768167256163,-3.15876357151415e-13,49753.2546574085];
    
%     p=[1.76e-3,0,-12.5,0,4.97e4];           %m0=4*pi*1e-7
%     H_conj = 44;
%     m_conj = polyval(p,H_conj);
%     p1 = polyder(p);
%     dmdH_conj = polyval(p1,H_conj);
%     dmdH_conj1 = -sign(H_conj)*m_conj*H_conj/H_conj^2;
%     c = dmdH_conj/dmdH_conj1;
%     %c = 1;
% 
%     if abs(H)>H_conj
%         dm = -c*sign(H)*m_conj*H_conj/H^2;
%     else
%         dm = polyval(p1,H);       
%     end
    
%% 3
    %m1 = 0;
%% 4
    %m1 = -20*H;
%% 5
%      if abs(H)>10000
%          H=sign(H)*10000;
%      end
     
    for i = 2:length(H_array)
        if H<H_array(i)
            H1 = H_array(i-1);
            H2 = H_array(i);
            dm1 = dmdH_array(i-1);
            dm2 = dmdH_array(i);
            break;
        end
    end
%     [v, i] = min(abs(H_array-H));
%     if (v <= H)
%         H1 = H_array(i);
%         H2 = H_array(i+1);
%         dm1 = dmdH_array(i);
%         dm2 = dmdH_array(i+1);
%      else
%          H1 = H_array(i);
%          H2 = H_array(i-1);
%          dm1 = dmdH_array(i);
%          dm2 = dmdH_array(i-1);
%      end
    
    
    dH = H2-H1;
    ddm = dm2-dm1;
    c = (H - H1)/dH;
    dm = dm1 + ddm*c;

%% 6
    %dm = -20*H;
end