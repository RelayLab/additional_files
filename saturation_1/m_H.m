function [m]=m_H(H,H_array,m_array)

% k1 = 0.00001;
% k2 = 1.9;
% k3 = 19000;
% k4 = 1000;
% k5 = 10000;
 m0 = (4*pi*1e-7);

%% 1
%     if abs(H)<37
%         m = 0.05/m0;
%     else
%         m = 1.85.*sign(H)./H/m0;
%     end
    
%%  2
    %p=[1.66510056949964e-12,7.17355216449709e-25,-1.03526515276425e-07,-3.26948558744223e-20,0.00176226272167064,2.78664895508172e-16,-12.4768167256163,-3.15876357151415e-13,49753.2546574085];           %m0=4*pi*1e-7
    
    
%     p=[1.76e-3,0,-12.5,0,4.97e4];
%     
% 
%     H_conj = 44;
%     m_conj = polyval(p,H_conj);
%  
%     if abs(H)>H_conj
%         m = m_conj*H_conj/abs(H)+140;
%     else
%         m = polyval(p,H);
%     end

%% 3
    %m = 49700;
%% 4
    %m = -10*(H^2)+5e4;
    
%% 5

%      if abs(H)>10000
%          H=sign(H)*10000;
%      end
    for i = 2:length(H_array)
        if H<H_array(i)
            H1 = H_array(i-1);
            H2 = H_array(i);
            m1 = m_array(i-1);
            m2 = m_array(i);
            break;
        end
    end
    
%     [v, i] = min(abs(H_array-H));
%     q = H_array(i);
%     if (abs(H_array(i)) <= abs(H))    
%         H1 = H_array(i);
%         H2 = H_array(i+1);
%         m1 = m_array(i);
%         m2 = m_array(i+1);
%     else
%         H1 = H_array(i-1);
%         H2 = H_array(i);
%         m1 = m_array(i-1);
%         m2 = m_array(i);
%     end
    
    dH = H2-H1;
    dm = m2-m1;
    c = (H - H1)/dH;
    m = m1 + dm*c;
%% 6
    %m = -10*H^2+5e4;

end