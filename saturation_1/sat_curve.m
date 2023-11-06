%решение диф. уравнений
%function [dm1dH,m3,B3,H]=sat_curve()
function []=sat_curve()

    clear;clc;
    
    m0 = 4*pi*1e-7;
    B0 = 0.6;
    H0 = 10;
    m = B0/H0/m0;
    
    A1 = [
        -500 -1.87 ;
        -400 -1.855;
        -300 -1.84 ;
        -200 -1.825;
        -100 -1.8  ;
        -80  -1.785;
        -60  -1.76 ;
        -40  -1.65 ;
        -30  -1.5  ;
        -20  -1.2  ;
        -10  -0.6  ;
        0   0      ;
        10 	0.6    ;
        20 	1.2    ;
        30 	1.5    ;
        40 	1.65   ;
        60 	1.76   ;
        80 	1.785  ;
        100 1.8    ;
        200 1.825  ;
        300 1.84   ;
        400 1.855  ;
        500 1.87  ];
    
    H1 = A1(:,1);
    B1 = A1(:,2);
    
%% m1 m2
    m1 = B1./H1/m0;
    for i = 1:length(m1)
        if isnan(m1(i))
            m1(i) = m1(i-1);
        end
    end
    
%% dmdH
    H_max = 10000;
    %H = -500:0.01:500;
    H = -H_max:0.1:H_max;
    H_conj = 44;
    
    %p1 = polyfit(H1,m1,8);
    p1=[1.76e-3,0,-12.5,0,4.97e4];
    
    m1_conj = polyval(p1,H_conj);
    dp1 = polyder(p1);
    dm1dH_conj = polyval(dp1,H_conj);
    dm1dH_conj1 = -sign(H_conj)*m1_conj*H_conj/H_conj^2;
    c1 = dm1dH_conj/dm1dH_conj1;
%   
    for i = 1:length(H)    
        if abs(H(i))>H_conj
 %           m4(i) = sign(H(i))*m1_conj*H_conj/H(i)+150;
            m4(i) = sign(H(i))*m1_conj*H_conj/H(i);
        else
            m4(i) = polyval(p1,H(i));                   
        end
    end
    
    for i = 1:length(H)
        %if abs(H(i))>500
        %    dm1dH(i) = -sign(H(i))*m1_conj*H_conj/H(i)^2+150;
        %elseif abs(H(i))>H_conj
        if abs(H(i))>H_conj
            dm1dH(i) = -sign(H(i))*m1_conj*H_conj/H(i)^2;
        else
            dm1dH(i) = 1/c1*polyval(dp1,H(i));                   
        end
    end
    
    %dm1dH(i) = -sign(H(i))*m1_conj*H_conj/H(i)^2;
    %m3 = 1.55e3;
    m3 = 0.17e3;
    for i = 2:length(H)
        m3(i) = m3(i-1) + dm1dH(i)*(H(i)-H(i-1));
    end
    B3 = m3*m0.*H;
     
    for i = 2:length(H)
        dm4dH(i) = (m4(i)-m4(i-1))/(H(i)-H(i-1));
    end
    
%% plot
    subplot (3,1,1)
    plot(H1,B1,'.r',...
         H ,B3,'b'  ,...
         ...%H ,m4.*m0.*H,'g',...
        'MarkerSize',15,...
        'LineWidth',1.3);
    
    grid on
    grid minor
    title ('B(H)');
    legend ('GOST', 'restored');
    legend ('show');
    %xlim([-H_max H_max]);
    %ylim([-2.5 2.5]);
    xlim([0 500]);
    ylim([0 2]);
    ylabel('B, Ts');
    xlabel('H, A*m');
    
%%    
    subplot (3,1,2)
    plot(H1,m1,'r',...
         H ,m3,'b',...,
         ...%H ,m4,'g',...,
         'MarkerSize',15,...
         'LineWidth',1.3);
    grid on
    grid minor
    title ('m(H)');
    legend ('GOST', 'restored');
    legend ('show');
    ylabel('m, ');
    xlabel('H, A*m');
    xlim([-H_max H_max]);
    ylim([0 6e4]);
    
    
%%

    subplot (3,1,3)
    plot(H,dm1dH,'b',...
         H,dm4dH,'g--',...
        'LineWidth',1.3);
    grid on
    title ('dm/dt');
    grid on
    grid minor
    legend ('GOST, restored');
    legend ('show');
    ylabel('dm/dt, ');
    xlabel('H, A*m');
    xlim([-H_max H_max]);
    %xlim([-44.5 -43.5]);
    %ylim([-1000 1000]);
    
end