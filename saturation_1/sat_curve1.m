%—Ä–µ—à–µ–Ω–∏–µ –¥–∏—Ñ. —É—Ä–∞–≤–Ω–µ–Ω–∏–π
function [dm3dH,m3,B3,H]=sat_curve1()
%function []=sat_curve1()

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
    H_max = 40000;
    %H_max = 200;
    H_sat = 500;
    H = -H_max:0.1:H_max;
    H_conj = 70;
    
    p1 = polyfit(H1,B1,17);
         %1 2 3 4 5 6  7         8   9       10   11      12    13    14   15      16   17    18
    %p1 = [0,0,0,0,0,0,-1.343e-21,0,4.826e-17,0,-7.470e-13,0,5.658e-09,0,-2.312e-05,0,6.672e-2,0];
    %p1=[1.76e-3,0,-12.5,0,4.97e4];
    
%     m1_conj = polyval(p1,H_conj);
%     dp1 = polyder(p1);
%     dm1dH_conj = polyval(dp1,H_conj);
%     dm1dH_conj1 = -sign(H_conj)*m1_conj*H_conj/H_conj^2;
%     c1 = dm1dH_conj/dm1dH_conj1;
      c = 0.03;
      B1_conj = polyval(p1,H_conj);
      B2_conj = c*log(H_conj);
%   

    for i = 1:length(H)    
        %if abs(H(i))>H_sat
        %    B3(i) = 1;
        %elseif abs(H(i))>H_conj
        if abs(H(i))>H_conj
            %m3(i) = sign(H(i))*m1_conj*H_conj/H(i);
            B3(i) = sign(H(i))*(c*log(abs(H(i)))+(B1_conj-B2_conj));
        else
            B3(i) = polyval(p1,H(i));                   
        end
    end
    
    m3 = B3./H/m0;
    for i=2:length(H)
        if isnan(m3(i))||isinf(m3(i))
            m3(i) = m3(i-1);
        end
    end
    for i = 2:length(H)
        dm3dH(i) = (m3(i)-m3(i-1))/(H(i)-H(i-1));
    end
    
    for i = 2:length(H1)
        dm1dH(i) = (m1(i)-m1(i-1))/(H1(i)-H1(i-1));
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
    legend ('√Œ—“', 'ÙÓÏÛÎ‡');
    legend ('show');
    xlim([-H_max H_max]);
    %xlim([-80 -60]);
    ylim([-2.5 2.5]);
    %xlim([-500 500]);
    %ylim([-2 2]);
    ylabel('B, “Î');
    xlabel('H, A/Ï');
    
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
    legend ('√Œ—“', 'ÙÓÏÛÎ‡');
    legend ('show');
    ylabel('m, √Ì/Ï');
    xlabel('H, A/Ï');
    xlim([-H_max H_max]);
    %xlim([-80 -60]);
    ylim([0 6e4]);
    
%%
    subplot (3,1,3)
    plot(H,dm3dH,'b',...
        'LineWidth',1.3);
    grid on
    title ('dm/dt');
    grid on
    grid minor
    legend ('ÙÓÏÛÎ‡');
    legend ('show');
    ylabel('dm/dt, √Ì/(Ï*Ò)');
    xlabel('H, ¿/Ï');
    %xlim([-80 -60]);
    xlim([-H_max H_max]);
end