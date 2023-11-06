%СЂРµС€РµРЅРёРµ РґРёС„. СѓСЂР°РІРЅРµРЅРёР№
function [dm3dH,m3,B3,H]=sat_curve2()
%function []=sat_curve1()

    clear;clc;
    
    m0 = 4*pi*1e-7;
    B0 = 0.6;
    H0 = 10;
    m = B0/H0/m0;
    
    A1 = [
        -14  -12.1;
        -13  -12  ;
        -12  -11.9;
        -11  -11.8;
        -10  -11.7;
        -9   -11.5;    
        -8   -11.2;
        -7   -11  ;
        -6   -10.5;
        -5   -10.2;
        -4   -9.7 ;
        -3   -9   ;
        -2   -8   ;
        -1   -5.3 ;
        0   0   ;
        1   5.3 ;
        2   8   ;
        3   9   ;
        4   9.7 ;
        5   10.2;
        6   10.5;
        7   11  ;
        8   11.2;
        9   11.5;
        10  11.7;
        11  11.8;
        12  11.9;
        13  12  ;
        14  12.1];
    
    A1 = [10*A1(:,1) 0.1*A1(:,2)];
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
    %H_max = 500;
    H_sat = 500;
    H = -H_max:0.1:H_max;
    H_conj = 35;
    
    p1 = polyfit(H1,B1,19);
         %1 2 3 4 5 6  7         8   9       10   11      12    13    14   15      16   17    18
    %p1 = [0,0,0,0,0,0,-1.343e-21,0,4.826e-17,0,-7.470e-13,0,5.658e-09,0,-2.312e-05,0,6.672e-2,0];
    %p1=[1.76e-3,0,-12.5,0,4.97e4];
    
%     m1_conj = polyval(p1,H_conj);
%     dp1 = polyder(p1);
%     dm1dH_conj = polyval(dp1,H_conj);
%     dm1dH_conj1 = -sign(H_conj)*m1_conj*H_conj/H_conj^2;
%     c1 = dm1dH_conj/dm1dH_conj1;
      c = 0.1;
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
    
    B3 = B3 + m0*H;
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
    legend ('опыт', 'формула');
    legend ('show');
    xlim([-H_max H_max]);
    %xlim([-80 -60]);
    %ylim([-1.5 1.5]);
    %xlim([-500 500]);
    %ylim([-2 2]);
    ylabel('B, Тл');
    xlabel('H, A/м');
    
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
    legend ('опыт', 'формула');
    legend ('show');
    ylabel('m, Гн/м');
    xlabel('H, A/м');
    xlim([-H_max H_max]);
    %xlim([-80 -60]);
    %ylim([0 6e4]);
    
%%
    subplot (3,1,3)
    plot(H,dm3dH,'b',...
        'LineWidth',1.3);
    grid on
    title ('dm/dH');
    grid on
    grid minor
    legend ('формула');
    legend ('show');
    ylabel('dm/dH, Гн/А');
    xlabel('H, А/м');
    %xlim([-80 -60]);
    xlim([-H_max H_max]);
end