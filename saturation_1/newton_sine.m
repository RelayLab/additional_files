
function []=newton()

    clear;clc;
    clf;
    warning('off','MATLAB:nearlySingularMatrix')

    [dmdH_graph,m_graph,B_graph,H_graph]=sat_curve2;
    %B_graph = 0; H_graph=0; m_graph = 0; dmdH_graph = 0;
    
    %исходные параметры
    h       = 0.25e-4;           
    t_max   = 0.06;          
    t       = 0:h:t_max;
    n_max   = length(t);
    i_max   = 20;
    epsilon = 1e-3;
    a       = 0.001;
    
    %%%%%%%%%%%%%%%%%%%%ток i1%%%%%%%%%%%%%%%%%%%%%%%%
    %Ta  = 0.128;
    Ta  = 0.02;
    Imax= 25*500*sqrt(2);
    %Imax= 23.145e3*sqrt(2);
    i1_array = zeros(size(t,1));
    t_start = 0.00;
    for i = 1:n_max
        t_elapsed = t(i)-t_start;
        if t(i)<t_start
            i1_array(i) = 0;
        else
            i1_array(i) = -Imax*(   cos(314.15*t(i))-0*exp(-t_elapsed /Ta)   );
        end
    end
    
    %%%%%%%%%%%%%%%%%%%константы%%%%%%%%%%%%%%%%%%%%%%%%%
    R2 = 0.6;
    S = 20e-4;
    len = 0.445;
    w1 = 1;
    w2 = 100;
%     R2  = 7.51+6.3;
%     S   = 13.125*1e-4;
%     len = 1.437;
%     w1  = 1;
%     w2  = 1997;
    
    m0  = 4*pi*1e-7;
    m   = m_H(0,H_graph,m_graph);
    L1m = m*m0*w1^2*S/len;
    L2m = m*m0*w2^2*S/len;
    L1s = L1m/m;
    L2s = L2m/m;
    Ln  = 0;
    L1  = L1m+L1s;
    L2  = L2m+L2s;
    M   = w1*w2*m*m0*S/len;
   
    %%%%%%%%%%%%%%%%%%%%нулевые исходные данные%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p1     =0;
    p1s    =0;
    p1m    =0;
    p2     =0;
    p2m    =0;
    e2     =0;
    i2     =0;
    p0     =0;
    H0     =0;
    dp1m_dt=0;
    dp2_dt =0;    
    p1m_0  =0;
    p2_0   =0;
    B0     =zeros(1,n_max);

    %1    2    3    4   5    6     7   8     9   10      11      12    13    14
    %p1  p1s  p1m   p2  p2m  e2    i2  p0    H0  dp1m_dt dp2_dt  L1m   L2m   m
    X = [
        p1     ;
        p1s    ;
        p1m    ;
        p2     ;
        p2m    ;
        e2     ;
        i2     ;
        p0     ;
        H0     ;
        dp1m_dt;
        dp2_dt ;
        L1m    ;
        L2m    ;
        m      ];    
    
      
      %p1    p1s    p1m   e2         p2m  dp1m_dt  i2    di2_dt  p0     H0    L1m        L2m     m
    %[2.048 2.047 2.048  3303371.009 0.9  357.08  7.935 14487   4.096  49855 5.7176e-05  228    49815]
    err = epsilon*[
        1;%p1  - (L1m+L1s) * i1                ;
        1;%p1  - p1s - p1m + p2m               ;
        1e6;%e2  + w2*dp1m_dt                    ;
        1e6;%-e2  - dp2_dt - i2*R2               ;
        1;%p2 - (L2m+L2s) * i2                 ;
        1;%p1m - L1m * i1                      ;
        1;%p2m - L2m/w2 * i2                   ;
        1;%p1m - p1m_0 - h*dp1m_dt             ;
        1;%p2  - p2_0  - h*dp2_dt              ;
        1e-3;%p0 - p1m + p2m                      ;
        1e-3;%p0 - H0*S*m*m0                      ;
        1e-5;%L1m - m*m0*w1*S/len                 ;
        1e2;%L2m - m*m0*w2^2*S/len               ;
        1e4];%m - m_H(H0,H_graph,m_graph)         ]

    
    
    
    ii = 0;
    %далее делаем n_max шагов метода Эйлера
    %на каждом шаге решаем систему диф. ур, предварительно записав
    for n = 2:n_max
        i1 = i1_array(n-1);
        for i = 1:i_max

            %1    2    3    4   5    6     7   8     9   10      11       12    13    14
            %p1  p1s  p1m   p2  p2m  e2    i2  p0    H0  dp1m_dt dp2_dt  L1m   L2m   m
            W = [
                p1  - (L1m+L1s) * i1                ;
                p1  - p1s - p1m + p2m               ;
                e2  + w2*dp1m_dt                    ;
                -e2  - dp2_dt - i2*R2               ;
                p2 - (L2m+L2s) * i2                 ;
                p1m - L1m/w1 * i1                   ;
                p2m - L2m/w2 * i2                   ;
                p1m - p1m_0 - h*dp1m_dt             ;
                p2  - p2_0  - h*dp2_dt              ;
                p0 - p1m + p2m                      ;
                p0 - H0*S*m*m0                      ;
                L1m - m*m0*w1^2*S/len               ;
                L2m - m*m0*w2^2*S/len               ;
                m - m_H(H0,H_graph,m_graph)         ];

            dmdH            = dm_dH(H0,H_graph,dmdH_graph);
            dmdH_array(n)   = dmdH;
            
            %Матрица Якоби
                %1    2    3    4   5       6       7    8        9    10       11    12    13    14
                %p1  p1s  p1m   p2  p2m     e2      i2   p0       H0  dp1m_dt dp2_dt  L1m   L2m   m
            dW_dt = [ 
                1     0    0    0    0      0       0    0        0     0        0    -i1   0     0          ;%p1  - (L1m+L1s) * i1
                1    -1   -1    0    1      0       0    0        0     0        0     0    0     0          ;%p1  - p1s - p1m + p2m 
                0     0    0    0    0      1       0    0        0     w2       0     0    0     0          ;%e2  + w2*dp1m_dt 
                0     0    0    0    0     -1      -R2   0        0     0       -1     0    0     0          ;%e2  - dp2_dt - i2*R2 
                0     0    0    1    0      0 -(L2m+L2s) 0        0     0        0     0   -i2    0          ;%p2 - (L2m+L2s) * i2
                0     0    1    0    0      0       0    0        0     0        0  -i1/w1  0     0          ;%p1m - L1m/w1 * i1    
                0     0    0    0    1      0   -L2m/w2  0        0     0        0     0 -i2/w2   0          ;%p2m - L2m/w2 * i2   
                0     0    1    0    0      0       0    0        0    -h        0     0    0     0          ;%p1m - p1m_0 - h*dp1m_dt
                0     0    0    1    0      0       0    0        0     0       -h     0    0     0          ;%p2  - p2_0  - h*dp2_dt
                0     0   -1    0    1      0       0    1        0     0        0     0    0     0          ;%p0 - p1m + p2m      
                0     0    0    0    0      0       0    1     -S*m*m0  0        0     0    0 -H0*S*m0       ;%p0 - H0*S*m*m0
                0     0    0    0    0      0       0    0        0     0        0     1    0 -m0*w1^2*S/len ;%L1m - m*m0*w1^2*S/len   
                0     0    0    0    0      0       0    0        0     0        0     0    1 -m0*w2^2*S/len ;%L2m - m*m0*w2^2*S/len     
                0     0    0    0    0      0       0    0     -dmdH    0        0     0    0     1         ];%m - m_H(H0,m0) 
             
            if i == 1
                X(:,n) = X(:,n-1) - dW_dt\W;
            else
                X(:,n) = X(:,n) - dW_dt\W;               
            end
            
            %1    2    3    4   5       6       7    8        9    10       11    12    13    14
            %p1  p1s  p1m   p2  p2m     e2      i2   p0       H0  dp1m_dt dp2_dt  L1m   L2m   m
            p1      = X(1,n);
            p1s     = X(2,n);
            p1m     = X(3,n);
            p2      = X(4,n);
            p2m     = X(5,n);
            e2      = X(6,n);
            i2      = X(7,n);
            p0      = X(8,n);
            H0      = X(9,n);
            dp1m_dt = X(10,n);
            dp2_dt  = X(11,n);
            L1m     = X(12,n);
            L2m     = X(13,n);
            m       = X(14,n);
            B0(1,n) = H0*m*m0;

            %ещё раз вычисляем W, чтобы определить невязку после итерации
            W = [
                p1  - (L1m+L1s) * i1                ;
                p1  - p1s - p1m + p2m               ;
                e2  + w2*dp1m_dt                    ;
                -e2  - dp2_dt - i2*R2               ;
                p2 - (L2m+L2s) * i2                 ;
                p1m - L1m/w1 * i1                   ;
                p2m - L2m/w2 * i2                   ;
                p1m - p1m_0 - h*dp1m_dt             ;
                p2  - p2_0  - h*dp2_dt              ;
                p0 - p1m + p2m                      ;
                p0 - H0*S*m*m0                      ;
                L1m - m*m0*w1^2*S/len               ;
                L2m - m*m0*w2^2*S/len               ;
                m - m_H(H0,H_graph,m_graph)         ];
        
             %*****НАЧАЛО***********проверка на окончание итераций******************
            has_solved = true;
            for k = 1:size(X,1)
                if abs(W(k))>abs(err(k))
                    if i == i_max
                        %continue
                        error('Step n=%d, equation=%d, t=%d. Max iteration is reached.',n,k,n*h);
                    end
                    has_solved = false;
                    break;
                end
            end
            if has_solved == true
                if (i>ii)
                    ii = i;
                end
                break
            end
            %*****КОНЕЦ***********проверка на окончание итераций******************
        
            
        end
        
        
        p2_0  = p2;
        p1m_0 = p1m;
    end
    
    %max(iii)
    
    x_min = 0;
    x_max = t_max;
    no = 6;
    noi = 1;
    
    subplot(no,1,noi);
    plot(t,i1_array/1000,'r',...
        'LineWidth',1.3);
    title('i_1(t)');
    ylabel('i_1, кА');
    grid;
    xlim([x_min x_max]);
    %ylim([-30 30]);
    %
    noi = noi + 1;
    subplot(no,1,noi);
    plot(t,X(7,:),'g',...
        ...%t,i1_array*w1/w2,'r'
         'LineWidth',1.3);
    title('i_2(t)');
    ylabel('i_2, А');
    grid;
    xlim([x_min x_max]);
    %ylim([-300 300]);
    
    noi = noi + 1;
    subplot(no,1,noi);
    plot(t,X(14,:),'b',...
        'LineWidth',1.3);
    title('{\mu}(t)');
    ylabel('{\mu}, Гн/м');
    grid;
    ylim([0 5e4]);
    
    noi = noi + 1;
    subplot(no,1,noi);
    plot(t,B0,'b',...
        'LineWidth',1.3);
    title('B_0(t)');
    ylabel('B_0, Тл')
    grid;
    ylim([-2 2]);
    
    noi = noi + 1;
    subplot(no,1,noi);
    plot(t,X(9,:)/1000,'b',...
        'LineWidth',1.3);
    title('H_0(t)');
    grid;
    ylim([-20 20]);
    ylabel('H_0, кА/м')
    
    %1e6 1e-5  200
    noi = noi + 1;
    subplot(no,1,noi);
    plot(t,i1_array*w1/w2-X(7,:),'m',...
        'LineWidth',1.3);
        
    title('i_0(t)');
    %legend ('опыт', 'расчёт');
    grid;
    ylabel('i_0, А')
    xlabel('t, с')
    ylim([-100 100]);
    disp('done')
    ii
warning('on','MATLAB:nearlySingularMatrix')
end
