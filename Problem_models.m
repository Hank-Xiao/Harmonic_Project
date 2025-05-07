
function Fit=Problem_models(x,Number)

switch Number

    case 3 %拉伸/压缩弹簧设计(Tension/compression spring design）            
        % x=0.447*cos(19.75*2*pi*t+x(1))+x(2)*cos(2*19.75*2*pi*t+x(3))+x(4)*cos(3*19.75*2*pi*t+x(5))+x(6)*cos(4*19.75*2*pi*t+x(7));
        % y=0.5*cos(19.75*2*pi*t+pi/9)+1*cos(2*19.75*2*pi*t+pi/6)+2*cos(3*19.75*2*pi*t+pi/3)+3*cos(4*19.75*2*pi*t+pi/2);
        % x=99.96*cos(50.1*2*pi*t+x(1))+x(2)*cos(250*2*pi*t+x(3))+x(4)*cos(350*2*pi*t+x(5))+x(6)*cos(650*2*pi*t+x(7))+x(8)*cos(x(9)*2*pi*t+x(10));%原信号
        % w0 = 100;
        % y=380*cos(w0*pi*t+pi/18)+20*cos(3*w0*pi*t+pi/4)+16.4*cos(5*w0*pi*t+5*pi/9)+...
        %      +12*cos(7*w0*pi*t+5*pi/6)+10*cos(9*w0*pi*t+7*pi/6);
        % %+3.2*cos(0.5*w0*pi*t+pi/9)+2.8*cos(3.6*w0*pi*t+pi*5/36);
        %      % 1.8*cos(6.4*w0+pi*t+2*pi/3)+1.4*cos(7.6*w0*pi*t+pi);
        % x=380.01*cos(100*pi*t+x(1))+19.9999*cos(3*100*pi*t+x(2))+16.3988*cos(5*100*pi*t+x(3))+...
        %         +11.9888*cos(7*100*pi*t+x(4))+10.0002*cos(9*100*pi*t+x(5));
        % f = [25 50 150 180 250 320 350 380 450]; % Frequency of signal components (Hz)
        % A = [3.2 380 20 2.8 16.4 1.8 12 1.4 10]; % Amplitude of signal components
        % ph = deg2rad([20 60 30 25 100 120 150 180 210]); % Phase of signal components (rad)
        
        f = [50 150 250 350 450 550]; % Frequency of signal components (Hz)
        A = [220 6.3 5.5 3.2 1.5 0.9]; % Amplitude of signal components
        ph = deg2rad([30 44.5 60 20 88.3 90]); % Phase of signal components (rad)
        
        Fs = 1024; % Sampling freq

        L = 4095; % Length of signal = 2*N-1
        
        N = (L+1)/2;
        t=(-N:N)/Fs;
        y = zeros(1,length(t));
        
        for k = 1:length(f)
            y = y + A(k)*cos(2*pi*f(k)*t+ph(k));
        end
        
  
        %amp=[2.20E+02	1.40E+00	6.30E+00	1.00E+00	5.50E+00	6.00E-01	3.20E+00	2.00E-01	1.50E+00	1.00E-01	9.00E-01];
        % x=amp(1)*cos(x(12)*2*pi*t+x(1))+amp(2)*cos(x(13)*2*pi*t+x(2))+amp(3)*cos(x(14)*2*pi*t+x(3))+amp(4)*cos(x(15)*2*pi*t+x(4))+amp(5)*cos(x(16)*2*pi*t+x(5))+...
        %     +amp(6)*cos(x(17)*2*pi*t+x(6))+amp(7)*cos(x(18)*2*pi*t+x(7))+amp(8)*cos(x(19)*2*pi*t+x(8))+ amp(9)*cos(x(20)*2*pi*t+x(9))+amp(10)*cos(x(21)*2*pi*t+x(10))+...
        %     amp(11)*cos(x(22)*2*pi*t+x(11));
        x=220*cos(50.05*2*pi*t+x(1))+6.3*cos(50*3*2*pi*t+x(2))+5.5*cos(50*5*2*pi*t+x(3))+...
            +3.2*cos(50*7*2*pi*t+x(4))+1.5*cos(50*9*2*pi*t+x(5))+...
            0.9*cos(50*11*2*pi*t+x(6));
        % x= x(10)*cos(x(19)*2*pi*t+x(1))+x(11)*cos(x(20)*2*pi*t+x(2))+x(12)*cos(x(21)*2*pi*t+x(3))+x(13)*cos(x(22)*2*pi*t+x(4))+x(14)*cos(x(23)*2*pi*t+x(5))+...
        %     +x(15)*cos(x(24)*2*pi*t+x(6))+x(16)*cos(x(25)*2*pi*t+x(7))+x(17)*cos(x(26)*2*pi*t+x(8))+x(18)*cos(x(27)*2*pi*t+x(9)); 
        % Add Gaussian noise to signal x
        % noise = 0.01*randn(size(x));
        % x =x + noise;
        % x(16)*cos(x(17)*w0+pi*t+x(18))+x(19)*cos(x(20)*w0*pi*t+x(21));
        %y加噪声-x x为估计函数

        %%新的范例
        % x= 0.447*cos(19.75*2*pi*t+x(1))+x(2)*cos(2*19.75*2*pi*t+x(3))+x(4)*cos(3*19.75*2*pi*t+x(5))+x(6)*cos(4*19.75*2*pi*t+x(7));
        % y= 0.5*cos(19.75*2*pi*t+pi/9)+1*cos(2*19.75*2*pi*t+pi/6)+2*cos(3*19.75*2*pi*t+pi/3)+3*cos(4*19.75*2*pi*t+pi/2);

        %二范数
         f=norm(x-y,2);

end

    function R=GetInequality(equation)
        if equation<=0
            R=0;
        else
            R=1;
        end
    end
end
