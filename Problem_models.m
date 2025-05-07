
function Fit=Problem_models(x,Number)

switch Number

    case 1 %压力容器设计(Pressure vessel design)
        lambda=10^20;
        g(1)=-x(1)+0.0193*x(3);
        g(2)=-x(2)+0.00954*x(3);
        g(3)= -pi*x(3)^2*x(4)-(4/3)*pi*x(3)^3+1296000;
        g(4)= x(4)-240;
        f=0.6224*x(1)*x(3)*x(4)+1.7781*x(2)*x(3)^2+3.1661*x(1)^2*x(4)+19.84*x(1)^2*x(3);
        Penalty=0;
        for i=1:length(g)
            Penalty= Penalty+ lambda*g(i)^2*GetInequality(g(i));
        end
        Fit=f+Penalty;

    case 2 %滚动轴承设计(Rolling element bearing design)
        x(3)=round(x(3));
        lambda=10^20;
        D=160;
        d=90;
        Bw=30;
        T=D-d-2*x(2);
        phio=2*pi-2*acos(((((D-d)/2)-3*(T/4))^2+(D/2-T/4-x(2))^2-(d/2+T/4)^2)...
            /(2*((D-d)/2-3*(T/4))*(D/2-T/4-x(2))));
        g(1)=-phio/(2*asin(x(2)/x(1)))+x(3)-1;
        g(2)=-2*x(2)+x(6)*(D-d);
        g(3)= -x(7)*(D-d)+2*x(2);
        g(4)=(0.5-x(9))*(D+d)-x(1);
        g(5)=-(0.5+x(9))*(D+d)+x(1);
        g(6)=-x(1)+0.5*(D+d);
        g(7)= -0.5*(D-x(1)-x(2))+x(8)*x(2);
        g(8)=x(10)*Bw-x(2);
        g(9)=0.515-x(4);
        g(10)=0.515-x(5);
        Penalty=0;
        for i=1:length(g)
            Penalty= Penalty+ lambda*g(i)^2*GetInequality(g(i));
        end
        gama=x(2)/x(1);
        fc=37.91*((1+(1.04*((1-gama/1+gama)^1.72)*((x(4)*(2*x(5)-1)/x(5)*...
            (2*x(4)-1))^0.41))^(10/3))^-0.3)*((gama^0.3*(1-gama)^1.39)/...
            (1+gama)^(1/3))*(2*x(4)/(2*x(4)-1))^0.41;
        if x(2)<=25.4
            f=-fc*x(3)^(2/3)*x(2)^1.8;
        else
            f=-3.647*fc*x(3)^(2/3)*x(2)^1.4;
        end
        Fit=f+Penalty;

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

        % %%Log-cosh函数
        % f1=0;
        % f2=0;
        % for i=1:length(t)
        %     if(t(i)<=1/12*0.02 || (t(i)>=5/12*0.02 && t(i)<=7/12*0.02) || (t(i)>=11/12*0.02 && t(i)<=0.02))
        %         f1=f1+(x(i)-y(i))^2/2;
        %     end
        %     if((t(i)>=1/12*0.02 && t(i)<=5/12*0.02) || (t(i)>=7/12*0.02 && t(i)<=11/12*0.02))
        %        f2=f2+abs(x(i)-y(i))-log10(2);
        %     end
        % end
        % f=f1+f2;

        % %% Quantile Loss 分位数损失
        % lambda=0.7;
        % N=200;
        % f1=0;
        % f2=0;
        % for i=1:length(t)
        %     if(abs(x(i)-y(i))>0)
        %         f1=1/N*(1-lambda)*abs(x(i)-y(i))+f1;
        %     end
        %     if(abs(x(i)-y(i))<=0)
        %         f2=1/N*(lambda)*abs(x(i)-y(i))+f2;
        %     end
        % end
        % f=f1+f2;

        % % % Huber-Loss函数
        % f1=0;
        % f2=0;
        % delta=0.05;
        % for i=1:length(t)
        %     if(abs(x(i)-y(i))>delta)
        %         f1=delta*abs(y(i)-x(i))-1/2*delta^2+f1;
        %     end
        %     if(abs(x(i)-y(i))<=delta)
        %         % f1=abs(x(i)-y(i))+f1;
        %         f2=1/2*(y(i)-x(i)).^2+f2;
        %     end
        % end
        % f=f1+f2;

        % for i=1:length(t)
        %     if(t(i)<=1/12*0.02 || (t(i)>=5/12*0.02 && t(i)<=7/12*0.02) || (t(i)>=11/12*0.02 && t(i)<=0.02))
        %         f2=norm((x(i)-y(i)),2)+f2;
        %     end
        %     if((t(i)>=1/12*0.02 && t(i)<=5/12*0.02) || (t(i)>=7/12*0.02 && t(i)<=11/12*0.02))
        %         % f1=abs(x(i)-y(i))+f1;
        %         f1=5*norm(x(i)-y(i),2)+f1;
        %     end
        % end
        % f=f1+f2;
        % g(1)=w1-5;
        % g(2)=w2-5;
        % g(3)=w3-5;
        % g(4)=w4-5;
        % g(5)=w5-5;
        % g(6)=-f1+70;
        % g(7)=f1-75;
        Fit=f;
        % t=0:0.2/1023:0.2;
        % x=100*cos(50*2*pi*t)+30*cos(100*2*pi*t)+80*cos(150*2*pi*t)+15*cos(200*2*pi*t)+20*cos(250*2*pi*t)...
        %     +8*cos(300*2*pi*t)+10*cos(350*2*pi*t)+2*cos(400*2*pi*t)+5*cos(450*2*pi*t)+cos(500*2*pi*t); %参考项,10次谐波
        % y=100*cos(50*1.2*2*pi*t)+30*cos(100*0.9*2*pi*t)+80*cos(150*1.15*2*pi*t)+15*cos(200*0.98*2*pi*t)+20*cos(1.02*250*2*pi*t)...
        %     +8*cos(300*1.04*2*pi*t)+10*cos(350*2*pi*t)+2*cos(400*1.1*2*pi*t)+5*cos(450*2*pi*t)+cos(0.99*500*2*pi*t); %参考项,10次谐波
% plot(t,x,'b');
% grid on
% hold on
% plot(t,y,'ro-');
        

    case 4 %悬臂梁设计(Cantilever beam design)
        lambda=10^20;
        g(1)=61/x(1)^3+37/x(2)^3+19/x(3)^3+7/x(4)^3+1/x(5)^3-1;
        f=0.0624*(x(1)+x(2)+x(3)+x(4)+x(5));
        Penalty=0;
        for i=1:length(g)
            Penalty= Penalty+ lambda*g(i)^2*GetInequality(g(i));
        end
        Fit=f+Penalty;

    case 5 %轮系设计(Gear train design)
        x=round(x);
        Fit=(1/6.931-(x(3)*x(2)/(x(1)*x(4))))^2;

    case 6 %三杆桁架设计(Three bar truss design) 
        lambda=10^20;
        l=100;
        P=2;
        o=2;
        g(1)=((sqrt(2)*x(1)+x(2))/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*P-o;
        g(2)=(x(2)/(sqrt(2)*x(1)^2+2*x(1)*x(2)))*P-o;
        g(3)=(1/(sqrt(2)*x(2)+x(1)))*P-o;
        f=(2*sqrt(2)*x(1)+x(2))*l;
        Penalty=0;
        for i=1:length(g)
            Penalty= Penalty+ lambda*g(i)^2*GetInequality(g(i));
        end
        Fit=f+Penalty;

end

    function R=GetInequality(equation)
        if equation<=0
            R=0;
        else
            R=1;
        end
    end
end