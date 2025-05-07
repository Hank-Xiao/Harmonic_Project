function [Lu,Up,Dim]=Problem_range(Number)


switch Number

    case 1 %压力容器设计(Pressure vessel design)

        Lu=[0 0 10 10];
        Up=[99 99 200 200];
        Dim=4;

    case 2 %滚动轴承设计(Rolling element bearing design)
        D=160;
        d=90;
        Lu=[0.5*(D+d) 0.15*(D-d) 4 0.515 0.515 0.4 0.6 0.3 0.02 0.6];
        Up=[0.6*(D+d) 0.45*(D-d) 50 0.6 0.6 0.5 0.7 0.4 0.1 0.85];
        Dim=10;

    case 3 %拉伸/压缩弹簧设计(Tension/compression spring design)
        % Lu=[0 0 0 0 0 0 0 0 0 0 0 0 0 130 0 0 270 0 0 330 0];
        % Up=[2*pi 25 2*pi 25 2*pi 25 2*pi 25 2*pi 25 50 2*pi 25 230 2*pi 25 330 2*pi 25 430 2*pi];
        % Dim=21;
        %3.2 380 20 2.8 16.4 1.8 12 1.4 10
        %25 50 150 180 250 320 350 380 450
        Lu=[0 0 0 0 0 0];
        Up=[2*pi 2*pi 2*pi 2*pi 2*pi 2*pi];
        Dim=6; 
        % Up=[pi 5 pi 5 pi 5 pi];
        % Dim=7;
        % Lu=[0 0 0 0 0 0 0 0 0];
        % Up=[pi 25 pi 25 pi 25 pi 25 pi];
        % Dim=9;

    case 4 %悬臂梁设计(Cantilever beam design)
        Lu=0.01;
        Up=100;
        Dim=5;


    case 5 %轮系设计(Gear train design)
        Lu=12;
        Up=60;
        Dim=4;

    case 6 %三杆桁架设计(Three bar truss design) 
        Lu=0;
        Up=1;
        Dim=2;

end