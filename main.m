clc;
clear;
digits(45);

Problem_name=[{'Pressure vessel design'},{'Rolling element bearing design'},...
    {'Tension/compression spring design'},{'Cantilever beam design'},...
    {'Gear train design'},{'Three bar truss design'}];

MaxIt=1000; %最大迭代次数
PopSize=60; %种群规模，注意：DBO的种群只能取30的倍数

Number=3; %问题编号，自行改动：1 2 3 4 5 6

[BestX,BestF,Curve]=DBO(Number,MaxIt,PopSize);

% if Number==2 || Number==3
%     BestX(2)=round(BestX(2));
% elseif  Number==5
%     BestX=round(BestX);
% end


format long e;
disp(['The best fitness of ',Problem_name{Number},' is: ']);
fprintf('%40.30f\n', BestF);
display(['The best solution is: ', regexprep(num2str(BestX),'\s*',';')]);


% 画图(收敛曲线)

if Number==2
    Curve=-Curve;
    %滚动轴承设计(Rolling element bearing design)是最大化问题，处理成了最小化问题，因此收敛曲线需要再加一个负号

    plot(Curve,'R','LineWidth',2);
else
    semilogy(Curve,'R','LineWidth',2);
end

hold on
xlabel('Iterations');
ylabel('Fitness');
title('损失函数');
grid on;
legend('二范数');