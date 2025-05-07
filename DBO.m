function [Best_pos,Best_score,curve] = DBO(Number,T,N)
%% Problem Definition

[lb,ub,D]=Problem_range(Number);

VarSize=[1 D];   % ���߱��������С

k=0.1;
b=0.3;
S=0.5;

%% ��ʼ����Ⱥ

habitat.Position=[];
habitat.Cost=[];

pop=repmat(habitat,N,1);

for i=1:N
    pop(i).Position=rand(1,D).*(ub-lb)+lb;
    pop(i).Cost=Problem_models(pop(i).Position,Number);
end

% ������Ӧ��ֵ����
[~, SortOrder]=sort([pop.Cost]); % sort(A) �����򣨴�С���󣩶� A ��Ԫ�ؽ�������
pop=pop(SortOrder);

% ��������ÿ��������ֵ
BestCost=zeros(T,1);
BestCost(1)=pop(1).Cost;


%% SO Main Loop

for it=2:T

    newpop=pop;  %������ʱ��Ⱥ

    if it==2
        oldpop=pop;
    end


    
    for i=1:N

        if i<=(6/30)*N  %ball-rolling dung beetle
            if rand<0.9
                % Algorithm1
                if rand()>0.5
                    alpha=1;
                else
                    alpha=-1;
                end
                % (1)ʽ
                newpop(i).Position=pop(i).Position+alpha*k*oldpop(i).Position+b*(pop(i).Position-pop(N).Position);
            else
                thea=unifrnd(0,pi);
                if thea~=0&&thea~=pi/2&&thea~=pi % Algorithm2
                    % (2)ʽ��
                    newpop(i).Position=pop(i).Position+tan(thea)*(pop(i).Position-oldpop(i).Position);
                end
            end
        end
        if i>(6/30)*N&&i<=(12/30)*N % brood ball
            % Algorithm3
            R=1-it/T;
            Lb=max(pop(1).Position*(1-R),lb); % (3)ʽ
            Ub=min(pop(1).Position*(1+R),ub); % (3)ʽ
            b1=rand(1,D);
            b2=rand(1,D);
            newpop(i).Position=pop((6/30)*N+1).Position+b1.*(pop(i).Position-Lb)+b2.*(pop(i).Position-Ub); % (4)ʽ

            for j=1:D
                if newpop(i).Position(j)>Ub(j)
                    newpop(i).Position(j)=Ub(j);
                end
                if newpop(i).Position(j)<Lb(j)
                    newpop(i).Position(j)=Lb(j);
                end
            end
        end
        if i>(12/30)*N&&i<=(19/30)*N % small dung beetle
            R=1-it/T;
            Lb=max(pop(1).Position*(1-R),lb); % (5)ʽ
            Ub=min(pop(1).Position*(1+R),ub); % (5)ʽ
            C1=normrnd(0,1);
            C2=rand(1,D);
            newpop(i).Position=pop(i).Position+C1*(pop(i).Position-Lb)+C2.*(pop(i).Position-Ub); % (6)ʽ
            % l=2*rand(1,D)-1;
            % r=1;
            % newpop(i).Position=pop(i).Position+exp(r*l).*(cos(2*pi*l)).*C1.*(pop(i).Position-Lb)+exp(r*l).*(cos(2*pi*l)).*C2.*(pop(i).Position-Ub); % (6)ʽ
        end
        if i>(19/30)*N % some thieves
            g=normrnd(0,1,VarSize);
            newpop(i).Position=pop(1).Position+S*g.*((pop((19/30)*N+1).Position-pop(i).Position)+(pop(1).Position)-pop(i).Position);
        end

         if it <= T/4
               for j=1:D
                    % Horizontal crossover search (HCS)
                    e1=rand(); e2=rand(); c1=rand(); c2=rand();
                    k=round(unifrnd(1,N)); % ���ѡ��һ������
                    newpop(i).Position(j)=e1*newpop(i).Position(j)+(1-e1)*newpop(k).Position(j)+c1*(newpop(i).Position(j)-newpop(k).Position(j)); % �ο������е�Eq. (17)
                    newpop(k).Position(j)=e2*newpop(k).Position(j)+(1-e2)*newpop(i).Position(j)+c2*(newpop(k).Position(j)-newpop(i).Position(j)); % �ο������е�Eq. (18)
                    if j==1
                        % Vertical crossover search (VCS)
                        e=rand();
                        newpop(i).Position(1)=e*newpop(i).Position(1)+(1-e)*newpop(i).Position(2); % �ο������е�Eq. (19)
                    end
                end
            end

        %��֤ÿ��������ȡֵ��Χ����

        newpop(i).Position=Bound_limit(newpop(i).Position,ub,lb);

        for j=1:D
            % �ο������е�Eq. (16)
            if newpop(i).Position(j)<lb(j)
                newpop(i).Position(j) = max(lb(j)*Levy(1), lb(j));
            end
            if newpop(i).Position(j)>ub(j)
                newpop(i).Position(j) = min(ub(j)*Levy(1), ub(j));
            end
        end

        %��������
        newpop(i).Cost=Problem_models(newpop(i).Position,Number);

        %̰��ѡ��
        if newpop(i).Cost<pop(i).Cost
            pop(i).Position=newpop(i).Position;
            pop(i).Cost=newpop(i).Cost;
        end


    end

    oldpop=newpop;

    [~, SortOrder]=sort([pop.Cost]);
    pop=pop(SortOrder);

    BestCost(it)=pop(1).Cost;
    BestSol=pop(1).Position;

    %��ʾÿ���ҵ�������ֵ
    format long e;
    disp(['Iteration' num2str(it)]);
    fprintf('Best Cost is��%40.30f\n',BestCost(it));

end
%% Results ��������ս����
curve = BestCost;
format long e;
Best_score =min(BestCost);
Best_pos = BestSol;
end