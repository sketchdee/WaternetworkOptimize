%处理约束%
m=12;
n=12;

FL=[53.3;2176;26;0;0;0;0;59;385;45;102;2846.30];
%简化模型

%不自身回用
Aeq1=[];
beq1=[];
for i=1:1:m
    M1=zeros(m,n);
    M1(i,i)=1;
    Aeq1=[Aeq1;Tool.M2V(M1)];
    beq1=[beq1;0];
end

%2个终端单元
for i=1:1:n
    M1=zeros(m,n);
    M1(10,i)=1;
    Aeq1=[Aeq1;Tool.M2V(M1)];
    beq1=[beq1;0];
end
for i=1:1:n
    M1=zeros(m,n);
    M1(11,i)=1;
    Aeq1=[Aeq1;Tool.M2V(M1)];
    beq1=[beq1;0];
end

%两个处理单元
for i=1:1:m
    if(i~=4)
        M1=zeros(m,n);
        M1(i,5)=1;
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
    if(i==4)
        M1=zeros(m,n);
        M1(4,5)=-1/0.431818;
        M1(1:12,4)=ones(12,1);
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
end
for i=1:1:m
    if(i~=6)
        M1=zeros(m,n);
        M1(i,7)=1;
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
    if(i==6)
        M1=zeros(m,n);
        M1(6,7)=-1/0.008;
        M1(1:12,6)=ones(12,1);
        Aeq1=[Aeq1;Tool.M2V(M1)];
        beq1=[beq1;0];
    end
end

%流量守恒约束
Aeq2=[];
beq2=[];
for i=1:1:m-1
    M1=zeros(m,n);
    M1(i,:)=ones(1,n)*(-1);
    M1(:,i)=ones(m,1);
    M1(i,i)=0;
    Aeq2=[Aeq2;Tool.M2V(M1)];
    beq2=[beq2;FL(i,1)];
end

%锅炉用水约束
Q_boiler=110;
for i=1:1:m
    if(i~=6)
        M1=zeros(m,n);
        M1(i,8)=1;
        Aeq2=[Aeq2;Tool.M2V(M1)];
        beq2=[beq2;0];
    end
    if(i==6)
        M1=zeros(m,n);
        M1(i,8)=1;
        Aeq2=[Aeq2;Tool.M2V(M1)];
        beq2=[beq2;110];
    end
end

%生活用水只使用新鲜水和除盐水
for i=1:1:m
    if(i==6)
        continue;
    end
    if(i==12)
        continue;
    end
    M1=zeros(m,n);
    M1(i,1)=1;
    Aeq2=[Aeq2;Tool.M2V(M1)];
    beq2=[beq2;0];
end

Aeq=[Aeq1;Aeq2];
beq=[beq1;beq2];
ub_M=ones(m,n)*500;
ub_M(:,2)=ones(m,1)*(2500);
ub=Tool.M2V(ub_M)';
lb=Tool.M2V(ones(m,n)*-1e-10)';

Aineq=[];
bineq=[];
%最大入口流量约束
%生活用水小于100
Fmax=[100;3500; 300;300;300;1000;1000;1000;1000];
for i=1:1:9
    M1=zeros(m,n);
    M1(i,:)=ones(n,1);
    Aineq=[Aineq;Tool.M2V(M1)];
    bineq=[bineq;Fmax(i,1)];
end


%去除多余的等式约束
Aeq_beq=[Aeq,beq];
if(rank(Aeq)~=rank(Aeq_beq))
    disp('存在冲突的等式约束，无可行解，请检查等式约束！');
end
T= rref([Aeq_beq]);
r=rank(T);
Aeq=T(1:r,1:144);
beq=T(1:r,145);

options = optimoptions('fmincon');
options = optimoptions(options,...
    'ConstraintTolerance',1e-10,...
    'Display', 'iter',...
    'FiniteDifferenceType','central',...
    'FunctionTolerance',1e-10,...
    'MaxFunctionEvaluations',2e5,...
    'MaxIterations', 20000,...
    'ObjectiveLimit',-1e100,...
    'OptimalityTolerance',1e-10,...
    'StepTolerance',1e-500,...
    'UseParallel',false,...
    'Algorithm', 'interior-point',...
    'TolConSQP',1e-100);

x0=ones(1,144)*6.37;
%[x,fval,exitflag,output] = fmincon(@ObjectFunction,x0,Aineq,bineq,Aeq,beq,lb,ub,[],options);
%result=Tool.V2M(x,m,n);

%%流量守恒测试函数
function [Z]=ObjectFunction1(x)
M= Tool.V2M(x,12,12);
Z=0;
for i=1:1:11
    Z=Z+M(12,i);
end
end