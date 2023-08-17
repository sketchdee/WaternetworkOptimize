%% ObjectFunction.m 
%This file contain the objective funcion 
%      whitch will be used in optiomal algrithim. 
%
function [Z]=ObjectFunction(x)
%为了并行计算支持，将数据嵌入到脚本中，防止并发读取造成死锁（报错）
%----------------杂质负荷数据--------------------
MF=[
314.28	237.6	0
76750.32	20220	265165
11129.5	0	63771
15085.4	-6600	-169840
0	0	383040
-33737.5	0	-25250
2149.9	0	264
244	0	78336
285263.3	65877	-283077    
];
%---------------------end--------------------------

%------------------浓度数据------------------------
Cw=[18.9,0,133];
CinMax=[
100	0	300
300	10	2000
300	10	2000
300	500	2000
300	100	30
300	0	300
0.1	0	300
0.1	0	500
1000	1000	3000
];
CoutMax=[
300	30	1000
1600	100	3000
1600	100	3000
500	0	30
500	300	4500
0.1	0	500
20000	0	1000
10	0	3000
20000	20000	300
];
%---------------------end--------------------------

%---------------向量转矩阵形式------------------
M= Tool.V2M(x,12,12);
%最大允许误差，计算机浮点不存在完全的精确
tol=1/10000000;
%---------------------end--------------------------

%-------------计算单元进出口浓度----------------
A=ones(9,9);
%构造杂质守恒矩阵
for i=1:1:9
    A(i,:)=(M(1:9,i)*(-1))';
    A(i,i)=M(i,1:12)*ones(12,1);
end
%水网络合理性检查
if(rank(A)<9)
    disp('Invalid water network');
end
%计算1-9单元的进出口浓度，并对出现负值的调整杂质移除率。
%化水、综水、脱硫三个负杂质负荷修正
b1=MF(:,1)+M(12,1:9)'*Cw(1);
Ainv=inv(A);
C1_out=Ainv*b1;
%+----+----+----+----+
%|unit| Cl | SS | SO4|
%+----+----+----+----+
%| 4  |  × |  √ |  √ |
%+----+----+----+----+
%| 6  |  √ |  × |  √ |
%+----+----+----+----+
%| 9  |  × |  × |  √ |
%+----+----+----+----+
%+-----------+-------------------+-----------+
%|   Cl      |         SS        |    SO4    |
%+-----------+-------------------+-----------+
% 15085.4            -6600          -169840
%杂质1
%如果出现负浓度，则修正杂质负荷使得出口浓度为0
if(C1_out(6,1)<0)
    C16min=0;
    A_sub=A;
    A_sub(6,:)=[];
    A_sub(:,6)=[];
    b_sub=b1([1:5,7:9],1)-C16min*A(6,6);
    %计算其他浓度
    C1_out_s=inv(A_sub)*b_sub;
    C1_out(1:5,1)=C1_out_s(1:5,1);
    C1_out(6,1)=C16min;
    C1_out(7:9,1)=C1_out_s(6:8,1);
    b1(6,1)=A(6,:)*C1_out;
end
C1sum=(C1_out'*M(1:9,1:9)+M(12,1:9)*Cw(1))';
Qsum=(ones(1,12)*M(1:12,1:9))';
C1_in=C1sum./Qsum;

%杂质2
b2=MF(:,2)+M(12,1:9)'*Cw(2);
C2_out=inv(A)*b2;
%如果出现负浓度，则修正杂质负荷使得出口浓度为0
if(C2_out(4,1)<0)
    C24min=0;
    A_sub=A;
    A_sub(4,:)=[];
    A_sub(:,4)=[];
    b_sub=b2([1:3,5:9],1)-C24min*A(4,4);
    %计算其他浓度
    C2_out_s=inv(A_sub)*b_sub;
    C2_out =[C2_out_s(1:3,1);C24min;C2_out_s(4:8,1);];
    b2(4,1)=A(6,:)*C2_out;
end
C2sum=(C2_out'*M(1:9,1:9)+M(12,1:9)*Cw(2))';
C2_in=C2sum./Qsum;

%杂质3
b3=MF(:,3)+M(12,1:9)'*Cw(3);
C3_out=inv(A)*b3;
%如果出现负浓度，则修正杂质负荷使得出口浓度为0
if(C3_out(4,1)<0||C3_out(6,1)<0||C3_out(9,1)<0)
    C34min=0;
    C36min=0;
    C39min=0;
    A_sub=A;
    A_sub([4,6,9],:)=[];
    A_sub(:,[4,6,9])=[];
    b_sub=b3;
    b_sub([4,6,9],:)=[];
    rang=[1,2,3,5,7,8];
    b_sub=b_sub-A(rang,4)*C34min-A(rang,6)*C36min-A(rang,9)*C39min;
    C3_out_s=inv(A_sub)*b_sub;
    C3_out(rang)=C3_out_s;
    C3_out([4,6,9],:)=[0;0;0];
    b3([4,6,9],:)=A([4,6,9],:)*C3_out;

    %计算实际移除负荷
    %确定最大移质位置
    MF_3=b3-M(12,1:9)'*Cw(3);
    rang_b=[];
    if(MF_3(4,1)<MF(4,3))
        rang_b=[rang_b,4];
    end
    if(MF_3(6,1)<MF(6,3))
        rang_b=[rang_b,6];
    end
    if(MF_3(9,1)<MF(9,3))
       rang_b=[rang_b,9];
    end
    if(~isempty(rang_b))
        rang_c=setdiff([4,6,9],rang_b);
        A_sub=A;
        A_sub(rang_c,:)=[];
        A_sub(:,rang_c)=[];
        b_sub=MF(:,3)+M(12,1:9)'*Cw(3);
        b_sub(rang_c,:)=[];
        C3_out_s=inv(A_sub)*b_sub;
        rang_out_c=setdiff(1:9,rang_c);
        C3_out(rang_out_c)=C3_out_s;
        C3_out(rang_c,:)=ones(size(rang_c,2),1)*0;   %%TODO: the minimize c in each unit must be spetified.
        b3=A*C3_out;
    end
end
C3sum=(C3_out'*M(1:9,1:9)+M(12,1:9)*Cw(3))';
C3_in=C3sum./Qsum;
%---------------------end--------------------------

%---------------进出口浓度超出限制加罚---------------
%罚系数
a=1000;
Z=0;
%分别对三种杂质进行加罚
for i=1:1:9
    if(C1_in(i)<=CinMax(i,1)+tol)
        Z=Z+0;
    else
        Z=Z+a*abs(C1_in(i)-CinMax(i,1));
    end
    
    if(C2_in(i)<=CinMax(i,2)+tol)
        Z=Z+0;
    else
        Z=Z+a*abs(C2_in(i)-CinMax(i,2));
    end
     
     if(C3_in(i)<=CinMax(i,3)+tol)
        Z=Z+0;
    else
        Z=Z+a*abs(C3_in(i)-CinMax(i,3));
     end
    
     if(C1_out(i)<=CoutMax(i,1) + tol)
         Z=Z+0;
     else
         Z=Z+a*abs(C1_out(i)-CoutMax(i,1));
     end
     
     if(C2_out(i)<=CoutMax(i,2)+tol)
         Z=Z+0;
     else
         Z=Z+a*abs(C2_out(i)-CoutMax(i,2));
     end
     
     if(C3_out(i)<=CoutMax(i,3)+tol)
         Z=Z+0;
     else
         Z=Z+a*abs(C3_out(i)-CoutMax(i,3));
     end
end
%---------------------end--------------------------

%-------------其他用水系统进口约束--------------
%其他用水系统进水约束
C11_inMax=[300;1000;2000];
F11sum=ones(1,9)*M(1:9,11)+M(12,11);
C11_in1=(C1_out'*M(1:9,11)+M(12,11)*Cw(1))/F11sum;
C11_in2=(C2_out'*M(1:9,11)+M(12,11)*Cw(2))/F11sum;
C11_in3=(C3_out'*M(1:9,11)+M(12,11)*Cw(3))/F11sum;
if(C11_in1<=C11_inMax(1)+tol)
    Z=Z+0;
else
    Z=Z+a*abs(C11_inMax(1)-C11_in1);
end
if(C11_in2<=C11_inMax(2)+tol)
    Z=Z+0;
else
    Z=Z+a*abs(C11_inMax(2)-C11_in2);
end
if(C11_in3<=C11_inMax(3)+tol)
    Z=Z+0;
else
    Z=Z+a*abs(C11_inMax(3)-C11_in3);
end
%---------------------end--------------------------

%---------------目标函数排水量最小---------------
for i=1:1:11
    Z=Z+M(i,12);
end
%---------------------end--------------------------

%-----------------水流数目优化--------------------
Z=Z+sum(sum(round(M,2)~=0));
%---------------------end--------------------------
end