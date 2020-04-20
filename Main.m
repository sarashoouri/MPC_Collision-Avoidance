clear all 
clear('yalmip')
%% parameters values
Q=0.1*[0 0 0 0;
    0 1 0 0;
    0 0 0 0;
    0 0 0 0];
rou1=100;
rou2=1000;
qx=0.1;
x0=[0;
    0;
    5;
    0];
xL=[0;
    -4.5;
    0;
    -90*pi/180];
xU=[1000;
    4.5;
    30;
    90*pi/180];
uL=[-5.4;
    -0.6];
uU=[2.7;
    0.6];
u2d=0.05;
u1d=0.5;
% set the obstacles (considering that obstacles are circular).
obstacle=[10;0;1;25;1;1.5;40;0;1;55;-1;1.5];
safetyR=2; % Safety region
T1=linspace(0,22,23);
T2=linspace(23,27,10);
T3=linspace(28,72,45);
T4=linspace(73,77,10);
T5=linspace(78,1000,923);
Tx=cat(2,T1,T2,T3,T4,T5);
tic
N=5;
M=100;
R=1;
%carry out MPC on a simulation horizon of M
feas=zeros(1,M);
xOpt=zeros(size(Q,2),M+1);xOpt(:,1)=x0;
uOpt=zeros(size(R,2),M);
JOpt=zeros(1,M);
pred=zeros(size(Q,2),N+1,M+1);
%%
K=1;
for i = 1:M
[feas(i), xIter,uIter, JIter] =MPCS(u2d,u1d,qx,rou1,rou2,Tx,i,Q, R, N, x0, xL, xU, uL, uU,obstacle,safetyR);
if feas(i)==false
    for t= i+1:M
        feas(t)=false;
    end
    return
end
x0=xIter{end}(:,2);
xOpt(:,i+1)=x0;
JOpt(i)=JIter{end};
uOpt(i)=uIter{end}(1);
for t=1:N+1
    pred(:,t,i)=xIter{end}(:,t);
end
feas = logical(feas);
K=K+1;
end
predErr=zeros(2,M-N+1);
DifferenceX1=zeros(1,N-1);
DifferenceX2=zeros(1,N-1);
for i=1:M-N+1
    for t=1:N-1
        DifferenceX1(t)= pred(1,t+2,i)-xOpt(1,i+t+1);
        DifferenceX2(t)= pred(2,t+2,i)-xOpt(2,i+t+1);
    end
predErr(1,i)=norm(DifferenceX1);
predErr(2,i)=norm(DifferenceX2);
end
%%
h=figure
plot(xOpt(1,:),xOpt(2,:),'-o')
xlabel('X coordinates ');
ylabel('Y coordinates ');
title('Obstacle Collision Avoidance for First Car')
hold on
for i=1:M+1
    plot(pred(1,:,i),pred(2,:,i),'--')
end
% title(['N=' num2str(N)])
legend('Closed loop trajectory','N step prediction at each time')
viscircles([obstacle(1) obstacle(2)],safetyR+obstacle(3),'LineStyle','--','color','b')
viscircles([obstacle(1) obstacle(2)],obstacle(3))
viscircles([obstacle(4) obstacle(5)],safetyR+obstacle(6),'LineStyle','--','color','b')
viscircles([obstacle(4) obstacle(5)],obstacle(6))
viscircles([obstacle(7) obstacle(8)],safetyR+obstacle(9),'LineStyle','--','color','b')
viscircles([obstacle(7) obstacle(8)],obstacle(9))
viscircles([obstacle(10) obstacle(11)],safetyR+obstacle(12),'LineStyle','--','color','b')
viscircles([obstacle(10) obstacle(11)],obstacle(12))
axis equal
hold off
pubgraph(h,20,3,'w')