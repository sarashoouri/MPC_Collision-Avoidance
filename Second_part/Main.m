clear all;
%%
load('MPC.mat'); % from first part
xled = xOpt(1,:);
yled = xOpt(2,:);
vled = xOpt(3,:);
thetaled = xOpt(4,:);
TS = 0.2;
lr = 3;
lf = 3;
nx = 101;
timestep = 100; 
zref = [xled; yled; vled; thetaled];
N = size(xled,2)-1;
P = 3;
x0 = 0; 
y0 = 0; 
v0 = 0; 
theta0 = 0;
t0 = 0; 
z0B = [x0; y0; v0; theta0]; 
nzB = size(z0B,1); 
a0 = 0; 
deltaf0 = 0; 
u0B = [a0; deltaf0];
nuB = size(u0B,1);
zB = sdpvar(nzB,N+1);
assign(zB(:,1), z0B);
uB = sdpvar(nuB, N);
assign(uB(:,1), u0B);
QB = 12*eye(nzB);
QB(3,3) = 3; 
RB = [0 0; 0 0];
usaveB = zeros(nuB,N);
zsaveB = zeros(nzB,N);
zsaveB(:,1) = z0B;
z0C = [x0; y0; v0; theta0]; 
nzC = size(z0C,1); 
u0C = [a0; deltaf0];
nuC = size(u0C,1);
nzC = size(z0C,1); 

zC = sdpvar(nzC,N+1);
assign(zC(:,1), z0C);
uC = sdpvar(nuC, N);
assign(uC(:,1), u0C);
QC = 12*eye(nzC);
QC(3,3) = 4; 
RC = [0.01 0; 0 0];
usaveC = zeros(nuC,N);
zsaveC = zeros(nzC,N);
zsaveC(:,1) = z0C;
zmin = [-10;-100;0;-pi];
zmax = [1000;100;30;pi];
umin = [-1.5;-60*pi/180];
umax = [4;60*pi/180];
%%
for i=1:N
    if i<=size(yled)-10
        sumy = sum(yled(i:i+10));
    else 
        sumy = sum(yled(i+5:N));
    end
    
    threshold = 3;
    if sumy>=threshold
        safe = [-0.03*(lr+lf);0*(lr+lf);0;0];
    elseif sumy<=-threshold
        safe = [-0.03;0*(lr+lf);0;0];
    else 
        safe= [0;0;0;0];
    end

    bar_zrefB = zref(:,i);
    objB = 0;
    if i<=10
        usaveB(:,i) = zeros(nuB,1);
        zsaveB(:,i+1) = zeros(nzB,1);
    else
        bar_zrefB = zref(:,i-10);
        for j =1:P
            if abs(zsaveB(1,i)-zref(1,i-10))<=3*(lr+lf) || abs(zsaveB(2,i)-zref(2,i-10)<=3*(lr+lf))
                objB = objB + (zB(:,j)-(bar_zrefB+safe))'*QB*(zB(:,j)-(bar_zrefB+safe)) + uB(:,j)'*RB*uB(:,j);
            else
                objB = objB + (zB(:,j)-(bar_zrefB+safe))'*QB*(zB(:,j)-(bar_zrefB+safe)) + uB(:,j)'*RB*uB(:,j);
            end

            bar_zrefB(1) = bar_zrefB(1)+TS*bar_zrefB(3)*cos(bar_zrefB(4));
            bar_zrefB(2) = bar_zrefB(2)+TS*bar_zrefB(3)*sin(bar_zrefB(4));
        end
        constraintB = zB(:,1) == zsaveB(:,i);
        for j = 1:P
            betaB(j) = atan((lr*tan(uB(2,j)))/(lf+lr));
            constraintB = [constraintB, zmin <= zB(:,j) <= zmax,...
                umin <= uB(:,j) <= umax,...
                zB(1,j+1) == zB(1,j)+TS*zB(3,j)*cos(zB(4,j)+betaB(j)),...
                zB(2,j+1) == zB(2,j)+TS*zB(3,j)*sin(zB(4,j)+betaB(j)),...
                zB(3,j+1) == zB(3,j)+TS*uB(1,j),...
                zB(4,j+1) == zB(4,j)+TS*zB(3,j)*sin(betaB(j))/lr];
        end
        options = sdpsettings('verbose',0);
        sol = optimize(constraintB, objB, options); 
        usaveB(:,i) = uB(:,1);
        zsaveB(:,i+1) = zB(:,2);
    end
    if i<=size(yled)-10
        sumy = sum(yled(i:i+10));
    else 
        sumy = sum(yled(i+5:N));
    end
    
    threshold = 3;
    if sumy>=threshold
        safe = [-0.03*(lr+lf);0*(lr+lf);0;0];
    elseif sumy<=-threshold
        safe = [-0.03;0*(lr+lf);0;0];
    else 
        safe= [0;0;0;0];
    end
    bar_zrefC = zsaveB(:,i);
    objC = 0;
    if i<=20
        usaveC(:,i) = zeros(nuC,1);
        zsaveC(:,i+1) = zeros(nzC,1);
    else
        bar_zrefC = zsaveB(:,i-10);
        for j =1:P
            if abs(zsaveC(1,i)-zsaveB(1,i-10))<=3*(lr+lf) || abs(zsaveC(2,i)-zsaveB(2,i-10)<=3*(lr+lf))
                objC = objC + (zC(:,j)-(bar_zrefC+safe))'*QC*(zC(:,j)-(bar_zrefC+safe)) + uC(:,j)'*RC*uC(:,j);
            else
                objC = objC + (zC(:,j)-(bar_zrefC))'*QC*(zC(:,j)-(bar_zrefC)) + uC(:,j)'*RB*uC(:,j);
            end

            bar_zrefC(1) = bar_zrefC(1)+TS*bar_zrefC(3)*cos(bar_zrefC(4));
            bar_zrefC(2) = bar_zrefC(2)+TS*bar_zrefC(3)*sin(bar_zrefC(4));
        end
        constraintC = zC(:,1) == zsaveC(:,i);
        for j = 1:P
            betaC(j) = atan((lr*tan(uC(2,j)))/(lf+lr));
            constraintC = [constraintC, zmin <= zC(:,j) <= zmax, umin <= uB(:,j) <= umax,zC(1,j+1) == zC(1,j)+TS*zC(3,j)*cos(zC(4,j)+betaC(j)),...
                zC(1,j+1) == zC(1,j)+TS*zC(3,j)*cos(zC(4,j)+betaC(j)),...
                zC(2,j+1) == zC(2,j)+TS*zC(3,j)*sin(zC(4,j)+betaC(j)),...
                zC(3,j+1) == zC(3,j)+TS*uC(1,j),...
                zC(4,j+1) == zC(4,j)+TS*zC(3,j)*sin(betaC(j))/lr];
        end
        options = sdpsettings('verbose',0);
        solC = optimize(constraintC, objC, options); 
        usaveC(:,i) = uC(:,1);
        zsaveC(:,i+1) = zC(:,2);
    end
end
%%
time = linspace(0,timestep,nx-1);
h=figure;
plot(time,usaveB(1,:))
title("acceleration of carB")
ylabel('acceleration ');
xlabel('Time [s]')
pubgraph(h,20,3,'w')

%%
h=figure;
plot(time,usaveC(1,:))
title("acceleration input of carC")
ylabel('acceleration ');
xlabel('Time [s]')
pubgraph(h,20,3,'w')
%%
h=figure;
plot(xled(1,:), yled(1,:),'-o')
hold on;
plot(zsaveB(1,:),zsaveB(2,:), '-x')
hold on;
plot(zsaveC(1,:),zsaveC(2,:), '-x')
xlabel('Time [s]')
legend("carA","carB","carC");
title("carA's vs carB's carC's positions")
pubgraph(h,20,3,'w')
hold off;
%%
t = linspace(1,100,100);

h=figure;
plot(t,zsaveB(1,1:100))
hold on
plot(t,zsaveC(1,1:100))
legend("car B","car C")
xlabel('Time [s]')
title("X Coordinates for Car B and Car C")
pubgraph(h,20,3,'w')
%%
h=figure;
plot(t,zsaveB(2,1:100))
hold on
plot(t,zsaveC(2,1:100))
xlabel('Time [s]')
title("Y Coordinates for Car B and Car C")
legend("car B","car C")
pubgraph(h,20,3,'w')
 %%
h=figure;
plot(t,zsaveB(3,1:100))
hold on
plot(t,zsaveC(3,1:100))
xlabel('Time [s]')
title("Speed for Car B and Car C")
legend("car B","car C")
pubgraph(h,20,3,'w')
%%
h=figure;
plot(t,zsaveB(4,1:100))
hold on
plot(t,zsaveC(4,1:100))
xlabel('Time [s]')
title("Heading angle for Car B and Car C")
legend("car B","car C")
pubgraph(h,20,3,'w')    
