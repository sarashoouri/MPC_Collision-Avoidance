function[feas, xIter,uIter, JIter] = MPCS(u2d,u1d,qx,rou1,rou2,Tx,M,Q, R, N, x0, xL, xU, uL, uU,obstacle,safetyR)
x = sdpvar(4,N+1); %state variable
u = sdpvar(1,N); 
objective = x(:,N+1)'*Q*x(:,N+1)+qx*(x(1,N)+x(3,N)*cos(x(4,N))*0.2-x(1,N+1))^2;%check the report
objective = objective+x(:,1)'*Q*x(:,1)+u(1)'*R*u(1);
for i=2:N
    objective=objective+x(:,i)'*Q*x(:,i)+u(i)'*R*u(i)+qx*(x(1,i-1)+x(3,i-1)*cos(x(4,i-1))*0.2-x(1,i))^2;
end
constraints = [];
constraints = [constraints x(:,1)==x0];
constraints = [constraints xL<=x(:,N+1)<=xU];
for i = 1:N
    constraints = [constraints  xL<=x(:,i)<=xU -0.6<=u(i)<=0.6 ];
    constraints = [constraints  x(:,i+1) == modelS(x(:,i),u(i))];
    if i <= N-1       
       constraints = [constraints -0.05<=(u(i+1)-u(i))<=0.05];
    end
end
options = sdpsettings('verbose',0,'solver','IPOPT');
optimize(constraints, objective, options);
diagnostics=optimize(constraints, objective, options);
if diagnostics.problem == 0 
    xIter{1}=value(x);
    uIter{1}=value(u); 
    JIter{1}=value(objective);
    feas=true;
   
else
    feas=false;
    xIter{1}=[];
    uIter{1}=[];
    JIter{1}=[];
    return
end
% applying SCP algorithm
for xj=xIter{1}(1,:)
    if ((xj<=obstacle(1)+obstacle(3)) && (xj>=obstacle(1)-obstacle(3)))||...
            ((xj<=obstacle(4)+obstacle(6)) && (xj>=obstacle(4)-obstacle(6)))||...
            ((xj<=obstacle(7)+obstacle(9)) && (xj>=obstacle(7)-obstacle(9)))||...
            ((xj<=obstacle(10)+obstacle(12)) && (xj>=obstacle(10)-obstacle(12)))
        k=1;
        tic
        while 1
            x = sdpvar(4,N+1);
            slack =sdpvar(1,N+1);
            u = sdpvar(1,N);
            objective = x(:,N+1)'*Q*x(:,N+1)+qx*(x(1,N)+x(3,N)*cos(x(4,N))*0.2-x(1,N+1))^2+slack(N+1)*rou2;
            objective = objective+x(:,1)'*Q*x(:,1)+u(1)'*R*u(1)+slack(1)*rou2;
            for i=2:N
                objective=objective+x(:,i)'*Q*x(:,i)+u(i)'*R*u(i)+qx*(x(1,i-1)+x(3,i-1)*cos(x(4,i))*0.2-x(1,i))^2+slack(i)*rou2;
            end
            constraints = [];
            constraints = [constraints x(:,1)==x0 slack(1)>=0];
            constraints = [constraints xL<=x(:,N+1)<=xU slack(N+1)>=0];          
            for i = 1:N
                constraints = [constraints  xL<=x(:,i)<=xU -0.6<=u(i)<=0.6 slack(i)>=0];
                constraints = [constraints  x(:,i+1) == bikeFE_Noa(x(:,i),u(i))];
                constraints = [constraints  (xIter{k}(1,i)-obstacle(1))^2+(xIter{k}(2,i)-obstacle(2))^2+(2*xIter{k}(1,i)-2*obstacle(1))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(2))*(x(2,i)-xIter{k}(2,i))>=(obstacle(3)+safetyR)^2-slack(i)];
                constraints = [constraints  (xIter{k}(1,i)-obstacle(4))^2+(xIter{k}(2,i)-obstacle(5))^2+(2*xIter{k}(1,i)-2*obstacle(4))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(5))*(x(2,i)-xIter{k}(2,i))>=(obstacle(6)+safetyR)^2-slack(i)];
                constraints = [constraints  (xIter{k}(1,i)-obstacle(7))^2+(xIter{k}(2,i)-obstacle(8))^2+(2*xIter{k}(1,i)-2*obstacle(7))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(8))*(x(2,i)-xIter{k}(2,i))>=(obstacle(9)+safetyR)^2-slack(i)];
                constraints = [constraints  (xIter{k}(1,i)-obstacle(10))^2+(xIter{k}(2,i)-obstacle(11))^2+(2*xIter{k}(1,i)-2*obstacle(10))*(x(1,i)-xIter{k}(1,i))+(2*xIter{k}(2,i)-2*obstacle(11))*(x(2,i)-xIter{k}(2,i))>=(obstacle(12)+safetyR)^2-slack(i)];                
                if i <= N-1       
                   constraints = [constraints -0.05<=(u(i+1)-u(i))<=0.05];
                end
            end
            options = sdpsettings('verbose',0,'solver','IPOPT');          
            optimize(constraints, objective, options);
            diagnostics=optimize(constraints, objective, options);
            if diagnostics.problem == 0 
                xIter{k+1}=value(x);
                uIter{k+1}=value(u); 
                JIter{k+1}=value(objective);
                feas=true;
            else
                feas=false;
                return
            end
            if abs(JIter{k+1}-JIter{k})<1
                
                break
            end
            if k>25
                break
            end
            k=k+1;
        end
        toc
        break
    else
        continue
    end
end
end