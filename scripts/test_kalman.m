% System matrices 
A=[-1 2;-3 0];
B=[1;1];
C=[1 2];
if rank(obsv(A,C))==size(A,1)
    disp('The system is fully observable')
else
    disp('The system is not fully observable')
end

% Covariance matrices
Qtilde=[0.05 0;0 0.05];
Qtilde_tuning=Qtilde;
Rtilde=10;
P0=10*eye(2);

% Initial conditions
x0=[2;2];
%x0hat=x0+(randn(1,2)*chol(P0))';
x0hat=[5.7505;6.8004];
T=20;

tspan=[0 T];
[t1,x_tot] = ode45(@(t,x)LTIdynamics_with_observer_finite(t,x,A,B,C,sin(t),P0,Qtilde,Rtilde,T), tspan, [x0;x0hat]);




figure(2)
subplot(2,1,1)
plot(t1,x_tot(:,1),'k')
hold on
plot(t1,x_tot(:,3),'b')
hold on
grid on 
box on
xlabel('t')
ylabel('x_1')
legend('Real state','Kalman estimate')
axes('position',[0.715,0.595,0.185,0.1])
box on % put box around new pair of axes
indexOfInterest = (t1>=15) & (t1 <= 20);
plot(t1(indexOfInterest),x_tot(indexOfInterest,1),'k')
hold on
plot(t1(indexOfInterest),x_tot(indexOfInterest,3),'b')
grid on
box on
zoom('on')
subplot(2,1,2)
plot(t1,x_tot(:,2),'k')
hold on
plot(t1,x_tot(:,4),'b')
grid on 
box on
xlabel('t')
ylabel('x_2')
legend('Real state','Kalman estimate')
axes('position',[0.715,0.12,0.185,0.12])
box on % put box around new pair of axes
indexOfInterest = (t1>=15) & (t1 <= 20);
plot(t1(indexOfInterest),x_tot(indexOfInterest,2),'k')
hold on
plot(t1(indexOfInterest),x_tot(indexOfInterest,4),'b')
grid on
box on
zoom('on')



% ----------------- Functions -----------------

% Dynamics of the system with observer------------------------------------------
function [dx_dt_tot, u]=LTIdynamics_with_observer_finite(t,xtot,A,B,C,u,P0,Qtilde,Rtilde,T)
x=xtot(1:2);
Sigmax=chol(Qtilde);
Sigmay=chol(Rtilde);
vx=randn(1,2)*Sigmax;
vy=randn(1,1)*Sigmay;

dx_dt=A*x+B*u+vx';
y=C*x+vy';

if t>0
    tspan=[0 t];
    [~,Pvec] = ode45(@(t1,z)riccati_Kalman(t1,z,size(P0),A,C,Qtilde,Rtilde),tspan,P0);
    for i=size(Pvec,1)
        Pvec_mat_out=reshape(Pvec(i,:),size(P0));
        L_KF=Pvec_mat_out*C'*inv(Rtilde);
    end

else
    Pvec_mat_out{1}=P0;
    L_KF=Pvec_mat_out{1}*C'*inv(Rtilde);
end



xhat=xtot(3:4);
dxhat_dt_KF=A*xhat+B*u+L_KF*(y-C*xhat);
dx_dt_tot=[dx_dt; dxhat_dt_KF];
end


% Riccati equation for Kalman filter------------------------------------------
function dPvec_dt = riccati_Kalman(~,Pvec,sizeP,A,C,Qtilde,Rtilde)
Pvec_mat=reshape(Pvec,sizeP);
dPvec_dt_mat=A*Pvec_mat+Pvec_mat*A'+Qtilde-Pvec_mat*C'*inv(Rtilde)*C*Pvec_mat;
dPvec_dt=dPvec_dt_mat(:);
end

