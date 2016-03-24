function [ sol, t ] = RK5( A,B,X0,dt,t_initial,t_final,order )
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
%% Inputs
% A              : state matix of the diffrential equations  (nxn)
% B              : B matix of the diffrential equations (nx1)
% X0            : initial condition matix (nx1)
% dt             : time step
% t_initial    : initial time
% t_final      : final time
% order       : n
%% Outputs
% sol              : the solution of diffrential equation [x, dx/dt, d2x/dt2, ... , d(order)x/dt(order) ; t]
% t                  : solution time
%% Function body
n=(t_final-t_initial)/dt;
X(1:order,1)=X0;
t(1)=t_initial;
for m=2:ceil(n)+1
    k1=(A*X(1:order,m-1)+B)*dt;
    k2=(A*(X(1:order,m-1)+k1/4)+B)*dt;
    k3=(A*(X(1:order,m-1)+k1*3/32+k2*9/32)+B)*dt;
    k4=(A*(X(1:order,m-1)+k1*1932/2197+k2*-7200/2197+k3*7296/2197)+B)*dt;
    k5=(A*(X(1:order,m-1)+k1*439/216+k2*-8+k3*3680/513+k4*-845/4104)+B)*dt;
    k6=(A*(X(1:order,m-1)+k1*-8/27+k2*2+k3*-3544/2565+k4*1859/4104+k5*-11/40)+B)*dt;
    X(1:order,m)=X(1:order,m-1)+(k1*16/135+k2*0+k3*6656/12825+k4*28561/56430+k5*-9/50+k6*2/55);
    t(m)=t(m-1)+dt;
end
sol=[X];
end


