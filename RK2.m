function [ sol, t ] = RK2( A,B,X0,dt,t_initial,t_final,order )
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
    k2=(A*(X(1:order,m-1)+k1)+B)*dt;
    X(1:order,m)=X(1:order,m-1)+(k1+k2)*1/2;
    t(m)=t(m-1)+dt;
end
sol=[X];
end