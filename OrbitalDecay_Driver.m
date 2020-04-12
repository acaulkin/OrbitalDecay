%Andrew Caulkins
%Cubesat Deorbit
%__________________________________________________________________________
clc
clear
close all

%Given Parameters:
global R_e M_e G m d mu rho_0 lambda omega_0 I k h_min h_0

R_e = 6371000;          %Radius of Earth [m]
M_e = 5.976*10^24;      %Mass of Earth [kg]
G = 6.67*10^-11;        %Gravitational Contstant [(m^3)/(kg*s^2)]
m = 1.2;                %Cubesat Mass [kg]
d = 0.1;                %Cubesat Diameter [m]
mu = 1.5*10^-5;         %Air Dynamic Viscosity [N*s/m^2]
rho_0 = 1.225;          %Sea-Level Density [kg/m^3]
lambda = -0.00015;      %Density Lapse Rate [1/m]
h_min = 80000;          %Cubesat Deorbit Height [m]
omega_0 = 3.0;          %Cubesat Initial Pitch Rate [rad/s]
I = 0.004;              %Cubesat Moment of Inertia [kg*m^2]
k = [ 0 0 1];           %k - unit vector [unitless]
h_0 = 300*1000;       %Initial Height Above Earth [meters]

%Defining Initial State:
u_0 = [0, R_e + h_0,-sqrt(G*M_e/(R_e+h_0)), 0, 0, omega_0];
% T = sqrt((4*pi^2*(h_0+R_e))/(G*M_e));

%Determining the Different Time and State Values using Different Methods:
%  [tn,un] = satFE(u_0,100,50);             %Forward Euler
%  [tn2,un2] = satBDF2(u_0,100,50);            %BDF2
%  [tn3,un3] = satRK4(u_0,100,50);     %RK4 Method
% % [tn4,un4] = satBE(u_0,10,1.0);           %AJ Special
%Determining Forbenius Norm of the Error:
E  = [sat3f(0.01) sat3f(0.1)];
vpa(abs(log10((E(1))/E(2))),2)  %must be ~ 2.0

%Question 3 Plot:
[tnQuest3,unQuest3] = satRK4(u_0,1000,2000);

%Altitude is equal to the norm of the position vector - Radius of earth:
Alt = sqrt(unQuest3(:,1).^2 + unQuest3(:,2).^2);

%Plotting:
figure(1)
plot(tnQuest3,Alt-R_e,'b')
title('Altitude as a Function of Time for RK4 Method')
xlabel('Time [s]')
ylabel('Altitude [m]')


% % %Convergence Study [Question 4]:

%Determine EXACT position using small time-step RK4 Method:
[tnExact,unExact] = satRK4(u_0,100,6400);  %RK4 Method  *Question 4)
Exact = unExact(end,1:2);


%Different Time Steps for Forward Euler:
[tnEuler1,unEuler1] = satFE(u_0,100,50);
EulerFinalState1 = unEuler1(end,1:2);

[tnEuler2,unEuler2] = satFE(u_0,100,100);
EulerFinalState2 = unEuler2(end,1:2);

[tnEuler3,unEuler3] = satFE(u_0,100,200);
EulerFinalState3 = unEuler3(end,1:2);

[tnEuler4,unEuler4] = satFE(u_0,100,400);
EulerFinalState4 = unEuler4(end,1:2);

[tnEuler5,unEuler5] = satFE(u_0,100,800);
EulerFinalState5 = unEuler5(end,1:2);

[tnEuler6,unEuler6] = satFE(u_0,100,1600); 
EulerFinalState6 = unEuler6(end,1:2);
%Determing Euclidean Distance for Forward Euler:
FEDist1 = sqrt((EulerFinalState1(1) - Exact(1))^2 + (EulerFinalState1(2) - Exact(2))^2);
FEDist2 = sqrt((EulerFinalState2(1) - Exact(1))^2 + (EulerFinalState2(2) - Exact(2))^2);
FEDist3 = sqrt((EulerFinalState3(1) - Exact(1))^2 + (EulerFinalState3(2) - Exact(2))^2);
FEDist4 = sqrt((EulerFinalState4(1) - Exact(1))^2 + (EulerFinalState4(2) - Exact(2))^2);
FEDist5 = sqrt((EulerFinalState5(1) - Exact(1))^2 + (EulerFinalState5(2) - Exact(2))^2);
FEDist6 = sqrt((EulerFinalState6(1) - Exact(1))^2 + (EulerFinalState6(2) - Exact(2))^2);

%Vector Format:
FEDist = [FEDist1 FEDist2 FEDist3 FEDist4 FEDist5 FEDist6];

%Different Time-Steps for RK4
[tnRK4,unRK41] = satRK4(u_0,100,50);
RK4FinalState1 = unRK41(end,1:2);

[tnRK42,unRK42] = satRK4(u_0,100,100);
RK4FinalState2 = unRK42(end,1:2);

[tnRK43,unRK43] = satRK4(u_0,100,200);
RK4FinalState3 = unRK43(end,1:2);

[tnRK44,unRK44] = satRK4(u_0,100,400);
RK4FinalState4 = unRK44(end,1:2);

[tnRK45,unRK45] = satRK4(u_0,100,800);
RK4FinalState5 = unRK45(end,1:2);

[tnRK46,unRK46] = satRK4(u_0,100,1600); 
RK4FinalState6 = unRK46(end,1:2);
%Determing Euclidean Distance for RK4 Method:
RK4Dist1 = sqrt((RK4FinalState1(1) - Exact(1))^2 + (RK4FinalState1(2) - Exact(2))^2);
RK4Dist2 = sqrt((RK4FinalState2(1) - Exact(1))^2 + (RK4FinalState2(2) - Exact(2))^2);
RK4Dist3 = sqrt((RK4FinalState3(1) - Exact(1))^2 + (RK4FinalState3(2) - Exact(2))^2);
RK4Dist4 = sqrt((RK4FinalState4(1) - Exact(1))^2 + (RK4FinalState4(2) - Exact(2))^2);
RK4Dist5 = sqrt((RK4FinalState5(1) - Exact(1))^2 + (RK4FinalState5(2) - Exact(2))^2);
RK4Dist6 = sqrt((RK4FinalState6(1) - Exact(1))^2 + (RK4FinalState6(2) - Exact(2))^2);

%Vector Format:
RK4Dist = abs([RK4Dist1 RK4Dist2 RK4Dist3 RK4Dist4 RK4Dist5 RK4Dist6]);

%TimeSteps/etc...
N = [50 100 200 400 800 1600];
delta_T = [(100/50) (100/100) (100/200) (100/400) (100/800) (100/1600)];
% 
% %Different Time-Steps for BKFD2
% [tnBD,unBD1] = satBDF2(u_0,100,50);
% BDFinalState1 = unBD1(end,1:2);
% 
% [tnBD2,unBD2] = satBDF2(u_0,100,100);
% BDFinalState2 = unBD2(end,1:2);
% 
% [tnBD3,unBD3] = satBDF2(u_0,100,200);
% BDFinalState3 = unBD3(end,1:2);
% 
% [tnBD4,unBD4] = satBDF2(u_0,100,400);
% BDFinalState4 = unBD4(end,1:2);
% % 
% % [tnBD5,unBD5] = satBDF2(u_0,100,800);
% % BDFinalState5 = unBD5(end,1:2);
% % 
% % [tnBD6,unBD6] = satBDF2(u_0,100,1600); 
% % BDFinalState6 = unBD6(end,1:2);
% 
% %Determing Euclidean Distance for RK4 Method:
% BDDist1 = sqrt((BDFinalState1(1) - Exact(1))^2 + (BDFinalState1(2) - Exact(2))^2);
% BDDist2 = sqrt((BDFinalState2(1) - Exact(1))^2 + (BDFinalState2(2) - Exact(2))^2);
% BDDist3 = sqrt((BDFinalState3(1) - Exact(1))^2 + (BDFinalState3(2) - Exact(2))^2);
% BDDist4 = sqrt((BDFinalState4(1) - Exact(1))^2 + (BDFinalState4(2) - Exact(2))^2);
% % BDDist5 = sqrt((BDFinalState5(1) - Exact(1))^2 + (BDFinalState5(2) - Exact(2))^2);
% % BDDist6 = sqrt((BDFinalState6(1) - Exact(1))^2 + (BDFinalState6(2) - Exact(2))^2);
% 
% %Vector Format:
% % BDDist = abs([BDDist1 BDDist2 BDDist3 BDDist4 BDDist5 BDDist6]);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
%Slopes:
FESlope = diff(log(FEDist))./diff(log(delta_T));
RK4Slope = diff(log(RK4Dist))./diff(log(delta_T));
% BDSlope = diff(log(BDDist))./diff(log(delta_T2));

%Plotting Euclidean Distance vs. N:
figure(2)
loglog(delta_T,FEDist)
hold on
loglog(delta_T,RK4Dist)
% loglog(delta_T,BDDist)
title('Logarithmic Convergence Study')
xlabel('Number of Time Steps')
ylabel('Euclidean Distance')
legend('Forward Euleran ','RK4')
axis equal


%Question 5:
%The purpose of this question is to determine the point of satellite
%de-orbit (defined as 80,000 meters) using a time integration method of my
%choice.

%Using RK4:
[TimeRK4,StateRK4] = satRK4(u_0,100e3,(70e3));
PosRK4 = StateRK4(:,1:2);
ALT = zeros(size(PosRK4,1),1);
for O=1:size(PosRK4,1)
    ALT(O) = sqrt(PosRK4(O,1).^2 + PosRK4(O,2).^2) - R_e;
    if ALT(O) <= h_min
        TimeOfDeOrbit = TimeRK4(O);
        break;
    end
end



%Plotting Altitude vs time:
figure(2)
plot(TimeRK4,ALT)
title('Altitude as a Function of Time to determine De-Orbit')
xlabel('Time [s]')
ylabel('Altitude [m]')
legend('De-Orbit time: 22.6940 Hours')
axis equal



%Checking if time value is correct:
figure(3)
plot(PosRK4(:,1),PosRK4(:,2),'b')
hold on
plot(PosRK4(O,1),PosRK4(O,2),'r*')
title('Orbit of Satellite in Cartesian Space')
legend('Flight Path','Point of De-Orbit')
xlabel('x Axis [m]')
ylabel('y-Axis [m]')
axis equal


%Plotting Orbit:
% figure(2)
% plot(un(:,1),un(:,2),'r') %FE Method
% hold on
% plot(un2(:,1),un2(:,2),'b') %BFD2 Method
% plot(un3(:,1),un3(:,2),'g') %RK4 Method
% title('Plot of Orbit')
% grid minor
% legend('Forward Euler', 'BDF2','RK4')
% axis equal
function [tn,un] = satFE(u_0,T,N )
% PURPOSE:
% This function solves the satellite orbit problem using forward Euler
% INPUTS:
% u_0 : initial state [x, y, xdot, ydot, theta, thetadot]
% T : end time (start is 0) [ s ]
% N : Number of time-steps
% OUTPUTS:
% tn : (N+1)x1 vector of time values
% un : (N+1)x4 vector of states
%___________________________________________________


%Initialize Arrays:
dt = T/N; % Number of time-steps
tn = linspace(0,T,N+1);
un = zeros(N+1,6);

un(1,:) = u_0; % initial state


%Forward Euler Time Integration:
for i=1:N
    tn(i+1) = tn(i) + dt;
    un(i+1,:) = un(i,:) + dt*fesat(un(i,:));
end

end

%Function to compute the new state, f(u):
function f = fesat(u)
global lambda R_e k M_e m d mu rho_0 G I
 
Pos = [u(1) u(2) 0];                        %Position
V = [u(3) u(4) 0];                          %Velocity
rho = rho_0*exp(lambda*(norm(Pos) - R_e));     %Reynold's Number
r = norm(Pos);
Re = rho*norm(V)*d/mu;

%Accelerations:
%Drag Acceleration:
C_D = (0.1)/(1+Re) + 0.05*(1+cos(u(5)));
D = -((1/2)*(rho)*norm(V)*V*d^2*C_D)/m;
%Lift Acceleration:
C_L = (1.0)/(1+Re^0.25)*sin(u(5));
L = ((1/2)*rho*norm(V)*cross(V,k)*d^2*C_L)/m;
%Acceleration Due to Earth:
a_earth = [((-G*M_e*Pos(1))/(r)^3), ((-G*M_e*Pos(2))/(r)^3)];
%Angular Acceleration:
C_M = (0.01)/(1+Re^1.2)*u(6);
M = (1/2)*(rho)*norm(V)^2*d^3*k*C_M;  
    

%NEW STATE:
f(1) = u(3);                            %x_dot
f(2) = u(4);                            %y_dot
f(3) =  a_earth(1) + D(1) + L(1);       %x_ddot
f(4) = a_earth(2) + D(2) + L(2);        %y_ddot
f(5) = u(6);                            %theta_dot
f(6) = -dot(M,k)/I;                     %theta_ddot

end



%RK4 Method:
function [tn2,un2] = satRK4(u0,T,N)% function 
dt = T/N;     % Number of Time-Steps
tn2 = linspace(0,T,N+1);


un2(1,:) = u0; % initial state
% RK4 time integration
for n=1:N% RK4
    f0 = fesat(un2(n,:));
    f1 = fesat(un2(n,:) + 0.5*dt*f0);
    f2 = fesat(un2(n,:) + 0.5*dt*f1);
    f3 = fesat(un2(n,:) +     dt*f2);
    un2(n+1,:) = un2(n,:) + dt/6*(f0+2*f1+2*f2+f3);
end
end 





function err =  sat3f(epsilon)
R_e = 6371000;          %Radius of Earth [m]
r = R_e + 80000;
% initial condition
u0 = [0.6*r,0.8*r,-4000,5000,0.7,2.1];



% loop over columns: build analytical/numerical matrices
for j=1:6,
  uprime = u0;
  uprime(j) = uprime(j) + epsilon;
  A = 0.5*(Linear(u0) + Linear(uprime));
  adfdu(:,j) = A(:,j); % analytical
  ndfdu(:,j) = (fesat(uprime) - fesat(u0))/(epsilon); % numerical
end

% Frobenius norm of the error between the matrices
err = norm(adfdu-ndfdu, 'fro');

end


%AJ Special Method:
function [tn,un] = satBE(u0,T,N)% function 
dt = T/N;     % Number of Time-Steps
tn = linspace(0,T,N+1);
un = zeros(N+1,6);
un(1,:) = u0; % initial state
A = Linear(un(1,:));
I = eye(6);

 %Applying AJ Special:
 for n=1:N
     tn(n+1) = tn(n) + dt;
     un(n+1,:) = ((I - dt.*A)^-1) * un(n,:)';
 end
 end 



% function [t,u] = BDF1(u0,T,dt)
%  % PURPOSE:
%  % Performs trapezoidal integration of a nonlinear problem
%  % INPUTS:
%  % N : number of time steps
%  % OUTPUTS:
%  % t : vector of N+1 time node values
%  % u : 2 by N+1 solution vector
% 
%  N = T/dt; % time step
%  t = linspace(0,T,N+1); % vector of time nodes
%  u = zeros(N+1,6); % state vector
%  u(1,:) = u0; % initial condition
%  Rnorm = 1e-8; % Newton?Raphson residual tolerance
% 
%  for n = 1:N % loop over time steps
% 
%  fn = fesat(u (n,:)); % fvalue at u^n
%  w = u(n,:); % Newton?Raphson initialization
% 
%  for k=1:3 % allow max 10 Newton?Raphson iterations
%  R = w - u(n,:) - dt*(fesat(w) + fn);
%  if (norm(R) < Rnorm) break; end
%  R_w = eye(6) - dt.*Linear(w);
%  dw = -R_w \ R';
%  w = w + dw';
%  end
% 
%  % make sure the residual tolerance was met
%  if (norm(R) > Rnorm),error('residual norm not met!')
%  end
%  
%  u (n+1,:) = w; % store solution as next state
% 
%  end
% end

%Second Order Backward Difference:
function [t,u] = satBDF2(u0,T,N)

dt = T/N; % time step
t = linspace(0,T,N+1); % vector of time nodes
u = zeros(N+1,6); % state vector
u(1,:) = u0; % initial condition
Rnorm = 1e-8; % Newton?Raphson residual tolerance

 for n = 1:N % loop over time steps   
    if (n==1) 
        fn = fesat(u (n,:)); % F(u) value at u^n
        w = u(n,:); % Newton-Raphson initialization
        for k=1:10 % allow max 10 Newton-Raphson iterations
        R = w - u(n,:) - dt*(fesat(w) + fn);
        if (norm(R) < Rnorm) break; end
        R_w = eye(6) - dt*Linear(w);
        dw = -R_w \ R';
        w = w + dw';
        end

    % make sure the residual tolerance was met
        if (norm(R) > Rnorm),error('residual norm not met!'); end
         u (n+1,:) = w; % store solution as next state  
    else 
       fn = fesat(u (n,:)); % F(u) value at u^n
        w = u(n,:); % Newton-Raphson initialization

        
        for k=1:3 % allow max 10 Newton-Raphson iterations     
            R =(3/2)*( w - (4/3)*u(n,:) + (1/3)*u(n-1,:) - (2/3)*dt*(fesat(w) + fn));
            if (norm(R) < Rnorm) break; end
            R_w = (3/2)*(eye(6) - (2/3)*dt*Linear(w));
            dw = -R_w\R';
            w = w + dw';
        end
        % make sure the residual tolerance was met
     if (norm(R) > Rnorm),error('residual norm not met!'); end;
     u (n+1,:) = w; % store solution as next state  
    end
 end
end




