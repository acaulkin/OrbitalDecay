%Andrew Caulkins
%Linearizaion
%__________________________________________________________________________




function A = Linear(U)

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
h_0 = 90*1000;          %Initial Height Above Earth [meters]







syms x y x_dot y_dot theta theta_dot
%Defining my state U:
U1 = [x y x_dot y_dot theta theta_dot];


Pos = [U1(1) U1(2) 0];                        %Position
V = [U1(3) U1(4) 0];                          %Velocity
rho = rho_0*exp(lambda*(norm(Pos) - R_e));     %Reynold's Number
r = norm(Pos);
Re = rho*norm(V)*d/mu;

%Drag Acceleration:
C_D = (0.1)/(1+Re) + 0.05*(1+cos(U1(5)));
D = -((1/2)*(rho)*norm(V)*V*d^2*C_D)/m;
%Lift Acceleration:
C_L = (1.0)/(1+Re^0.25)*sin(U1(5));
L = ((1/2)*rho*norm(V)*cross(V,k)*d^2*C_L)/m;
%Acceleration Due to Earth:
a_earth = [((-G*M_e*Pos(1))/(r)^3), ((-G*M_e*Pos(2))/(r)^3)];
%Angular Acceleration:
C_M = (0.01)/(1+Re^1.2)*U1(6);
M = (1/2)*(rho)*norm(V)^2*d^3*k*C_M;  


%NEW STATE:
F = [U1(3),U1(4),a_earth(1) + D(1) + L(1),a_earth(2) + D(2) + L(2),U1(6),-dot(M,k)/I] ;                



%Taking Derivatives WRT U values

%F1:
F1_U1 = diff(F(1),U1(1)); 
F1_U2 = diff(F(1),U1(2)); 
F1_U3 = diff(F(1),U1(3)); 
F1_U4 = diff(F(1),U1(4)); 
F1_U5 = diff(F(1),U1(5)); 
F1_U6 = diff(F(1),U1(6)); 

F1_Deriv = [F1_U1 F1_U2 F1_U3 F1_U4 F1_U5 F1_U6];


%F2:
F2_U1 = diff(F(2),U1(1)); 
F2_U2 = diff(F(2),U1(2)); 
F2_U3 = diff(F(2),U1(3)); 
F2_U4 = diff(F(2),U1(4)); 
F2_U5 = diff(F(2),U1(5)); 
F2_U6 = diff(F(2),U1(6)); 

F2_Deriv = [F2_U1 F2_U2 F2_U3 F2_U4 F2_U5 F2_U6];


%F3:
F3_U1 = diff(F(3),U1(1)); 
F3_U2 = diff(F(3),U1(2)); 
F3_U3 = diff(F(3),U1(3)); 
F3_U4 = diff(F(3),U1(4)); 
F3_U5 = diff(F(3),U1(5)); 
F3_U6 = diff(F(3),U1(6)); 

F3_Deriv = [F3_U1 F3_U2 F3_U3 F3_U4 F3_U5 F3_U6];


%F4:
F4_U1 = diff(F(4),U1(1)); 
F4_U2 = diff(F(4),U1(2)); 
F4_U3 = diff(F(4),U1(3)); 
F4_U4 = diff(F(4),U1(4)); 
F4_U5 = diff(F(4),U1(5)); 
F4_U6 = diff(F(4),U1(6)); 

F4_Deriv = [F4_U1 F4_U2 F4_U3 F4_U4 F4_U5 F4_U6];


%F5:
F5_U1 = diff(F(5),U1(1)); 
F5_U2 = diff(F(5),U1(2)); 
F5_U3 = diff(F(5),U1(3)); 
F5_U4 = diff(F(5),U1(4)); 
F5_U5 = diff(F(5),U1(5)); 
F5_U6 = diff(F(5),U1(6)); 

F5_Deriv = [F5_U1 F5_U2 F5_U3 F5_U4 F5_U5 F5_U6];


%F2:
F6_U1 = diff(F(6),U1(1)); 
F6_U2 = diff(F(6),U1(2)); 
F6_U3 = diff(F(6),U1(3)); 
F6_U4 = diff(F(6),U1(4)); 
F6_U5 = diff(F(6),U1(5)); 
F6_U6 = diff(F(6),U1(6)); 

F6_Deriv = [F6_U1 F6_U2 F6_U3 F6_U4 F6_U5 F6_U6];



%Linearization Matrix:
P = [F1_Deriv;F2_Deriv;F3_Deriv;F4_Deriv;F5_Deriv;F6_Deriv];

%Testing if Correct:
U_0 = [0, R_e + h_0,-sqrt(G*M_e/(R_e+h_0)), 0, 0, omega_0];

A = subs(P,U1,U);
A = double(A);


end

