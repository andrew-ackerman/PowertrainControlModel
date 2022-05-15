%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Powertrain Control 
%
% Model of I.C. engine dynamics for idle speed control.
%set(gcf,'color','w');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load engine geometric parameters and constant inputs

Vd = 2.4e-3;   % Displacement (m^3)
Z = 4;         % Number of Cylinders
Vm = 5.8e-3;   % Intake Manifold Volume (m^3)
J = 0.0789;    % Mass moment of inertia

p_amb = 1.0121*1e5;
T_amb = 302;
R=288;
gam = 1.35;

P0 = 26431;   % Initial MAP
N0 = 828;     % Initial RPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters (from steady-state calibration)

a = [1.69e-07,-1.136e-06,6.89e-06];  % Throttle
si = 0.812;   yi = 0.0633;           % Volumetric efficiency
P_IMEP = [0.0220,-1.1649];           % IMEP
par_sp = [-0.0017 -0.0277 1.36];     % Spark timing effect
par_fr = [7.4198e-7 -4.989e-4 11.3]; % Friction
par_IMEP0 = [1.2323e-4 2.1256];      % Base IMEP model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Conversion to Crank Angle Domain

% Coefficients for crank-angle based model
Kthr = p_amb/sqrt(R*T_amb)*sqrt(gam)*sqrt((2/(gam+1))^((gam+1)/(gam-1)));
Kp1 = R*T_amb/Vm;
Kp2 = si*Vd/(4*pi*Vm);
Kp3 = yi*Vd/(4*pi*Vm);
KT = 1e5*Vd/(4*pi);
Kfr1 = (30/pi)^2 * par_fr(1);
Kfr2 = (30/pi) * par_fr(2);
Kfr3 = par_fr(3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Equilibrium Condition (p_0 T_0 w_0)

setp(1) = 9.81; % Throttle Position
setp(2) = -25;  % Spark Timing
setp(3) = 10;   % Load Torque
X0 = [26424 21.3765773202354 83.9019428270409]; % Equilibrium Conditions

% New initial conditions for the NL model based on linearization
P0 = 26424;   % Initial MAP
N0 = 83.9019428270409*30/pi;     % Initial RPM

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linearization

% Coefficients for linearized model as shown in lecture
K1 = Kp1*Kthr*(2*a(1)*setp(1)+a(2));
K2 = Kp1*Kthr*(a(1)*setp(1)^2 + a(2)*setp(1) + a(3));
K3 = Kp2;
Kp = KT*par_IMEP0(1)*(par_sp(1)*setp(2)^2 + par_sp(2)*setp(2) + par_sp(3));    % Pressure multiplier
Kt = KT*(par_IMEP0(1)*X0(1) - par_IMEP0(2)) * (par_sp(1)*setp(2) + par_sp(2)); % Spark Timing multiplier
Kf = 2*Kfr1*X0(3)^2 + Kfr2*X0(3);
Ktheta = KT*(par_IMEP0(1)*X0(1)-par_IMEP0(2))*(2*par_sp(1)*setp(2)+par_sp(2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing the transfer function, poles and zeros of the system

Tau0 = pi/3; % Step size
CT1 = 1 + (0.5*K3*Tau0);
CT2 = -1 + (0.5*K3*Tau0);
CT3 = 0.5*K2*Tau0/(X0(3)^2);
CT4 = 0.5*K1*Tau0/(X0(3));
CT5 = 1 + ((0.5*Kf*Tau0)/(J*(X0(3)^2)));
CT6 = -1 + ((0.5*Kf*Tau0)/(J*(X0(3)^2)));
CT7 = 0.5*Tau0/(J*X0(3));


%% Computing the transfer function coefficients
num = Kp*CT4*CT7*[1 2 1];
den = [CT1*CT5 (CT1*CT6+CT2*CT5) CT2*CT6 Kp*CT3*CT7*[1 2 1]];
H_tf = tf(num,den,Tau0);

% Calculating poles and zeros
Poles = roots(den);
Zeros = roots(num);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices for the state-space representation SISO

M = 1/sqrt(Tau0) * [-CT1, -CT3; 0,  -CT5];
N = 1/sqrt(Tau0) * [ CT4, 0;    0,  CT7];
O = 1/sqrt(Tau0) * [ CT2, CT3;  0,  CT6];
P = 1/sqrt(Tau0) * [-CT4, 0;    0, -CT7];

Q = -O*M^(-1)*N+P;
R2 = O*M^(-1);
S = M^(-1);
U = -(M^(-1)*N);

Phi = [ R2(1,1), R2(1,2), 0, 0, Q(1,2)*Kp;...
        R2(2,1), R2(2,2), 0, 0, Q(2,2)*Kp;...
        S(1,1), S(1,2), 0, 0, U(1,2)*Kp;...
             0,      0, 1, 0,         0;...
             0,      0, 0, 1,         0 ]; 

Gamma =  [Q(1,1); Q(2,1); U(1,1); 0; 0];

H     =  [S(2,1), S(2,2), 0, 0, U(2,2)*Kp];       

D     =  [U(2,1)];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Matrices for the state-space representation MIMO
% Phi and H are same as the SISO case  

Gamma_m =  [Q(1,1), Q(1,2)*Ktheta; Q(2,1), Q(2,2)*Ktheta; U(1,1), U(1,2)*Ktheta; 0, 0; 0, 0];

D_m     =  [U(2,1), U(2,2)*Ktheta];

d = 3;

%% Problem 1: SISO Pole Placement

Phi_aS = [ R2(1,1), R2(1,2), 0, 0, Q(1,2)*Kp,0;...
        R2(2,1), R2(2,2), 0, 0, Q(2,2)*Kp,0;...
        S(1,1), S(1,2), 0, 0, U(1,2)*Kp,0;...
             0,      0, 1, 0,         0,0;...
             0,      0, 0, 1,         0,0;...
        S(2,1), S(2,2), 0, 0, U(2,2)*Kp,1 ]; 

Gamma_aS =  [Q(1,1); Q(2,1); U(1,1); 0; 0;0];

H_aS     =  [S(2,1), S(2,2), 0, 0, U(2,2)*Kp,0];       

D_aS     =  [U(2,1)];


p = [0.9+0.1i 0.9-0.1i 0.2+0.4i 0.2-0.4i .9,0.85];
K = place(Phi_aS,Gamma_aS,p);

sim('SISO_Model_P1');
dwplot = dw*9.549297+800; %Converts dw (rad/s) into rpm
figure(1)
subplot(3,1,1);
sgtitle('SISO Model Simulation with Integration')
stairs(theta,dwplot);
title('Engine Speed with Controller')
xlabel('Crank Angle (Radians)')
ylabel('Engine Speed (RPM)')
subplot(3,1,2);
alphaplot= dalpha+9.81;
stairs(theta,alphaplot);
title('Throttle Angle with Controller')
xlabel('Crank Angle (Radians)')
ylabel('Throttle Angle (deg)')
subplot(3,1,3)
stairs(theta,dTl);
title('Applied Torque Disturbance');
xlabel('Crank Angle (Radians)');
ylabel('Torque (Nm)');


%% Problem 2: MIMO

Phi_aM = [ R2(1,1), R2(1,2), 0, 0, Q(1,2)*Kp,0;...
        R2(2,1), R2(2,2), 0, 0, Q(2,2)*Kp,0;...
        S(1,1), S(1,2), 0, 0, U(1,2)*Kp,0;...
             0,      0, 1, 0,         0,0;...
             0,      0, 0, 1,         0,0;...
        S(2,1), S(2,2), 0, 0, U(2,2)*Kp,1 ]; 

Gamma_aM =  [Q(1,1), Q(1,2)*Ktheta; Q(2,1), Q(2,2)*Ktheta; U(1,1), U(1,2)*Ktheta; 0, 0; 0, 0;0,0];

H_aM     =  [S(2,1), S(2,2), 0, 0, U(2,2)*Kp,0,0];       

D_aM     =  [U(2,1), U(2,2)*Ktheta];

p = [0.7+0.05i 0.7-0.05i 0.05+0.1i 0.05-0.1i 0.5 -0.9];
K = place(Phi_aM,Gamma_aM,p);

sim('MIMO_Model_Prob2');
dwplot = dw*9.549297+800; %Converts dw (rad/s) into rpm
figure(2)
subplot(4,1,1);
sgtitle('MIMO Model Simulation with Integration')
stairs(theta,dwplot);
title('Engine Speed with Controller')
xlabel('Crank Angle (Radians)')
ylabel('Engine Speed (RPM)')
subplot(4,1,2);
alphaplot= dalpha+9.81;
stairs(theta,alphaplot);
title('Throttle Angle with Controller')
xlabel('Crank Angle (Radians)')
ylabel('Throttle Angle (deg)')
subplot(4,1,3);
sparkplot=25-spark;
stairs(theta,sparkplot);
title('Spark Advance Angle with Controller')
xlabel('Crank Angle (Radians)')
ylabel('Spark Angle (deg bTDC)')
subplot(4,1,4)
stairs(theta,dTl);
title('Applied Torque Disturbance');
xlabel('Crank Angle (Radians)');
ylabel('Torque (Nm)');