function [xdot, U, Forces, Moments] = remus100_withM(x,ui)

% function returns:
% xdot - time derivative of state vector
% speed U in m/s
% forces X,Y,Z,K,M,N
% function inputs:
% x - state vector 
% x = [ u v w p q r x y z phi theta psi ]' Euler angles
% or x = [ u v w p q r x y z eta eps1 eps2 eps3 ]' Unit quaternion
%   u:       surge velocity          (m/s)
%   v:       sway velocity           (m/s)
%   w:       heave velocity          (m/s)
%   p:       roll rate               (rad/s)
%   q:       pitch rate              (rad/s)
%   r:       yaw rate                (rad/s)
%   x:       North position          (m)
%   y:       East position           (m)
%   z:       downwards position      (m)
%   phi:     roll angle              (rad)       
%   theta:   pitch angle             (rad)
%   psi:     yaw angle               (rad)
% ui - control inputs: one tail rudder, two stern planes and a single-screw propeller
% ui = [ delta_r delta_s n ]'
% delta_r:   rudder angle (rad)
% delta_s:   stern plane angle (rad) 
% n:         propeller revolution (rpm)

% the length of the Remus 100 AUV is 1.6 m
% the cylinder diameter is 19 cm  
% the mass of the vehicle is 31.9 kg
% the maximum speed of 2.5 m/s is obtained when the propeller runs at 1525 rpm in zero currents.

if (nargin == 2), Vc = 0; betaVc = 0; w_c = 0; end  % no ocean currents
if (nargin == 4), w_c = 0; end             % no vertical ocean currents

if (length(ui) ~= 3),error('u-vector must have dimension 3!'); end
if (length(x) ~= 12 && length(x) ~= 13)
    error('x-vector must have dimension 12 or 13'); 
end
% REMUS Hydrodynamic Coefficients
% Daniel Sgarioto, DTA
% Nov 2006
global V scale
% Vehicle Parameters
%
% State vectors and control inputs
nu = x(1:6);            % velocities u,v,w,p,q,r
eta = x(7:12);          % positions x,y,z,phi,theta,psi
delta_r = ui(1);        % tail rudder (rad)
delta_s = ui(2);        % stern plane (rad)
n = ui(3)/60;           % propeller revolution (rps)
u = nu(1);
v = nu(2);
w = nu(3);

% AUV model parameters; Fossen (2021, Section 8.4.2) and Allen et al. (2000)
L_auv = 1.6;             % AUV length (m)
D_auv = 0.19;            % AUV diamater (m)
S = 0.7 * L_auv * D_auv; % planform area S = 70% of rectangle L_auv * D_auv
a = L_auv/2;             % spheroid semi-axes a and b
b = D_auv/2;                  
r44 = 0.3;               % added moment of inertia in roll: A44 = r44 * Ix
r_bg = [ 0 0 0.02 ]';    % CG w.r.t. to the CO
r_bb = [ 0 0 0 ]';       % CB w.r.t. to the CO

% Tail rudder (single)
CL_delta_r = 0.5;        % rudder lift coefficient (-)
A_r = 2 * 0.10 * 0.05;   % rudder area (m2)
x_r = -a;                % rudder x-position (m)

% Stern plane (double)
CL_delta_s = 0.7;        % stern-plane lift coefficient (-)
A_s = 2 * 0.10 * 0.05;   % stern-plane area (m2)
x_s = -a;                % stern-plane z-position (m)

% Low-speed linear damping matrix parameters
T1 = 20;                 % time constant in surge (s)
T2 = 20;                 % time constant in sway (s)
zeta4 = 0.3;             % relative damping ratio in roll
zeta5 = 0.8;             % relative damping ratio in pitch
T6 = 5;                  % time constant in yaw (s)

U0=V;
m = 30.48;
g = 9.81;
%

%
zg = 0.0196;
%
Ixx = 0.177;
Iyy = 3.45;
Izz = 3.45;
%
cdu = 0.2;
rho = 1030;
Af = 0.0285;
d = 0.191;
xcp = 0.321;
Cydb = 1.2;
%
mq = 0.3;
%
cL_alpha = 3.12;
Sfin = 0.00665;
xfin = -0.6827;
%
gamma = 1;
a_prop = 0.25;
w_prop = 0.2;
tau = 0.1;
Jm = 1;
%
l_prop = 0.8*0.0254;
d_prop = 5.5*0.0254;
A_prop = (pi/4)*(d_prop^2);
m_f = gamma*rho*A_prop*l_prop;
%
Kn = 0.5;
%
% Most are Prestero's estimates, but revised values of some linear
% coefficients are due to Fodrea. Thruster coeffs based on results 
% reported by Allen et al.
%
Xwq= -35.5;
Xqq= -1.93;
Xvr= 35.5;
Xrr= -1.93;
Yvv= -1310.0;
Yrr= 0.632;
Yuv= -28.6;
Yur= 5.22;
Ywp= 35.5;
Ypq= 1.93;
Yuudr= 9.64;
Zww= -1310.0;
Zqq= -0.632;
Zuw= -28.6;
Zuq= -5.22;
Zvp= -35.5;
Zrp= 1.93;
Zuuds= -9.64;
Kpp= -0.130;
Mww= 3.18;
Mqq= -188;
Muw= 24.0;
Muq= -2.0;
Mvp= -1.93;
Mrp= 4.86;
Muuds= -6.15;
Nvv= -3.18;
Nrr= -94.0;
Nuv= -24.0;
Nur= -2.0;
Nwp= -1.93;
Npq= -4.86;
Nuudr= -6.15;
Kpdot= -0.0704;
%
Xuu = -0.5*rho*cdu*Af;
Xu = -rho*cdu*Af*U0;
%
% Added Mass Coeffs
Xudot= -0.93;
%
Yvdot= -35.5;
Yrdot= 1.93;
%
Zwdot= -35.5;
Zqdot= -1.93;
%
Mwdot= -1.93;
Mqdot= -4.88;
%
Nvdot= 1.93;
Nrdot= -4.88;
%
% Added Mass Terms
%
Zwc = -15.7;
Zqc = 0.12;
%
Mwc = -0.403;
Mqc = -2.16;
%
% Added Mass Coupling Terms
Xqa = Zqdot*mq;
Zqa = -Xudot*U0;
Mwa = -(Zwdot - Xudot)*U0;
Mqa = -Zqdot*U0;
%
% Body Lift Contribution
%
Zwl = -0.5*rho*(d^2)*Cydb*U0;
%
Mwl = -0.5*rho*(d^2)*Cydb*xcp*U0;
%
% Fin Contribution
%
Zwf = -0.5*rho*cL_alpha*Sfin*U0;
Zqf = 0.5*rho*cL_alpha*Sfin*xfin*U0;
%
Mwf = 0.5*rho*cL_alpha*Sfin*xfin*U0;
Mqf = -0.5*rho*cL_alpha*Sfin*(xfin^2)*U0;
%
% Dive Plane Coeffs
Zw = Zwc + Zwl + Zwf;
Zq = 2.2;
%
Mw = -9.3;
Mq = Mqc +Mqa +Mqf;
%
% Steering Coeffs
Yv = Zw;
Yr = 2.2;
%
Nv = -4.47;
Nr = Mq;
%
% Control Surface Coeffs
Zds = -rho*cL_alpha*Sfin*(U0^2);
Mds = rho*cL_alpha*Sfin*xfin*(U0^2);
%
Ydr = -Zds/3.5;
Ndr = Mds/3.5;
%
% Thruster Coeffs
Tnn = 6.279e-004;
Tnu = 0; 
%
Qnn = -1.121e-005;
Qnu = 0;
%
Tnu0 = (1/(a_prop + 1))*Tnu;
Qnu0 = (1/(a_prop + 1))*Qnu;
%
df0 = (-Xu)/((1 - tau)*(1 + a_prop)*(1 - w_prop));
df = (-Xuu)/((1 - tau)*(1 + a_prop)*a_prop*((1 - w_prop)^2));

% Set total forces from equations of motion
% ----------------------------------------------------- -------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nprop = ui(3)/60;

%%
% Hydrodynamic Thrust & Torque
%
T = Tnn*abs(nprop)*nprop;
Q = Qnn*abs(nprop)*nprop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = x(4);
q = x(5);
r = x(6);
phi = eta(4);
theta = eta(5);
m=30.48;
g=9.81;
W = m*g;
B = W + (0.75*4.44822162);
L = 1.3327;

% Constants
mu = 63.446827;         % Lattitude for Trondheim, Norway (deg)
g_mu = gravity(mu);     % gravity vector (m/s2)



% Amplitude saturation of rudder angle, stern plane and propeller revolution
n_max = 1525;                                % maximum propeller rpm
max_ui = [deg2rad(30) deg2rad(30) n_max/60]';   % rad, rad, rps

if (abs(delta_r) > max_ui(1)), delta_r = sign(delta_r) * max_ui(1); end
if (abs(delta_s) > max_ui(2)), delta_s = sign(delta_s) * max_ui(2); end
if (abs(n)       > max_ui(3)), n = sign(n) * max_ui(3); end

% Angle of attack and vehicle speed
alpha = atan2(nu(3), nu(1));                  % angle of attack (rad) = taninv(w/u)
%alpha = -0.58*180/pi;
U  = sqrt( nu(1)^2 + nu(2)^2 + nu(3)^2 );     % speed (m/s) = sqrt(u^2+v^2+w^2)
beta = asin([nu(2),U]);                        % sideslip angle

% AUV model parameters; Fossen (2021, Section 8.4.2) and Allen et al. (2000)
L_auv = 1.6;             % AUV length (m)
D_auv = 0.19;            % AUV diamater (m)
S = 0.7 * L_auv * D_auv; % planform area S = 70% of rectangle L_auv * D_auv
a = L_auv/2;             % spheroid semi-axes a and b
b = D_auv/2;                  
r44 = 0.3;               % added moment of inertia in roll: A44 = r44 * Ix
r_bg = [ 0 0 0.02 ]';    % CG w.r.t. to the CO
r_bb = [ 0 0 0 ]';       % CB w.r.t. to the CO

% Parasitic drag coefficient CD_0, i.e. zero lift and alpha = 0
% F_drag = 0.5 * rho * Cd * (pi * b^2)   
% F_drag = 0.5 * rho * CD_0 * S
Cd = 0.42;                              % from Allen et al. (2000)
CD_0 = Cd * pi * b^2 / S;

% Propeller coeffs. KT and KQ are computed as a function of advance no.
% Ja = Va/(n*D_prop) where Va = (1-w)*U = 0.944 * U; Allen et al. (2000)
D_prop = 0.14;   % propeller diameter corresponding to 5.5 inches
t_prop = 0.1;    % thrust deduction number
Va = 0.944 * U;  % advance speed (m/s)

% Ja_max = 0.944 * 2.5 / (0.14 * 1525/60) = 0.6632
Ja_max = 0.6632;
        
% Single-screw propeller with 3 blades and blade-area ratio = 0.718.    
% >> [KT_0, KQ_0] = wageningen(0,1,0.718,3)
KT_0 = 0.4566;
KQ_0 = 0.0700;
% >> [KT_max, KQ_max] = wageningen(0.6632,1,0.718,3) 
KT_max = 0.1798;
KQ_max = 0.0312;
        
% Propeller thrust and propeller-induced roll moment
% Linear approximations for positive Ja values
% KT ~= KT_0 + (KT_max-KT_0)/Ja_max * Ja   
% KQ ~= KQ_0 + (KQ_max-KQ_0)/Ja_max * Ja  
      
if n > 0   % forward thrust
    % Force in x dirn on propeller
    X_prop = rho * D_prop^4 * (... 
        KT_0 * abs(n) * n + (KT_max-KT_0)/Ja_max * (Va/D_prop) * abs(n) );
    % Moment in roll dirn of propeller
    K_prop = rho * D_prop^5 * (...
        KQ_0 * abs(n) * n + (KQ_max-KQ_0)/Ja_max * (Va/D_prop) * abs(n) );           
            
else    % reverse thrust (braking)
        
    X_prop = rho * D_prop^4 * KT_0 * abs(n) * n; 
    K_prop = rho * D_prop^5 * KQ_0 * abs(n) * n;
            
end            


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mass Matrix
M = zeros(6,6);
M(1,1) = m - Xudot;
M(1,5) = m*zg;
M(1,6) = 0;
%
M(2,2) = m - Yvdot;
M(2,4) = -m*zg;
M(2,6) = - Yrdot;
%
M(3,3) = m - Zwdot;
M(3,4) = 0;
M(3,5) = - Zqdot;
%
M(4,2) = -m*zg;
M(4,3) = 0;
M(4,4) = Ixx - Kpdot;
%
M(5,1) = m*zg;
M(5,3) = - Mwdot;
M(5,5) = Iyy - Mqdot;
%
M(6,1) = 0;
M(6,2) = - Nvdot;
M(6,6) = Izz - Nrdot;
%
Minv = inv(M);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

C = m2c(M,nu);


% Dissipative forces and moments
D = Dmtrx([T1 T2 T6],[zeta4 zeta5],M,[W r_bg' r_bb']);

% 6x6 diagonal linear damping matrix
D(1,1) = D(1,1) * exp(-3*U);   % vanish at high speed where quadratic
D(2,2) = D(2,2) * exp(-3*U);   % drag and lift forces dominates
D(6,6) = D(6,6) * exp(-3*U);

tau_liftdrag = forceLiftDrag(D_auv,S,CD_0,alpha,U); % g0
tau_crossflow = crossFlowDrag(L_auv,D_auv,D_auv,nu); % w

% Kinematics
if (length(x) == 13)
    [J,R] = quatern(x(10:13)); % J(eta)
else
    [J,R] = eulerang(x(10),x(11),x(12)); % Transformation/Rotation matrices
end

% Restoring forces and moments 6x1 vector about CO 
% for a submerged body using the rotation matrix R as input
% g = gRvect(W,B,R,r_bg,r_bb); % g(eta)

% pitch and roll angles 
% theta = deg2rad(10);
% phi = deg2rad(30); 
% % g-vector 
% g = gvect(W,B,theta,phi,r_bg,r_bb); 

g = gRvect(W,B,R,r_bg,r_bb);
g = [...
   -(W-B) * R(3,1)
   -(W-B) * R(3,2)
   -(W-B) * R(3,3)
   -(r_bg(2)*W - r_bb(2)*B) * R(3,3) + (r_bg(3)*W - r_bb(3)*B) * R(3,2)
   -(r_bg(3)*W - r_bb(3)*B) * R(3,1) + (r_bg(1)*W - r_bb(1)*B) * R(3,3)
   -(r_bg(1)*W - r_bb(1)*B) * R(3,2) + (r_bg(2)*W - r_bb(2)*B) * R(3,1) ];

% Horizontal- and vertical-plane relative speed
U_rh = sqrt( nu(1)^2 + nu(2)^2 );  %u, v 
U_rv = sqrt( nu(1)^2 + nu(3)^2 );  %u, w

% Rudder and stern-plane drag
X_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r^2; 
X_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s^2;

% Rudder sway force 
Y_r = -0.5 * rho * U_rh^2 * A_r * CL_delta_r * delta_r;

% Stern-plane heave force
Z_s = -0.5 * rho * U_rv^2 * A_s * CL_delta_s * delta_s;

% Generalized propulsion force vector
tau = zeros(6,1);                                
tau(1) = (1-t_prop) * X_prop + X_r + X_s; % X
tau(2) = Y_r; % Y
tau(3) = Z_s; % Z
tau(4) = K_prop; % K
tau(5) = x_s * Z_s; % L = stern-plane z-position * stern-plane heave force
tau(6) = x_r * Y_r; % M = rudder x-position * rudder sway force 

% State-space model
xdot = [ M \ (tau + tau_liftdrag + tau_crossflow - C * nu - D * nu  - g)
         J * nu ]; 
% nu dot, eta dot
% velocities in body frame transformed to velocities in ned frame

%%%%%%%%%%%%%%%%%%% Forces: %%%%%%%%%%%%%%%%%%%%
tau=0.1;
%  X_ = (B-W)*sin(theta) + Xuu*u*abs(u) + (1 - tau)*Tnn*n^2;
  X_ = -(W-B)*sin(theta) + Xuu*u*abs(u) + (Xwq-m)*w*q + (Xqq)*q^2 ...
  + (Xvr+m)*v*r + (Xrr)*r^2 - m*zg*p*r + (1 - tau)*T;

% Y_ = (W-B)*cos(theta)*sin(phi) + Yuudr*u^2*delta_r ;
Y_ = (W-B)*cos(theta)*sin(phi) + Yvv*v*abs(v) + Yrr*r*abs(r) + Yuv*u*v...
+ (Ywp+m)*w*p + (Yur-m)*u*r - (m*zg)*q*r + (Ypq)*p*q ...
+ Yuudr*u^2*delta_r ;

% Z_ = (W-B)*cos(theta)*cos(phi) + Zww*w*abs(w) + Zuuds*u^2*delta_s ;
Z_ = (W-B)*cos(theta)*cos(phi) + Zww*w*abs(w) + Zqq*q*abs(q)+ Zuw*u*w ...
+ (Zuq+m)*u*q + (Zvp-m)*v*p + (m*zg)*p^2 + (m*zg)*q^2 ...
+ (Zrp)*r*p + Zuuds*u^2*delta_s ;

% K_ = -(zg*W)*cos(theta)*sin(phi) + Qnn*n*abs(n);
% 
% M_ = -(zg*W)*sin(theta) + Mww*w*abs(w) + Muw*u*w + Muuds*u^2*delta_s ;
% 
% N_ = Nvv*v*abs(v) + Nuv*u*v + Nuudr*u^2*delta_r ;

K_ = - (zg*W)*cos(theta)*sin(phi) ...
+ Kpp*p*abs(p) - (Izz-Iyy)*q*r - (m*zg)*w*p + (m*zg)*u*r + Q;

M_ = -(zg*W)*sin(theta) + Mww*w*abs(w) ...
+ Mqq*q*abs(q) + (Mrp - (Ixx-Izz))*r*p + (m*zg)*v*r - (m*zg)*w*q ...
+ (Muq)*u*q + Muw*u*w + (Mvp)*v*p ...
+ Muuds*u^2*delta_s ;

N_ = Nvv*v*abs(v) + Nrr*r*abs(r) + Nuv*u*v ...
+ (Npq - (Iyy-Ixx))*p*q + (Nwp)*w*p + (Nur)*u*r ...
+ Nuudr*u^2*delta_r ;

Forces = [X_ Y_ Z_]';
Moments = [K_ M_ N_]';

