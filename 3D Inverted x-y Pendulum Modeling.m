
% 3D Inverted x-y Pendulum Modeling

% AmirKabir University Of Technology (Tehran Polytechnic) - 15/11/1402

% Hamid Reza Heidari 

clc;
clear;
close all;

M = 1; % Mass of cart
m = 0.1; % Mass of Pendulum ( assume at Center Of Gravity )
l = 0.3; % Pendulum length from the Centre Of Gravity
g = 9.81; % Gravity Coefficient


% --------------- System Modeling ---------------------------


syms x1 x2 x3 x4 x5 x6 x7 x8 Fx Fy 

A = M*m*l* (sin(x5*(x6^2))) + M*m*l * (cos(x5))^2 * (sin(x5*(x8^2))) - M*m*g * (cos(x5))^2 * sin(x5) * cos(x7)^2 ...
+ m * cos(x5)^2 * sin(x7)^2 * Fx + M * Fx - m * cos(x5)^2 * sin (x5) * sin (x7) * Fy;

B = (M + m) * Fy - m * cos (x5)^2 * Fy + m * sin (x5) * cos(x5) * sin (x7) * Fx - M*m*g * cos(x5)^2 * sin (x7) * cos (x7) ...
+ M*m*l * cos(x5)^3 * sin (x7*(x8^2)) - M*m*l * cos (x5) * sin (x7*(x6^2)) ;


C = m * cos(x5) * cos(x7)^2 * Fx - (M + m) * cos(x5) * Fx + (M + m) * sin(x5) * sin(x7) * Fy + M*(m + M)*g * sin(x5) * cos(x7) ...
- M*m*l * sin(x5) * cos (x5) * cos (x7*(x6^2))^2 - M*(m + M)*l * cos(x5) * sin(x7*(x8^2)) ;

D = -m * sin(x5)^2 * cos(x7) * Fy - M * cos(x7) * Fy + m * sin(x5) * cos(x5) * sin(x7) * cos(x7) * Fx + 2*M*m*l * sin(x5*x6*x8) + M*(m + M)*g * sin (x7) + 2*(M^2)*l * sin (x5*x6*x8) ...
- 2 * M*m*l * sin(x5) * cos(x5)^2 * cos(x7*x6*x8)^2 - M*m*l * cos(x5)^3 * sin (x7) * cos (x7*(x8^2)) - M*m*l * cos(x5) * sin(x7) * cos(x7*(x6^2));

E = M^2+M*m * sin(x5)^2 * cos(x7)^2 +M*m * sin (x7) ;

x1_dot = x2; 
x2_dot = A / E;
x3_dot = x4;
x4_dot = B / E;
x5_dot = x6;
x6_dot = C / (E*l);
x7_dot = x8;
x8_dot = D / (E*l*cos(x5));


% --------------- Jacobian Matrix ---------------------------

IP_Jacobian_A = jacobian([x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot;x7_dot;x8_dot],[ x1 ;x2 ;x3 ;x4 ;x5 ;x6 ;x7 ;x8])
 
IP_Jacobian_B = jacobian([x1_dot;x2_dot;x3_dot;x4_dot;x5_dot;x6_dot;x7_dot;x8_dot],[Fx ; Fy])


% --------------- Equilibrium Points ---------------------------

x1=0;
x2=0;
x3=0;
x4=0;
x5=0;
x6=0;
x7=0;
x8=0;

Fx=0;
Fy=0;


% --------------- linearize System ---------------------------

A = zeros(8,8);

B = zeros(8,2);



for i=1:8
    for j=1:8
        A(i,j) = subs(IP_Jacobian_A (i,j));

        j = j+1;
    end
    i = i+1;
end


for i=1:8
    for j=1:2
        B(i,j) = subs(IP_Jacobian_B (i,j));

        j = j+1;
    end
    i = i+1;
end

A , B 

C = [1 0 0 0 0 0 0 0 ;0 0 1 0 0 0 0 0 ; 0 0 0 0 1 0 0 0 ; 0 0 0 0 0 0 1 0]

D = zeros(4,2)

% --------------- Observability & Controllability Check ---------------------------


Observ = obsv(A,C);
Control = ctrb(A,B);

O_rank = rank(Observ)
C_rank = rank(Control)


