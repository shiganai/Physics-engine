clear all
format compact

%% prepare L, Fx, Fy
syms m1 real
syms g real
syms l1 real

syms tau1 real

syms th1_pre(t)
syms x1(t)
syms y1(t)

syms th1 dth1 ddth1 real

syms_replaced = [
    th1_pre, diff(th1_pre, t), diff(th1_pre, t, t), ...
    ];

syms_replacing = [
    th1, dth1, ddth1, ...
    ];

p1 = [ x1, y1 ] + l1 * [cos(th1_pre), sin(th1_pre)];
p1_G = [ x1, y1 ] + 1/2 * l1 * [cos(th1_pre), sin(th1_pre)];

p_G_whole = ( m1 * p1_G ) / ( m1 );

p1 = formula(p1);
p1_G = formula(p1_G);
p_G_whole = formula( p_G_whole );

v1_G = formula( diff(p1_G, t) );
v_G_whole = formula( diff(p_G_whole, t) );

T1 = 1/2 * m1 * (v1_G * v1_G') + 1/2 * (1/12 * m1 * l1^2) * diff(th1_pre, t)^2;
U1 = m1 * g * p1_G(2);
L1 = T1 - U1;

% define force at fixed point
Fx = -functionalDerivative( L1, x1 );
Fy = -functionalDerivative( L1, y1 );

%% adopt that [ x1, y1 ] = [ 0, 0 ]
% prepare restricted variable
syms_assumed = [ x1, y1 ];
syms_assuming = sym( [ 0, 0 ] );

% replace
L1 = subs( L1, syms_assumed, syms_assuming );
p1_G = subs( p1_G, syms_assumed, syms_assuming );
p_G_whole = subs( p_G_whole, syms_assumed, syms_assuming );
v1_G = subs( v1_G, syms_assumed, syms_assuming );
v_G_whole = subs( v_G_whole, syms_assumed, syms_assuming );
Fx = subs( Fx, syms_assumed, syms_assuming );
Fy = subs( Fy, syms_assumed, syms_assuming );

%% variable to be used for verification

% define Momentum
M1 = m1 * v1_G;
dM1 = formula( diff( M1, t ) );

unclear_outer_force = [ Fx, Fy ];

% define Angular Momentum from p_G_whole
AM1_G = formula( ...
    cross( ...
        [ p1_G, 0 ] - [p_G_whole, 0], ...
        m1 * ( [ v1_G, 0 ] - [ v_G_whole, 0 ] ) ...
    ) );
AM1_G = AM1_G(3);

AM1_each_rot = (1/12 * m1 * l1^2) * diff(th1_pre, t);
dAM1 = formula( diff( AM1_each_rot, t ) + diff( AM1_G, t ) );

unclear_outer_force_torque = ...
    cross( ...
        [p_G_whole, 0] - [p_G_whole, 0], ...
        [0, -m1 * g, 0] ...
    ) + ...
    cross( ...
        [0, 0, 0] - [p_G_whole, 0], ...
        [Fx, Fy, 0] ...
    );
unclear_outer_force_torque = formula( unclear_outer_force_torque );
unclear_outer_force_torque = unclear_outer_force_torque(3);

%% showing : functionalDerivative = dLdq - ddtdLddq for x1

dLdq_x1 = simplify( diff( L1, x1 ) );
ddtdLddq_x1 = simplify( diff( L1, diff( x1, t ), t ) );

simplify( functionalDerivative( L1, x1 ) - ( dLdq_x1 - ddtdLddq_x1 ) )

%% showing : functionalDerivative = dLdq - ddtdLddq for y1

dLdq_y1 = simplify( diff( L1, y1 ) );
ddtdLddq_y1 = simplify( diff( L1, diff( y1, t ), t ) );

% simplify( functionalDerivative( L1, y1 ) - ( dLdq_y1 - ddtdLddq_y1 ) ) % output no zero due to conj
simplify( subs( functionalDerivative( L1, y1 ) - ( dLdq_y1 - ddtdLddq_y1 ), syms_replaced, syms_replacing ) )


%% showing : functionalDerivative = dLdq - ddtdLddq for th1
% -functionalDerivative is used. See https://www.mns.kyutech.ac.jp/~okamoto/education/mechanicsII/lagrangeeq_text.pdf

dLdq_th1 = simplify( diff( L1, th1_pre ) );
ddtdLddq_th1 = simplify( diff( L1, diff( th1_pre, t ), t ) );

simplify( functionalDerivative( L1, th1_pre ) - ( dLdq_th1 - ddtdLddq_th1 ) )

%% verify Newton eq

eq = [
    -functionalDerivative( L1, th1_pre ) == tau1;
    ];

eq = subs( eq, syms_replaced, syms_replacing );
dM1 = subs( dM1, syms_replaced, syms_replacing );
dAM1 = subs( dAM1, syms_replaced, syms_replacing );
Fx = subs( Fx, syms_replaced, syms_replacing );
Fy = subs( Fy, syms_replaced, syms_replacing );
unclear_outer_force = subs( unclear_outer_force, syms_replaced, syms_replacing );
unclear_outer_force_torque = subs( unclear_outer_force_torque, syms_replaced, syms_replacing );

variables = [ ddth1 ];

[A, B] = equationsToMatrix( eq, variables );
% det(A) % verify det(A) ~= 0
X = inv(A) * B;

dM1 = subs( dM1, variables, X' );
unclear_outer_force = subs( unclear_outer_force, variables, X' );

% this shold be equal to [0, -g*m1]
simplify( dM1 - unclear_outer_force )

dAM1 = subs( dAM1, variables, X' );
unclear_outer_force_torque = subs( unclear_outer_force_torque, variables, X' );

% this shold be equal to tau1
simplify( dAM1 - unclear_outer_force_torque )

Fx = subs( Fx, variables, X' );
Fy = subs( Fy, variables, X' );
FxFy = simplify( [Fx, Fy] );
% display when tau1 = 0
simplify( subs( FxFy, tau1, 0 ) )






































