clear all
format compact

%% prepare L
syms m1 real
syms l1 real
syms m2 real
syms l2 real

syms g real

syms tau1 tau2 real

syms th1_pre(t)
syms x1(t)
syms y1(t)

syms th2_pre(t)
syms x2(t)
syms y2(t)

syms th1 dth1 ddth1 real
syms th2 dth2 ddth2 real

syms_replaced = [
    th1_pre, diff(th1_pre, t), diff(th1_pre, t, t), ...
    th2_pre, diff(th2_pre, t), diff(th2_pre, t, t), ...
    ];

syms_replacing = [
    th1, dth1, ddth1, ...
    th2, dth2, ddth2, ...
    ];

% define postions
p1_top = [ x1, y1 ];
p1_bottom = p1_top + l1 * [cos(th1_pre), sin(th1_pre)];
p1_G = p1_top + 1/2 * l1 * [cos(th1_pre), sin(th1_pre)];
p2_top = [ x2, y2 ];
p2_bottom = p2_top + l2 * [cos(th2_pre), sin(th2_pre)];
p2_G = p2_top + 1/2 * l2 * [cos(th2_pre), sin(th2_pre)];

p_G_whole = ( m1 * p1_G + m2 * p2_G ) / ( m1 + m2 );

p1_top = formula(p1_top);
p1_bottom = formula(p1_bottom);
p1_G = formula(p1_G);
p2_top = formula(p2_top);
p2_bottom = formula(p2_bottom);
p2_G = formula(p2_G);
p_G_whole = formula( p_G_whole );

% define velocities
v1_G = formula( diff(p1_G, t) );
v2_G = formula( diff(p2_G, t) );
v_G_whole = formula( diff(p_G_whole, t) );

% define Inertia
I1 = (1/12 * m1 * l1^2);
I2 = (1/12 * m1 * l2^2);

% define L1, L2, L
T1 = 1/2 * m1 * (v1_G * v1_G') ...
    + 1/2 * I1 * diff(th1_pre, t)^2 ...
    ;
U1 = m1 * g * p1_G(2) ...
    ;
L1 = T1 - U1;

T2 = 1/2 * m2 * (v2_G * v2_G') ...
    + 1/2 * I2 * diff(th2_pre, t)^2 ...
    ;
U2 = m2 * g * p2_G(2) ...
    ;
L2 = T2 - U2;

L = L1 + L2;

%% solve and adopt restriction: [x2, y2] = p1_bottom
% MATLAB does not support x2( x1(t) ).
% meaning diff( x2, x1 ) cannot be calced.
% Thus, got to solve the restriction from bottom

% define force at p1_bottom ( = p2_top )
F2x = -functionalDerivative( L, x2 );
F2y = -functionalDerivative( L, y2 );

syms_assumed = [ x2, y2 ];
syms_assuming = p1_bottom;

% replace
L = subs( L, syms_assumed, syms_assuming );
p1_G = subs( p1_G, syms_assumed, syms_assuming );
p2_G = subs( p2_G, syms_assumed, syms_assuming );
p_G_whole = subs( p_G_whole, syms_assumed, syms_assuming );
v1_G = subs( v1_G, syms_assumed, syms_assuming );
v2_G = subs( v2_G, syms_assumed, syms_assuming );
v_G_whole = subs( v_G_whole, syms_assumed, syms_assuming );
F2x = subs( F2x, syms_assumed, syms_assuming );
F2y = subs( F2y, syms_assumed, syms_assuming );

%% solve and adopt restriction: [x1, y1] = [0, 0]
% MATLAB does not support x2( x1(t) ).
% meaning diff( x2, x1 ) cannot be calced.
% Thus, got to solve the restriction from bottom

% define force at p1_top ( = [0, 0] )
F1x = -functionalDerivative( L, x1 );
F1y = -functionalDerivative( L, y1 );

syms_assumed = [ x1, y1 ];
syms_assuming = sym( [0, 0] );

% replace
L = subs( L, syms_assumed, syms_assuming );
p1_G = subs( p1_G, syms_assumed, syms_assuming );
p2_G = subs( p2_G, syms_assumed, syms_assuming );
p_G_whole = subs( p_G_whole, syms_assumed, syms_assuming );
v1_G = subs( v1_G, syms_assumed, syms_assuming );
v2_G = subs( v2_G, syms_assumed, syms_assuming );
v_G_whole = subs( v_G_whole, syms_assumed, syms_assuming );
F2x = subs( F2x, syms_assumed, syms_assuming );
F2y = subs( F2y, syms_assumed, syms_assuming );
F1x = subs( F1x, syms_assumed, syms_assuming );
F1y = subs( F1y, syms_assumed, syms_assuming );

%% variable to be used for verification

% define Momentum
M = m1 * v1_G ...
    + m2 * v2_G ...
    ;
dM = formula( diff( M, t ) );

unclear_outer_force = [ F1x, F1y ];

% define Angular Momentum from p_G_whole
AM_G = ...
    formula( ...
        cross( ...
            [ p1_G, 0 ] - [p_G_whole, 0], ...
            m1 * ( [ v1_G, 0 ] - [ v_G_whole, 0 ] ) ...
        ) ...
        + cross( ...
            [ p2_G, 0 ] - [p_G_whole, 0], ...
            m2 * ( [ v2_G, 0 ] - [ v_G_whole, 0 ] ) ...
        ) ...
    );
AM_G = AM_G(3);

AM_each_rot = ...
    I1 * diff(th1_pre, t) ...
    + I2 * diff(th2_pre, t);
dAM = formula( diff( AM_each_rot, t ) + diff( AM_G, t ) );

unclear_outer_force_torque = ...
    cross( ...
        [p_G_whole, 0] - [p_G_whole, 0], ...
        [0, -m1 * g, 0] ...
    ) + ...
    cross( ...
        [0, 0, 0] - [p_G_whole, 0], ...
        [F1x, F1y, 0] ...
    );
unclear_outer_force_torque = formula( unclear_outer_force_torque );
unclear_outer_force_torque = unclear_outer_force_torque(3);

%% showing : functionalDerivative = dLdq - ddtdLddq for th1

dLdq = simplify( diff( L, th1_pre ) );
ddtdLddq = simplify( diff( L, diff( th1_pre, t ), t ) );

simplify( functionalDerivative( L, th1_pre ) - ( dLdq - ddtdLddq ) )

%% showing : functionalDerivative = dLdq - ddtdLddq for th2

dLdq = simplify( diff( L, th2_pre ) );
ddtdLddq = simplify( diff( L, diff( th2_pre, t ), t ) );

simplify( functionalDerivative( L, th2_pre ) - ( dLdq - ddtdLddq ) )

%% verify Newton eq
% -functionalDerivative is used. See https://www.mns.kyutech.ac.jp/~okamoto/education/mechanicsII/lagrangeeq_text.pdf

eq = [
    -functionalDerivative( L, th1_pre ) == tau1;
    -functionalDerivative( L, th2_pre ) == tau2;
    ];

eq = subs( eq, syms_replaced, syms_replacing );
dM = subs( dM, syms_replaced, syms_replacing );
dAM = subs( dAM, syms_replaced, syms_replacing );
unclear_outer_force = subs( unclear_outer_force, syms_replaced, syms_replacing );
unclear_outer_force_torque = subs( unclear_outer_force_torque, syms_replaced, syms_replacing );

variables = [ ddth1, ddth2 ];

[A, B] = equationsToMatrix( eq, variables );
% det(A) % verify det(A) ~= 0
X = inv(A) * B;

dM = subs( dM, variables, X' );
unclear_outer_force = subs( unclear_outer_force, variables, X' );

% this should be equal to [0, -g*(m1 + m2)]
simplify( dM - unclear_outer_force )

dAM = subs( dAM, variables, X' );
unclear_outer_force_torque = subs( unclear_outer_force_torque, variables, X' );

% this should be equal to tau1 + tau2
simplify( dAM - unclear_outer_force_torque )

%% make files
matlabFunction(X(1), X(2), 'file', 'find_ddth12.m', 'outputs', {'ddth1', 'ddth2'})







































