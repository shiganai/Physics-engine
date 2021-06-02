
clear all
tic

syms width_Body real

syms m_Hand real
syms length_Hand real
syms g real

syms l_X_Fixed l_Y_Fixed l_Z_Fixed real
syms l_F_X l_F_Y l_F_Z real

syms l_Alpha_Hand_Pre(t)
syms l_Beta_Hand_Pre(t)

%%
syms l_Alpha_Hand dl_Alpha_Hand ddl_Alpha_Hand real
syms l_Beta_Hand dl_Beta_Hand ddl_Beta_Hand real

syms_Replaced = [
    l_Alpha_Hand_Pre diff(l_Alpha_Hand_Pre, t) diff(l_Alpha_Hand_Pre, t, t), ...
    l_Beta_Hand_Pre diff(l_Beta_Hand_Pre, t) diff(l_Beta_Hand_Pre, t, t), ...
    ];

syms_Replacing = [
    l_Alpha_Hand dl_Alpha_Hand ddl_Alpha_Hand ...
    l_Beta_Hand dl_Beta_Hand ddl_Beta_Hand ...
    ];

%%
I_Hand = 1/12 * m_Hand * [
    length_Hand^2 + 0^2, 0, 0;
    0, 0^2 + length_Hand^2, 0;
    0, 0, 0^2 + 0^2;
    ];

% I_Hand = 1/12 * m_Hand * [
%     length_Hand^2 + radius_Hand^2, 0, 0;
%     0, radius_Hand^2 + length_Hand^2, 0;
%     0, 0, radius_Hand^2 + radius_Hand^2;
%     ];

%%
l_Arm_Bottom = [0, 0 + length_Hand, 0, 1];
l_G_Hand = [0, 0 + length_Hand/2, 0, 1];

l_Trans_Vec_Hand = ...
    [1, 0, 0, 0; 0, cos(l_Alpha_Hand_Pre), -sin(l_Alpha_Hand_Pre), 0; 0, sin(l_Alpha_Hand_Pre), cos(l_Alpha_Hand_Pre), 0; 0, 0, 0, 1]' ...
    * ...
    [cos(-l_Beta_Hand_Pre), -sin(-l_Beta_Hand_Pre), 0, 0; sin(-l_Beta_Hand_Pre), cos(-l_Beta_Hand_Pre), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]' ...
    * ...
    [1, 0, 0, l_X_Fixed; 0, 1, 0, l_Y_Fixed; 0, 0, 1, l_Z_Fixed; 0, 0, 0, 1]' ...
    * ...
    1;

l_Arm_Bottom = l_Arm_Bottom * l_Trans_Vec_Hand;
l_G_Hand = l_G_Hand * l_Trans_Vec_Hand;

l_Arm_Bottom = formula(l_Arm_Bottom);
l_G_Hand = formula(l_G_Hand);

l_Arm_Bottom = l_Arm_Bottom(1:3);
l_G_Hand = l_G_Hand(1:3);

%%
l_V_G_Hand = diff(l_G_Hand, t);

T = ...
    1/2 * m_Hand * (l_V_G_Hand * l_V_G_Hand')...
    + ...
    1/2 * ([diff(l_Alpha_Hand_Pre, t), diff(l_Beta_Hand_Pre, t), 0] * (I_Hand * [diff(l_Alpha_Hand_Pre, t), diff(l_Beta_Hand_Pre, t), 0]'))...
    + ...
    0;

U = ...
    m_Hand * g * l_G_Hand(3)...
    +...
    0;

L = T - U;

%%

equations = [
    -functionalDerivative(L, l_Alpha_Hand_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Arm_Bottom, l_Alpha_Hand_Pre)');
    -functionalDerivative(L, l_Beta_Hand_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Arm_Bottom, l_Beta_Hand_Pre)');
    ];

%%

equations = subs(equations, syms_Replaced, syms_Replacing);

variables = [ddl_Alpha_Hand, ddl_Beta_Hand];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

parallel.defaultClusterProfile('local');
c = parcluster();

job = createJob(c);
createTask(job, @matlabFunction, 1,{X(1), X(2), ...
    'file', 'find_Dds_Arm_L.m', 'outputs', ...
    {'ddl_Alpha_Hand', 'ddl_Beta_Hand'}});
submit(job)
job.Tasks

%%

ddl_Arm_Bottom = diff(l_Arm_Bottom, t, t);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, syms_Replaced, syms_Replacing);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, variables, X');

[coeffs_Ddl_Arm_Bottom(1, :), ~] = coeffs(ddl_Arm_Bottom(1), [l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Arm_Bottom(2, :), ~] = coeffs(ddl_Arm_Bottom(2), [l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Arm_Bottom(3, :), ~] = coeffs(ddl_Arm_Bottom(3), [l_F_X, l_F_Y, l_F_Z]);

size(coeffs_Ddl_Arm_Bottom)

job = createJob(c);
createTask(job, @matlabFunction, 1,{coeffs_Ddl_Arm_Bottom, ...
    'file', 'find_Coeffs_Ddl_Arm_Bottom.m', 'outputs', ...
    {'coeffs_Ddl_Arm_Bottom'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddl_Arm_Bottom, ...
    'file', 'find_Ddl_Arm_Bottom.m', 'outputs', ...
    {'ddl_Arm_Bottom'}});
submit(job)
job.Tasks

%%
l_Arm_Bottom = subs(l_Arm_Bottom, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{l_Arm_Bottom, ...
    'file', 'find_L_Arm_Bottom.m', 'outputs', ...
    {'l_Arm_Bottom'}});
submit(job)
job.Tasks









































