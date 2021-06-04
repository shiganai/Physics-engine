
clear all
tic

syms width_Body real

syms m_Hand real
syms length_Hand real
syms g real

syms r_X_Fixed r_Y_Fixed r_Z_Fixed real
syms r_F_X r_F_Y r_F_Z real
syms r_Tau_Alpha_Shoulder real

syms r_Alpha_Hand_Pre(t)
syms r_Beta_Hand_Pre(t)

%%
syms r_Alpha_Hand dr_Alpha_Hand ddr_Alpha_Hand real
syms r_Beta_Hand dr_Beta_Hand ddr_Beta_Hand real

syms_Replaced = [
    r_Alpha_Hand_Pre diff(r_Alpha_Hand_Pre, t) diff(r_Alpha_Hand_Pre, t, t), ...
    r_Beta_Hand_Pre diff(r_Beta_Hand_Pre, t) diff(r_Beta_Hand_Pre, t, t), ...
    ];

syms_Replacing = [
    r_Alpha_Hand dr_Alpha_Hand ddr_Alpha_Hand ...
    r_Beta_Hand dr_Beta_Hand ddr_Beta_Hand ...
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
r_Arm_Bottom = [0, 0+ length_Hand, 0, 1];
r_G_Hand = [0, 0 + length_Hand/2, 0, 1];

r_Trans_Vec_Hand = ...
    [1, 0, 0, 0; 0, cos(r_Alpha_Hand_Pre), -sin(r_Alpha_Hand_Pre), 0; 0, sin(r_Alpha_Hand_Pre), cos(r_Alpha_Hand_Pre), 0; 0, 0, 0, 1]' ...
    * ...
    [cos(r_Beta_Hand_Pre), -sin(r_Beta_Hand_Pre), 0, 0; sin(r_Beta_Hand_Pre), cos(r_Beta_Hand_Pre), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]' ...
    * ...
    [1, 0, 0, r_X_Fixed; 0, 1, 0, r_Y_Fixed; 0, 0, 1, r_Z_Fixed; 0, 0, 0, 1]' ...
    * ...
    1;

r_Arm_Bottom = r_Arm_Bottom * r_Trans_Vec_Hand;
r_G_Hand = r_G_Hand * r_Trans_Vec_Hand;

r_Arm_Bottom = formula(r_Arm_Bottom);
r_G_Hand = formula(r_G_Hand);

r_Arm_Bottom = r_Arm_Bottom(1:3);
r_G_Hand = r_G_Hand(1:3);

%%
r_V_G_Hand = diff(r_G_Hand, t);

T = ...
    1/2 * m_Hand * (r_V_G_Hand * r_V_G_Hand')...
    + ...
    1/2 * ([diff(r_Alpha_Hand_Pre, t), diff(r_Beta_Hand_Pre, t), 0] * (I_Hand * [diff(r_Alpha_Hand_Pre, t), diff(r_Beta_Hand_Pre, t), 0]'))...
    + ...
    0;

U = ...
    m_Hand * g * r_G_Hand(3)...
    +...
    0;

L = T - U;

%%

equations = [
    -functionalDerivative(L, r_Alpha_Hand_Pre) == -r_Tau_Alpha_Shoulder - ([r_F_X, r_F_Y, r_F_Z] * diff(r_Arm_Bottom, r_Alpha_Hand_Pre)');
    -functionalDerivative(L, r_Beta_Hand_Pre) == 0 - ([r_F_X, r_F_Y, r_F_Z] * diff(r_Arm_Bottom, r_Beta_Hand_Pre)');
    ];

%%

equations = subs(equations, syms_Replaced, syms_Replacing);

variables = [ddr_Alpha_Hand, ddr_Beta_Hand];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

parallel.defaultClusterProfile('local');
c = parcluster();

job = createJob(c);
createTask(job, @matlabFunction, 1,{X(1), X(2), ...
    'file', 'find_Dds_Arm_R.m', 'outputs', ...
    {'ddr_Alpha_Hand', 'ddr_Beta_Hand'}});
submit(job)
job.Tasks

%%

ddr_Arm_Bottom = diff(r_Arm_Bottom, t, t);
ddr_Arm_Bottom = subs(ddr_Arm_Bottom, syms_Replaced, syms_Replacing);
ddr_Arm_Bottom = subs(ddr_Arm_Bottom, variables, X');

[coeffs_Ddr_Arm_Bottom(1, :), ~] = coeffs(ddr_Arm_Bottom(1), [r_F_X, r_F_Y, r_F_Z]);

[coeffs_Ddr_Arm_Bottom(2, :), ~] = coeffs(ddr_Arm_Bottom(2), [r_F_X, r_F_Y, r_F_Z]);

[coeffs_Ddr_Arm_Bottom(3, :), ~] = coeffs(ddr_Arm_Bottom(3), [r_F_X, r_F_Y, r_F_Z]);

size(coeffs_Ddr_Arm_Bottom)

job = createJob(c);
createTask(job, @matlabFunction, 1,{coeffs_Ddr_Arm_Bottom, ...
    'file', 'find_Coeffs_Ddr_Arm_Bottom.m', 'outputs', ...
    {'coeffs_Ddr_Arm_Bottom'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddr_Arm_Bottom, ...
    'file', 'find_Ddr_Arm_Bottom.m', 'outputs', ...
    {'ddr_Arm_Bottom'}});
submit(job)
job.Tasks

%%
dr_Arm_Bottom = diff(r_Arm_Bottom, t);
dr_Arm_Bottom = subs(dr_Arm_Bottom, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{dr_Arm_Bottom, ...
    'file', 'find_Dr_Arm_Bottom.m', 'outputs', ...
    {'dr_Arm_Bottom'}});
submit(job)
job.Tasks

%%
r_Arm_Bottom = subs(r_Arm_Bottom, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{r_Arm_Bottom, ...
    'file', 'find_R_Arm_Bottom.m', 'outputs', ...
    {'r_Arm_Bottom'}});
submit(job)
job.Tasks

































