
clear all
tic

syms width_Body height_Body depth_Body real
syms m_Body real
syms g real

syms r_F_X r_F_Y r_F_Z real
syms l_F_X l_F_Y l_F_Z real

syms offset_Body_Pre(t)
syms alpha_Body_Pre(t)
syms beta_Body_Pre(t)
syms gamma_Body_Pre(t)
syms x_Head_Pre(t) y_Head_Pre(t) z_Head_Pre(t)

%%

syms alpha_Body dalpha_Body ddalpha_Body real
syms beta_Body dbeta_Body ddbeta_Body real
syms gamma_Body dgamma_Body ddgamma_Body real
syms x_Head dx_Head ddx_Head real
syms y_Head dy_Head ddy_Head real
syms z_Head dz_Head ddz_Head real

syms_Replaced = [
    alpha_Body_Pre diff(alpha_Body_Pre, t) diff(alpha_Body_Pre, t, t), ...
    beta_Body_Pre diff(beta_Body_Pre, t) diff(beta_Body_Pre, t, t), ...
    gamma_Body_Pre diff(gamma_Body_Pre, t) diff(gamma_Body_Pre, t, t), ...
    x_Head_Pre diff(x_Head_Pre, t) diff(x_Head_Pre, t, t), ...
    y_Head_Pre diff(y_Head_Pre, t) diff(y_Head_Pre, t, t), ...
    z_Head_Pre diff(z_Head_Pre, t) diff(z_Head_Pre, t, t), ...
    ];

syms_Replacing = [
    alpha_Body dalpha_Body ddalpha_Body ...
    beta_Body dbeta_Body ddbeta_Body ...
    gamma_Body dgamma_Body ddgamma_Body ...
    x_Head dx_Head ddx_Head ...
    y_Head dy_Head ddy_Head ...
    z_Head dz_Head ddz_Head ...
    ];

%%
I_Body = 1/12 * m_Body * [
    height_Body^2 + depth_Body^2, 0, 0;
    0, width_Body^2 + depth_Body^2, 0;
    0, 0, width_Body^2 + height_Body^2;
    ];

%%
head = [0, 0, 0, 1];
r_Shoulder = [width_Body/2, 0, 0, 1];
l_Shoulder = [-width_Body/2, 0, 0, 1];
r_Hip = [width_Body/2, height_Body, 0, 1];
l_Hip = [-width_Body/2, height_Body, 0, 1];
g_Body = (r_Shoulder + l_Shoulder + r_Hip + l_Hip)/4;

trans_Vec_Body = ...
    [cos(beta_Body_Pre), 0, sin(beta_Body_Pre), 0; 0, 1, 0, 0; -sin(beta_Body_Pre), 0, cos(beta_Body_Pre), 0; 0, 0, 0, 1]' ...
    * ...
    [cos(gamma_Body_Pre), -sin(gamma_Body_Pre), 0, 0; sin(gamma_Body_Pre), cos(gamma_Body_Pre), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]' ...
    * ...
    [1, 0, 0, 0; 0, cos(alpha_Body_Pre), -sin(alpha_Body_Pre), 0; 0, sin(alpha_Body_Pre), cos(alpha_Body_Pre), 0; 0, 0, 0, 1]' ...
    * ...
    [1, 0, 0, x_Head_Pre; 0, 1, 0, y_Head_Pre; 0, 0, 1, z_Head_Pre; 0, 0, 0, 1]' ...
    * ...
    1;

head = formula(head * trans_Vec_Body);
r_Shoulder = formula(r_Shoulder * trans_Vec_Body);
l_Shoulder = formula(l_Shoulder * trans_Vec_Body);
r_Hip = formula(r_Hip * trans_Vec_Body);
l_Hip = formula(l_Hip * trans_Vec_Body);
g_Body = formula(g_Body * trans_Vec_Body);

r_Shoulder = r_Shoulder(1:3);
l_Shoulder = l_Shoulder(1:3);
r_Hip = r_Hip(1:3);
l_Hip = l_Hip(1:3);
g_Body = g_Body(1:3);

%%
v_G_Body = diff(g_Body, t);

T = ...
    1/2 * m_Body * (v_G_Body * v_G_Body')...
    + ...
    1/2 * ([diff(alpha_Body_Pre, t), diff(beta_Body_Pre, t), diff(gamma_Body_Pre, t)] * (I_Body * [diff(alpha_Body_Pre, t), diff(beta_Body_Pre, t), diff(gamma_Body_Pre, t)]'))...
    + ...
    0;

U = ...
    m_Body * g * g_Body(3)...
    +...
    0;

L = T - U;

%%
equations = [
    -functionalDerivative(L, alpha_Body_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, alpha_Body_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, alpha_Body_Pre)');
    -functionalDerivative(L, beta_Body_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, beta_Body_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, beta_Body_Pre)');
    -functionalDerivative(L, gamma_Body_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, gamma_Body_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, gamma_Body_Pre)');
    -functionalDerivative(L, x_Head_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, x_Head_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, x_Head_Pre)');
    -functionalDerivative(L, y_Head_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, y_Head_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, y_Head_Pre)');
    -functionalDerivative(L, z_Head_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, z_Head_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, z_Head_Pre)');
    ];

%%

equations = subs(equations, syms_Replaced, syms_Replacing);

variables = [ddalpha_Body, ddbeta_Body, ddgamma_Body, ddx_Head, ddy_Head, ddz_Head];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

parallel.defaultClusterProfile('local');
c = parcluster();

job = createJob(c);
createTask(job, @matlabFunction, 1,{X(1), X(2), X(3), X(4), X(5), X(6), ...
    'file', 'find_Dds_Body.m', 'outputs', ...
    {'ddalpha_Body', 'ddbeta_Body', 'ddgamma_Body', 'ddx_Head', 'ddy_Head', 'ddz_Head'}});
submit(job)
job.Tasks

%%

ddoffset_Body = X(1);
ddalpha_Body = X(2);
ddbeta_Body = X(3);
ddgamma_Body = X(4);

%%

ddr_Shoulder = diff(r_Shoulder, t, t);
ddr_Shoulder = subs(ddr_Shoulder, syms_Replaced, syms_Replacing);
ddr_Shoulder = subs(ddr_Shoulder, variables, X');

[coeffs_Ddr_Shoulde(1, :), ~] = coeffs(ddr_Shoulder(1), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddr_Shoulde(2, :), ~] = coeffs(ddr_Shoulder(2), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddr_Shoulde(3, :), ~] = coeffs(ddr_Shoulder(3), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

size(coeffs_Ddr_Shoulde)

job = createJob(c);
createTask(job, @matlabFunction, 1,{coeffs_Ddr_Shoulde, ...
    'file', 'find_Coeffs_Ddr_Shoulder.m', 'outputs', ...
    {'coeffs_Ddr_Shoulder'}});
submit(job)
job.Tasks

ddl_Shoulder = diff(l_Shoulder, t, t);
ddl_Shoulder = subs(ddl_Shoulder, syms_Replaced, syms_Replacing);
ddl_Shoulder = subs(ddl_Shoulder, variables, X');

[coeffs_Ddl_Shoulde(1, :), ~] = coeffs(ddl_Shoulder(1), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Shoulde(2, :), ~] = coeffs(ddl_Shoulder(2), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Shoulde(3, :), ~] = coeffs(ddl_Shoulder(3), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

size(coeffs_Ddl_Shoulde)

job = createJob(c);
createTask(job, @matlabFunction, 1,{coeffs_Ddl_Shoulde, ...
    'file', 'find_Coeffs_Ddl_Shoulder.m', 'outputs', ...
    {'coeffs_Ddl_Shoulder'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddr_Shoulder, ...
    'file', 'find_Ddr_Shoulder.m', 'outputs', ...
    {'ddr_Shoulder'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddl_Shoulder, ...
    'file', 'find_Ddl_Shoulder.m', 'outputs', ...
    {'ddl_Shoulder'}});
submit(job)
job.Tasks

%%
r_Shoulder = subs(r_Shoulder, syms_Replaced, syms_Replacing);
l_Shoulder = subs(l_Shoulder, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{l_Shoulder, ...
    'file', 'find_L_Shoulder.m', 'outputs', ...
    {'l_Shoulder'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{r_Shoulder, ...
    'file', 'find_R_Shoulder.m', 'outputs', ...
    {'r_Shoulder'}});
submit(job)
job.Tasks

%%
r_Hip = subs(r_Hip, syms_Replaced, syms_Replacing);
l_Hip = subs(l_Hip, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{l_Hip, ...
    'file', 'find_L_Hip.m', 'outputs', ...
    {'l_Hip'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{r_Hip, ...
    'file', 'find_R_Hip.m', 'outputs', ...
    {'r_Hip'}});
submit(job)
job.Tasks



































