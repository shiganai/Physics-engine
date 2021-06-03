
clear all

syms m_Hand real
syms length_Hand real
syms g real

syms width_Body height_Body depth_Body real
syms m_Body real


syms l_X_Fixed l_Y_Fixed l_Z_Fixed real
syms l_F_X l_F_Y l_F_Z real

syms r_X_Fixed r_Y_Fixed r_Z_Fixed real
syms r_F_X r_F_Y r_F_Z real

syms l_Alpha_Hand_Pre(t)
syms l_Beta_Hand_Pre(t)

syms r_Alpha_Hand_Pre(t)
syms r_Beta_Hand_Pre(t)


syms offset_Body_Pre(t)
syms alpha_Body_Pre(t)
syms beta_Body_Pre(t)
syms gamma_Body_Pre(t)
syms x_Head_Pre(t) y_Head_Pre(t) z_Head_Pre(t)

%%
syms l_Alpha_Hand dl_Alpha_Hand ddl_Alpha_Hand real
syms l_Beta_Hand dl_Beta_Hand ddl_Beta_Hand real

syms r_Alpha_Hand dr_Alpha_Hand ddr_Alpha_Hand real
syms r_Beta_Hand dr_Beta_Hand ddr_Beta_Hand real

syms alpha_Body dalpha_Body ddalpha_Body real
syms beta_Body dbeta_Body ddbeta_Body real
syms gamma_Body dgamma_Body ddgamma_Body real
syms x_Head dx_Head ddx_Head real
syms y_Head dy_Head ddy_Head real
syms z_Head dz_Head ddz_Head real

syms_Replaced = [
    l_Alpha_Hand_Pre diff(l_Alpha_Hand_Pre, t) diff(l_Alpha_Hand_Pre, t, t), ...
    l_Beta_Hand_Pre diff(l_Beta_Hand_Pre, t) diff(l_Beta_Hand_Pre, t, t), ...
    r_Alpha_Hand_Pre diff(r_Alpha_Hand_Pre, t) diff(r_Alpha_Hand_Pre, t, t), ...
    r_Beta_Hand_Pre diff(r_Beta_Hand_Pre, t) diff(r_Beta_Hand_Pre, t, t), ...
    alpha_Body_Pre diff(alpha_Body_Pre, t) diff(alpha_Body_Pre, t, t), ...
    beta_Body_Pre diff(beta_Body_Pre, t) diff(beta_Body_Pre, t, t), ...
    gamma_Body_Pre diff(gamma_Body_Pre, t) diff(gamma_Body_Pre, t, t), ...
    x_Head_Pre diff(x_Head_Pre, t) diff(x_Head_Pre, t, t), ...
    y_Head_Pre diff(y_Head_Pre, t) diff(y_Head_Pre, t, t), ...
    z_Head_Pre diff(z_Head_Pre, t) diff(z_Head_Pre, t, t), ...
    ];

syms_Replacing = [
    l_Alpha_Hand dl_Alpha_Hand ddl_Alpha_Hand ...
    l_Beta_Hand dl_Beta_Hand ddl_Beta_Hand ...
    r_Alpha_Hand dr_Alpha_Hand ddr_Alpha_Hand ...
    r_Beta_Hand dr_Beta_Hand ddr_Beta_Hand ...
    alpha_Body dalpha_Body ddalpha_Body ...
    beta_Body dbeta_Body ddbeta_Body ...
    gamma_Body dgamma_Body ddgamma_Body ...
    x_Head dx_Head ddx_Head ...
    y_Head dy_Head ddy_Head ...
    z_Head dz_Head ddz_Head ...
    ];

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
restraint = [
    r_Arm_Bottom' == r_Shoulder';
    l_Arm_Bottom' == l_Shoulder';
    ];

restraint_Velocity = diff(restraint, t);
restraint_Velocity = subs(restraint_Velocity, syms_Replaced, syms_Replacing);

variables_Velocity = [dbeta_Body, dx_Head, dy_Head, dz_Head];
[A,B] = equationsToMatrix(restraint_Velocity(2:5), variables_Velocity);
% variables_Velocity = [dbeta_Body, dgamma_Body, dx_Head, dy_Head, dz_Head];
% [A,B] = equationsToMatrix(restraint_Velocity(1:5), variables_Velocity);
% variables_Velocity = [dalpha_Body, dbeta_Body, dgamma_Body, dx_Head, dy_Head, dz_Head];
% [A,B] = equationsToMatrix(restraint_Velocity, variables_Velocity);

sol_Velocity = inv(A) * B;

parallel.defaultClusterProfile('local');
c = parcluster();

job = createJob(c);
createTask(job, @matlabFunction, 1,{sol_Velocity(1), sol_Velocity(2), sol_Velocity(3), sol_Velocity(4), ...
    'file', 'find_Restrained_Velocity.m', 'outputs', ...
    {'dbeta_Body', 'dx_Head', 'dy_Head', 'dz_Head'}});
submit(job)
job.Tasks

% job = createJob(c);
% createTask(job, @matlabFunction, 1,{sol_Velocity(1), sol_Velocity(2), sol_Velocity(3), sol_Velocity(4), sol_Velocity(5), ...
%     'file', 'find_Restrained_Velocity.m', 'outputs', ...
%     {'dbeta_Body', 'dgamma_Body', 'dx_Head', 'dy_Head', 'dz_Head'}});
% submit(job)
% job.Tasks

% job = createJob(c);
% createTask(job, @matlabFunction, 1,{sol_Velocity(1), sol_Velocity(2), sol_Velocity(3), sol_Velocity(4), sol_Velocity(5), sol_Velocity(6), ...
%     'file', 'find_Restrained_Velocity.m', 'outputs', ...
%     {'dalpha_Body', 'dbeta_Body', 'dgamma_Body', 'dx_Head', 'dy_Head', 'dz_Head'}});
% submit(job)
% job.Tasks















































