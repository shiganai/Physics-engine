
clear all
tic

parallel.defaultClusterProfile('local');
c = parcluster();

syms t real

%% Start about r_Arm
syms r_Alpha_Hand real
syms r_Beta_Hand real
syms r_Tau_Alpha_Shoulder real
syms r_Tau_Beta_Shoulder real

%% rotate beta around z
r_Tauvec_Beta = symfun([0, 0, 1, 1], t);
r_Trans_Matrix_Beta = [cos(r_Beta_Hand), -sin(r_Beta_Hand), 0, 0; sin(r_Beta_Hand), cos(r_Beta_Hand), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]';

%% rotate alpha around x
r_Tauvec_Alpha = symfun([1, 0, 0, 1], t);
r_Trans_Matrix_Alpha = [1, 0, 0, 0; 0, cos(r_Alpha_Hand), -sin(r_Alpha_Hand), 0; 0, sin(r_Alpha_Hand), cos(r_Alpha_Hand), 0; 0, 0, 0, 1]';

r_Tauvec_Beta = r_Tauvec_Beta * r_Trans_Matrix_Alpha;

%% find r_Tau vector
r_Tauvec_Alpha = formula(r_Tauvec_Alpha);
r_Tauvec_Beta = formula(r_Tauvec_Beta);

r_Tauvec_Alpha = r_Tauvec_Alpha(1:3);
r_Tauvec_Beta = r_Tauvec_Beta(1:3);

r_Tau_Vec = r_Tau_Alpha_Shoulder * r_Tauvec_Alpha + r_Tau_Beta_Shoulder * r_Tauvec_Beta;

%% Start about l_Arm
syms l_Alpha_Hand real
syms l_Beta_Hand real
syms l_Tau_Alpha_Shoulder real
syms l_Tau_Beta_Shoulder real

%% rotate beta around z
l_Tauvec_Beta = symfun([0, 0, 1, 1], t);
l_Trans_Matrix_Beta = [cos(-l_Beta_Hand), -sin(-l_Beta_Hand), 0, 0; sin(-l_Beta_Hand), cos(-l_Beta_Hand), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]';

%% rotate alpha around x
l_Tauvec_Alpha = symfun([1, 0, 0, 1], t);
l_Trans_Matrix_Alpha = [1, 0, 0, 0; 0, cos(l_Alpha_Hand), -sin(l_Alpha_Hand), 0; 0, sin(l_Alpha_Hand), cos(l_Alpha_Hand), 0; 0, 0, 0, 1]';

l_Tauvec_Beta = l_Tauvec_Beta * l_Trans_Matrix_Alpha;

%% find l_Tau vector
l_Tauvec_Alpha = formula(l_Tauvec_Alpha);
l_Tauvec_Beta = formula(l_Tauvec_Beta);

l_Tauvec_Alpha = l_Tauvec_Alpha(1:3);
l_Tauvec_Beta = l_Tauvec_Beta(1:3);

l_Tau_Vec = l_Tau_Alpha_Shoulder * l_Tauvec_Alpha + l_Tau_Beta_Shoulder * l_Tauvec_Beta;

%% find external_Tau vector
external_Tau_Vec = r_Tau_Vec + l_Tau_Vec;

%% Start about body
syms width_Body height_Body depth_Body real
syms m_Body real
syms g real

syms r_F_X r_F_Y r_F_Z real
syms l_F_X l_F_Y l_F_Z real

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

%% rotate beta_Body around y
tauvec_Beta = symfun([0,1,0,1], t);
trans_Matrix_Beta = [cos(beta_Body_Pre), 0, sin(beta_Body_Pre), 0; 0, 1, 0, 0; -sin(beta_Body_Pre), 0, cos(beta_Body_Pre), 0; 0, 0, 0, 1]';

head = head * trans_Matrix_Beta;
r_Shoulder = r_Shoulder * trans_Matrix_Beta;
l_Shoulder = l_Shoulder * trans_Matrix_Beta;
r_Hip = r_Hip * trans_Matrix_Beta;
l_Hip = l_Hip * trans_Matrix_Beta;
g_Body = g_Body * trans_Matrix_Beta;

%% rotate gamma_Body around z
tauvec_Gamma = symfun([0,0,1,1], t);
trans_Matrix_Gamma = [cos(gamma_Body_Pre), -sin(gamma_Body_Pre), 0, 0; sin(gamma_Body_Pre), cos(gamma_Body_Pre), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]';

head = head * trans_Matrix_Gamma;
r_Shoulder = r_Shoulder * trans_Matrix_Gamma;
l_Shoulder = l_Shoulder * trans_Matrix_Gamma;
r_Hip = r_Hip * trans_Matrix_Gamma;
l_Hip = l_Hip * trans_Matrix_Gamma;
g_Body = g_Body * trans_Matrix_Gamma;

tauvec_Beta = tauvec_Beta * trans_Matrix_Gamma;

%% rotate beta_Body around x
tauvec_Alpha = symfun([1,0,0,1], t);
trans_Matrix_Alpha = [1, 0, 0, 0; 0, cos(alpha_Body_Pre), -sin(alpha_Body_Pre), 0; 0, sin(alpha_Body_Pre), cos(alpha_Body_Pre), 0; 0, 0, 0, 1]';

head = head * trans_Matrix_Alpha;
r_Shoulder = r_Shoulder * trans_Matrix_Alpha;
l_Shoulder = l_Shoulder * trans_Matrix_Alpha;
r_Hip = r_Hip * trans_Matrix_Alpha;
l_Hip = l_Hip * trans_Matrix_Alpha;
g_Body = g_Body * trans_Matrix_Alpha;

tauvec_Beta = tauvec_Beta * trans_Matrix_Alpha;
tauvec_Gamma = tauvec_Gamma * trans_Matrix_Alpha;

%% move origin
trans_Matrix_Origin = [1, 0, 0, x_Head_Pre; 0, 1, 0, y_Head_Pre; 0, 0, 1, z_Head_Pre; 0, 0, 0, 1]';

head = head * trans_Matrix_Origin;
r_Shoulder = r_Shoulder * trans_Matrix_Origin;
l_Shoulder = l_Shoulder * trans_Matrix_Origin;
r_Hip = r_Hip * trans_Matrix_Origin;
l_Hip = l_Hip * trans_Matrix_Origin;
g_Body = g_Body * trans_Matrix_Origin;

%%
%{
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
%}

%%
head = formula(head);
r_Shoulder = formula(r_Shoulder);
l_Shoulder = formula(l_Shoulder);
r_Hip = formula(r_Hip);
l_Hip = formula(l_Hip);
g_Body = formula(g_Body);
tauvec_Beta = formula(tauvec_Beta);
tauvec_Gamma = formula(tauvec_Gamma);
tauvec_Alpha = formula(tauvec_Alpha);

r_Shoulder = r_Shoulder(1:3);
l_Shoulder = l_Shoulder(1:3);
r_Hip = r_Hip(1:3);
l_Hip = l_Hip(1:3);
g_Body = g_Body(1:3);
tauvec_Beta = tauvec_Beta(1:3);
tauvec_Gamma = tauvec_Gamma(1:3);
tauvec_Alpha = tauvec_Alpha(1:3);

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

coeffs_Tau_Alpha_Body = external_Tau_Vec * tauvec_Alpha';
coeffs_Tau_Beta_Body = external_Tau_Vec * tauvec_Beta';
coeffs_Tau_Gamma_Body = external_Tau_Vec * tauvec_Gamma';

%%
equations = [
    -functionalDerivative(L, alpha_Body_Pre) == coeffs_Tau_Alpha_Body + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, alpha_Body_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, alpha_Body_Pre)');
    -functionalDerivative(L, beta_Body_Pre) == coeffs_Tau_Beta_Body + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, beta_Body_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, beta_Body_Pre)');
    -functionalDerivative(L, gamma_Body_Pre) == coeffs_Tau_Gamma_Body + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, gamma_Body_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, gamma_Body_Pre)');
    -functionalDerivative(L, x_Head_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, x_Head_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, x_Head_Pre)');
    -functionalDerivative(L, y_Head_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, y_Head_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, y_Head_Pre)');
    -functionalDerivative(L, z_Head_Pre) == 0 + ([l_F_X, l_F_Y, l_F_Z] * diff(l_Shoulder, z_Head_Pre)') + ([r_F_X, r_F_Y, r_F_Z] * diff(r_Shoulder, z_Head_Pre)');
    ];

equations = subs(equations, syms_Replaced, syms_Replacing);

%% Full forward dynamics
%{
variables = [ddalpha_Body, ddbeta_Body, ddgamma_Body, ddx_Head, ddy_Head, ddz_Head];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

job = createJob(c);
createTask(job, @matlabFunction, 1,{X(1), X(2), X(3), X(4), X(5), X(6), ...
    'file', 'FFD_Dds_Body.m', 'outputs', ...
    {'ddalpha_Body', 'ddbeta_Body', 'ddgamma_Body', 'ddx_Head', 'ddy_Head', 'ddz_Head'}});
submit(job)
job.Tasks

ddr_Shoulder = diff(r_Shoulder, t, t);
ddr_Shoulder = subs(ddr_Shoulder, syms_Replaced, syms_Replacing);
ddr_Shoulder = subs(ddr_Shoulder, variables, X');

[coeffs_Ddr_Shoulder_FFD(1, :), ~] = coeffs(ddr_Shoulder(1), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddr_Shoulder_FFD(2, :), ~] = coeffs(ddr_Shoulder(2), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddr_Shoulder_FFD(3, :), ~] = coeffs(ddr_Shoulder(3), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

size(coeffs_Ddr_Shoulder_FFD)

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddr_Shoulder_FFD(1,1), coeffs_Ddr_Shoulder_FFD(1,2), coeffs_Ddr_Shoulder_FFD(1,3), coeffs_Ddr_Shoulder_FFD(1,4), coeffs_Ddr_Shoulder_FFD(1,5), coeffs_Ddr_Shoulder_FFD(1,6), coeffs_Ddr_Shoulder_FFD(1,7), ...
    coeffs_Ddr_Shoulder_FFD(2,1), coeffs_Ddr_Shoulder_FFD(2,2), coeffs_Ddr_Shoulder_FFD(2,3), coeffs_Ddr_Shoulder_FFD(2,4), coeffs_Ddr_Shoulder_FFD(2,5), coeffs_Ddr_Shoulder_FFD(2,6), coeffs_Ddr_Shoulder_FFD(2,7), ...
    coeffs_Ddr_Shoulder_FFD(3,1), coeffs_Ddr_Shoulder_FFD(3,2), coeffs_Ddr_Shoulder_FFD(3,3), coeffs_Ddr_Shoulder_FFD(3,4), coeffs_Ddr_Shoulder_FFD(3,5), coeffs_Ddr_Shoulder_FFD(3,6), coeffs_Ddr_Shoulder_FFD(3,7), ...
    'file', 'FFD_Coeffs_Ddr_Shoulder.m', 'outputs', ...
    {...
    'A11','A12','A13','A14','A15','A16','A17',...
    'A21','A22','A23','A24','A25','A26','A27',...
    'A31','A32','A33','A34','A35','A36','A37',...
    }});
submit(job)
job.Tasks

ddl_Shoulder = diff(l_Shoulder, t, t);
ddl_Shoulder = subs(ddl_Shoulder, syms_Replaced, syms_Replacing);
ddl_Shoulder = subs(ddl_Shoulder, variables, X');

[coeffs_Ddl_Shoulder_FFD(1, :), ~] = coeffs(ddl_Shoulder(1), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Shoulder_FFD(2, :), ~] = coeffs(ddl_Shoulder(2), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Shoulder_FFD(3, :), ~] = coeffs(ddl_Shoulder(3), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z]);

size(coeffs_Ddl_Shoulder_FFD)

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Shoulder_FFD(1,1), coeffs_Ddl_Shoulder_FFD(1,2), coeffs_Ddl_Shoulder_FFD(1,3), coeffs_Ddl_Shoulder_FFD(1,4), coeffs_Ddl_Shoulder_FFD(1,5), coeffs_Ddl_Shoulder_FFD(1,6), coeffs_Ddl_Shoulder_FFD(1,7), ...
    coeffs_Ddl_Shoulder_FFD(2,1), coeffs_Ddl_Shoulder_FFD(2,2), coeffs_Ddl_Shoulder_FFD(2,3), coeffs_Ddl_Shoulder_FFD(2,4), coeffs_Ddl_Shoulder_FFD(2,5), coeffs_Ddl_Shoulder_FFD(2,6), coeffs_Ddl_Shoulder_FFD(2,7), ...
    coeffs_Ddl_Shoulder_FFD(3,1), coeffs_Ddl_Shoulder_FFD(3,2), coeffs_Ddl_Shoulder_FFD(3,3), coeffs_Ddl_Shoulder_FFD(3,4), coeffs_Ddl_Shoulder_FFD(3,5), coeffs_Ddl_Shoulder_FFD(3,6), coeffs_Ddl_Shoulder_FFD(3,7), ...
    'file', 'FFD_Coeffs_Ddl_Shoulder.m', 'outputs', ...
    {...
    'A11','A12','A13','A14','A15','A16','A17',...
    'A21','A22','A23','A24','A25','A26','A27',...
    'A31','A32','A33','A34','A35','A36','A37',...
    }});
submit(job)
job.Tasks
%}

%% Half forward dynamics
%{/
variables = [ddalpha_Body, ddbeta_Body, ddgamma_Body, ddx_Head, ddy_Head, ddz_Head];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

% no need
% job = createJob(c);
% createTask(job, @matlabFunction, 1,{X(1), X(2), X(3), X(4), X(5), X(6), ...
%     'file', 'find_Dds_Body_Tau_Beta_HFD.m', 'outputs', ...
%     {'ddalpha_Body', 'ddx_Head', 'ddy_Head', 'ddz_Head', 'r_Tau_Beta_Shoulder', 'l_Tau_Beta_Shoulder'}});
% submit(job)
% job.Tasks

ddr_Shoulder = diff(r_Shoulder, t, t);
ddr_Shoulder = subs(ddr_Shoulder, syms_Replaced, syms_Replacing);
ddr_Shoulder = subs(ddr_Shoulder, variables, X');

[coeffs_Ddr_Shoulder_HFD(1, :), ~] = coeffs(ddr_Shoulder(1), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z, r_Tau_Beta_Shoulder, l_Tau_Beta_Shoulder]);
[coeffs_Ddr_Shoulder_HFD(2, :), ~] = coeffs(ddr_Shoulder(2), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z, r_Tau_Beta_Shoulder, l_Tau_Beta_Shoulder]);
[coeffs_Ddr_Shoulder_HFD(3, :), ~] = coeffs(ddr_Shoulder(3), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z, r_Tau_Beta_Shoulder, l_Tau_Beta_Shoulder]);

size(coeffs_Ddr_Shoulder_HFD)

coeffs_Ddr_Shoulder_Force = coeffs_Ddr_Shoulder_HFD(:, 1:6);
coeffs_Ddr_Shoulder_Tau_Beta = coeffs_Ddr_Shoulder_HFD(:, 7:8);
coeffs_Ddr_Shoulder_Constant = coeffs_Ddr_Shoulder_HFD(:, 9);

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddr_Shoulder_Force(1,1), coeffs_Ddr_Shoulder_Force(1,2), coeffs_Ddr_Shoulder_Force(1,3), coeffs_Ddr_Shoulder_Force(1,4), coeffs_Ddr_Shoulder_Force(1,5), coeffs_Ddr_Shoulder_Force(1,6), ...
    coeffs_Ddr_Shoulder_Force(2,1), coeffs_Ddr_Shoulder_Force(2,2), coeffs_Ddr_Shoulder_Force(2,3), coeffs_Ddr_Shoulder_Force(2,4), coeffs_Ddr_Shoulder_Force(2,5), coeffs_Ddr_Shoulder_Force(2,6), ...
    coeffs_Ddr_Shoulder_Force(3,1), coeffs_Ddr_Shoulder_Force(3,2), coeffs_Ddr_Shoulder_Force(3,3), coeffs_Ddr_Shoulder_Force(3,4), coeffs_Ddr_Shoulder_Force(3,5), coeffs_Ddr_Shoulder_Force(3,6), ...
    'file', 'HFD_Coeffs_Ddr_Shoulder_Force.m', 'outputs', ...
    {...
    'A11','A12','A13','A14','A15','A16',...
    'A21','A22','A23','A24','A25','A26',...
    'A31','A32','A33','A34','A35','A36',...
    }});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddr_Shoulder_Tau_Beta(1,1), coeffs_Ddr_Shoulder_Tau_Beta(1,2), ...
    coeffs_Ddr_Shoulder_Tau_Beta(2,1), coeffs_Ddr_Shoulder_Tau_Beta(2,2), ...
    coeffs_Ddr_Shoulder_Tau_Beta(3,1), coeffs_Ddr_Shoulder_Tau_Beta(3,2), ...
    'file', 'HFD_Coeffs_Ddr_Shoulder_Tau_Beta.m', 'outputs', ...
    {...
    'A11','A12',...
    'A21','A22',...
    'A31','A32',...
    }});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddr_Shoulder_Constant(1,1), ...
    coeffs_Ddr_Shoulder_Constant(2,1), ...
    coeffs_Ddr_Shoulder_Constant(3,1), ...
    'file', 'HFD_Coeffs_Ddr_Shoulder_Constant.m', 'outputs', ...
    {...
    'A11',...
    'A21',...
    'A31',...
    }});
submit(job)
job.Tasks

ddl_Shoulder = diff(l_Shoulder, t, t);
ddl_Shoulder = subs(ddl_Shoulder, syms_Replaced, syms_Replacing);
ddl_Shoulder = subs(ddl_Shoulder, variables, X');

[coeffs_Ddl_Shoulder_HFD(1, :), ~] = coeffs(ddl_Shoulder(1), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z, r_Tau_Beta_Shoulder, l_Tau_Beta_Shoulder]);
[coeffs_Ddl_Shoulder_HFD(2, :), ~] = coeffs(ddl_Shoulder(2), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z, r_Tau_Beta_Shoulder, l_Tau_Beta_Shoulder]);
[coeffs_Ddl_Shoulder_HFD(3, :), ~] = coeffs(ddl_Shoulder(3), [r_F_X, r_F_Y, r_F_Z, l_F_X, l_F_Y, l_F_Z, r_Tau_Beta_Shoulder, l_Tau_Beta_Shoulder]);

size(coeffs_Ddl_Shoulder_HFD)

coeffs_Ddl_Shoulder_Force = coeffs_Ddl_Shoulder_HFD(:, 1:6);
coeffs_Ddl_Shoulder_Tau_Beta = coeffs_Ddl_Shoulder_HFD(:, 7:8);
coeffs_Ddl_Shoulder_Constant = coeffs_Ddl_Shoulder_HFD(:, 9);

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Shoulder_Force(1,1), coeffs_Ddl_Shoulder_Force(1,2), coeffs_Ddl_Shoulder_Force(1,3), coeffs_Ddl_Shoulder_Force(1,4), coeffs_Ddl_Shoulder_Force(1,5), coeffs_Ddl_Shoulder_Force(1,6), ...
    coeffs_Ddl_Shoulder_Force(2,1), coeffs_Ddl_Shoulder_Force(2,2), coeffs_Ddl_Shoulder_Force(2,3), coeffs_Ddl_Shoulder_Force(2,4), coeffs_Ddl_Shoulder_Force(2,5), coeffs_Ddl_Shoulder_Force(2,6), ...
    coeffs_Ddl_Shoulder_Force(3,1), coeffs_Ddl_Shoulder_Force(3,2), coeffs_Ddl_Shoulder_Force(3,3), coeffs_Ddl_Shoulder_Force(3,4), coeffs_Ddl_Shoulder_Force(3,5), coeffs_Ddl_Shoulder_Force(3,6), ...
    'file', 'HFD_Coeffs_Ddl_Shoulder_Force.m', 'outputs', ...
    {...
    'A11','A12','A13','A14','A15','A16',...
    'A21','A22','A23','A24','A25','A26',...
    'A31','A32','A33','A34','A35','A36',...
    }});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Shoulder_Tau_Beta(1,1), coeffs_Ddl_Shoulder_Tau_Beta(1,2), ...
    coeffs_Ddl_Shoulder_Tau_Beta(2,1), coeffs_Ddl_Shoulder_Tau_Beta(2,2), ...
    coeffs_Ddl_Shoulder_Tau_Beta(3,1), coeffs_Ddl_Shoulder_Tau_Beta(3,2), ...
    'file', 'HFD_Coeffs_Ddl_Shoulder_Tau_Beta.m', 'outputs', ...
    {...
    'A11','A12',...
    'A21','A22',...
    'A31','A32',...
    }});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Shoulder_Constant(1,1), ...
    coeffs_Ddl_Shoulder_Constant(2,1), ...
    coeffs_Ddl_Shoulder_Constant(3,1), ...
    'file', 'HFD_Coeffs_Ddl_Shoulder_Constant.m', 'outputs', ...
    {...
    'A11',...
    'A21',...
    'A31',...
    }});
submit(job)
job.Tasks
%}
%% Full Reverse dynamics
%{
ddr_Shoulder = diff(r_Shoulder, t, t);
ddr_Shoulder = subs(ddr_Shoulder, syms_Replaced, syms_Replacing);
ddl_Shoulder = diff(l_Shoulder, t, t);
ddl_Shoulder = subs(ddl_Shoulder, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddr_Shoulder, ...
    'file', 'FRD_Ddr_Shoulder.m', 'outputs', ...
    {'ddr_Shoulder'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddl_Shoulder, ...
    'file', 'FRD_Ddl_Shoulder.m', 'outputs', ...
    {'ddl_Shoulder'}});
submit(job)
job.Tasks

r_Shoulder = subs(r_Shoulder, syms_Replaced, syms_Replacing);
l_Shoulder = subs(l_Shoulder, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{l_Shoulder, ...
    'file', 'FRD_L_Shoulder.m', 'outputs', ...
    {'l_Shoulder'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{r_Shoulder, ...
    'file', 'FRD_R_Shoulder.m', 'outputs', ...
    {'r_Shoulder'}});
submit(job)
job.Tasks

r_Hip = subs(r_Hip, syms_Replaced, syms_Replacing);
l_Hip = subs(l_Hip, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{l_Hip, ...
    'file', 'FRD_L_Hip.m', 'outputs', ...
    {'l_Hip'}});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{r_Hip, ...
    'file', 'FRD_R_Hip.m', 'outputs', ...
    {'r_Hip'}});
submit(job)
job.Tasks
%}



































