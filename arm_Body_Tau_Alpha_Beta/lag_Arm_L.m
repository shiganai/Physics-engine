
clear all
tic

parallel.defaultClusterProfile('local');
c = parcluster();

syms width_Body real

syms m_Hand real
syms length_Hand real
syms g real

syms l_X_Fixed l_Y_Fixed l_Z_Fixed real
syms l_F_X l_F_Y l_F_Z real
syms l_Tau_Alpha_Shoulder real
syms l_Tau_Beta_Shoulder real

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
l_Arm_G = [0, 0 + length_Hand/2, 0, 1];

%% rotate beta around z
l_Tauvec_Beta = symfun([0, 0, 1, 1], t);
l_Trans_Matrix_Beta = [cos(-l_Beta_Hand_Pre), -sin(-l_Beta_Hand_Pre), 0, 0; sin(-l_Beta_Hand_Pre), cos(-l_Beta_Hand_Pre), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]';

l_Arm_Bottom = l_Arm_Bottom * l_Trans_Matrix_Beta;
l_Arm_G = l_Arm_G * l_Trans_Matrix_Beta;

%% rotate alpha around x
l_Tauvec_Alpha = symfun([1, 0, 0, 1], t);
l_Trans_Matrix_Alpha = [1, 0, 0, 0; 0, cos(l_Alpha_Hand_Pre), -sin(l_Alpha_Hand_Pre), 0; 0, sin(l_Alpha_Hand_Pre), cos(l_Alpha_Hand_Pre), 0; 0, 0, 0, 1]';

l_Arm_Bottom = l_Arm_Bottom * l_Trans_Matrix_Alpha;
l_Arm_G = l_Arm_G * l_Trans_Matrix_Alpha;
l_Tauvec_Beta = l_Tauvec_Beta * l_Trans_Matrix_Alpha;

%% move origin
l_Trans_Matrix_Origin = [1, 0, 0, l_X_Fixed; 0, 1, 0, l_Y_Fixed; 0, 0, 1, l_Z_Fixed; 0, 0, 0, 1]';

l_Arm_Bottom = l_Arm_Bottom * l_Trans_Matrix_Origin;
l_Arm_G = l_Arm_G * l_Trans_Matrix_Origin;


%%
%{
l_Trans_Vec_Hand = ...
    [1, 0, 0, 0; 0, cos(l_Alpha_Hand_Pre), -sin(l_Alpha_Hand_Pre), 0; 0, sin(l_Alpha_Hand_Pre), cos(l_Alpha_Hand_Pre), 0; 0, 0, 0, 1]' ...
    * ...
    [cos(-l_Beta_Hand_Pre), -sin(-l_Beta_Hand_Pre), 0, 0; sin(-l_Beta_Hand_Pre), cos(-l_Beta_Hand_Pre), 0, 0; 0, 0, 1, 0; 0, 0, 0, 1]' ...
    * ...
    [1, 0, 0, l_X_Fixed; 0, 1, 0, l_Y_Fixed; 0, 0, 1, l_Z_Fixed; 0, 0, 0, 1]' ...
    * ...
    1;

l_Arm_Bottom = l_Arm_Bottom * l_Trans_Vec_Hand;
l_Arm_G = l_Arm_G * l_Trans_Vec_Hand;
%}

%%
l_Arm_Bottom = formula(l_Arm_Bottom);
l_Arm_G = formula(l_Arm_G);
l_Tauvec_Alpha = formula(l_Tauvec_Alpha);
l_Tauvec_Beta = formula(l_Tauvec_Beta);

l_Arm_Bottom = l_Arm_Bottom(1:3);
l_Arm_G = l_Arm_G(1:3);
l_Tauvec_Alpha = l_Tauvec_Alpha(1:3);
l_Tauvec_Beta = l_Tauvec_Beta(1:3);

l_V_G_Hand = diff(l_Arm_G, t);

T = ...
    1/2 * m_Hand * (l_V_G_Hand * l_V_G_Hand')...
    + ...
    1/2 * ([diff(l_Alpha_Hand_Pre, t), diff(l_Beta_Hand_Pre, t), 0] * (I_Hand * [diff(l_Alpha_Hand_Pre, t), diff(l_Beta_Hand_Pre, t), 0]'))...
    + ...
    0;

U = ...
    m_Hand * g * l_Arm_G(3)...
    +...
    0;

L = T - U;

l_Tau_Vec = (-l_Tau_Alpha_Shoulder) * l_Tauvec_Alpha + (-l_Tau_Beta_Shoulder) * l_Tauvec_Beta;

coeffs_Tau_L_Alpha = l_Tau_Vec * l_Tauvec_Alpha';
coeffs_Tau_L_Beta = l_Tau_Vec * l_Tauvec_Beta';

%%

equations = [
    -functionalDerivative(L, l_Alpha_Hand_Pre) == coeffs_Tau_L_Alpha - ([l_F_X, l_F_Y, l_F_Z] * diff(l_Arm_Bottom, l_Alpha_Hand_Pre)');
    -functionalDerivative(L, l_Beta_Hand_Pre) == coeffs_Tau_L_Beta - ([l_F_X, l_F_Y, l_F_Z] * diff(l_Arm_Bottom, l_Beta_Hand_Pre)');
    ];

equations = subs(equations, syms_Replaced, syms_Replacing);

%% Full forward dynamics
%{
variables = [ddl_Alpha_Hand, ddl_Beta_Hand];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

job = createJob(c);
createTask(job, @matlabFunction, 1,{X(1), X(2), ...
    'file', 'FFD_Dds_Arm_L.m', 'outputs', ...
    {'ddl_Alpha_Hand', 'ddl_Beta_Hand'}});
submit(job)
job.Tasks

ddl_Arm_Bottom = diff(l_Arm_Bottom, t, t);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, syms_Replaced, syms_Replacing);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, variables, X');

[coeffs_Ddl_Arm_Bottom(1, :), ~] = coeffs(ddl_Arm_Bottom(1), [l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Arm_Bottom(2, :), ~] = coeffs(ddl_Arm_Bottom(2), [l_F_X, l_F_Y, l_F_Z]);

[coeffs_Ddl_Arm_Bottom(3, :), ~] = coeffs(ddl_Arm_Bottom(3), [l_F_X, l_F_Y, l_F_Z]);

size(coeffs_Ddl_Arm_Bottom)

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Arm_Bottom(1,1), coeffs_Ddl_Arm_Bottom(1,2), coeffs_Ddl_Arm_Bottom(1,3), coeffs_Ddl_Arm_Bottom(1,4), ...
    coeffs_Ddl_Arm_Bottom(2,1), coeffs_Ddl_Arm_Bottom(2,2), coeffs_Ddl_Arm_Bottom(2,3), coeffs_Ddl_Arm_Bottom(2,4), ...
    coeffs_Ddl_Arm_Bottom(3,1), coeffs_Ddl_Arm_Bottom(3,2), coeffs_Ddl_Arm_Bottom(3,3), coeffs_Ddl_Arm_Bottom(3,4), ...
    'file', 'FFD_Coeffs_Ddl_Arm_Bottom.m', 'outputs', ...
    {...
    'A11','A12','A13','A14',...
    'A21','A22','A23','A24',...
    'A31','A32','A33','A34',...
    }});
submit(job)
job.Tasks
%}

%% Half forward dynamics
%{/
variables = [ddl_Alpha_Hand, l_Tau_Beta_Shoulder];

[A, B] = equationsToMatrix(equations, variables);
toc
tic
X = inv(A)*B;
toc

job = createJob(c);
createTask(job, @matlabFunction, 1,{X(1), X(2), ...
    'file', 'HFD_Dds_Arm_L.m', 'outputs', ...
    {'ddl_Alpha_Hand', 'l_Tau_Beta_Shoulder'}});
submit(job)
job.Tasks

ddl_Arm_Bottom = diff(l_Arm_Bottom, t, t);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, syms_Replaced, syms_Replacing);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, variables, X');

target_Variables = [l_F_X, l_F_Y, l_F_Z, 1];
coeffs_Ddl_Arm_Bottom = sym(zeros(3, size(target_Variables, 2)));

[coeffs_Tmp, terms_Tmp] = coeffs(ddl_Arm_Bottom(1), target_Variables(1:end-1));
if ~isequal(size(terms_Tmp), size(target_Variables))
    for ii = 1:size(target_Variables, 2)
        if any(terms_Tmp == target_Variables(ii))
            coeffs_Ddl_Arm_Bottom(1, ii) = coeffs_Tmp(terms_Tmp == target_Variables(ii));
        end
    end
else
    coeffs_Ddl_Arm_Bottom(1, :) = coeffs_Tmp;
end

[coeffs_Tmp, terms_Tmp] = coeffs(ddl_Arm_Bottom(2), target_Variables(1:end-1));
if ~isequal(size(terms_Tmp), size(target_Variables))
    for ii = 1:size(target_Variables, 2)
        if any(terms_Tmp == target_Variables(ii))
            coeffs_Ddl_Arm_Bottom(2, ii) = coeffs_Tmp(terms_Tmp == target_Variables(ii));
        end
    end
else
    coeffs_Ddl_Arm_Bottom(2, :) = coeffs_Tmp;
end

[coeffs_Tmp, terms_Tmp] = coeffs(ddl_Arm_Bottom(3), target_Variables(1:end-1));
if ~isequal(size(terms_Tmp), size(target_Variables))
    for ii = 1:size(target_Variables, 2)
        if any(terms_Tmp == target_Variables(ii))
            coeffs_Ddl_Arm_Bottom(3, ii) = coeffs_Tmp(terms_Tmp == target_Variables(ii));
        end
    end
else
    coeffs_Ddl_Arm_Bottom(3, :) = coeffs_Tmp;
end

size(coeffs_Ddl_Arm_Bottom)

coeffs_Ddl_Arm_Bottom_Force = coeffs_Ddl_Arm_Bottom(:, 1:3);
coeffs_Ddl_Arm_Bottom_Constant = coeffs_Ddl_Arm_Bottom(:, 4);

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Arm_Bottom_Force(1,1), coeffs_Ddl_Arm_Bottom_Force(1,2), coeffs_Ddl_Arm_Bottom_Force(1,3), ...
    coeffs_Ddl_Arm_Bottom_Force(2,1), coeffs_Ddl_Arm_Bottom_Force(2,2), coeffs_Ddl_Arm_Bottom_Force(2,3), ...
    coeffs_Ddl_Arm_Bottom_Force(3,1), coeffs_Ddl_Arm_Bottom_Force(3,2), coeffs_Ddl_Arm_Bottom_Force(3,3), ...
    'file', 'HFD_Coeffs_Ddl_Arm_Bottom_Force.m', 'outputs', ...
    {...
    'A11','A12','A13',...
    'A21','A22','A23',...
    'A31','A32','A33',...
    }});
submit(job)
job.Tasks

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_Ddl_Arm_Bottom_Constant(1,1), ...
    coeffs_Ddl_Arm_Bottom_Constant(2,1), ...
    coeffs_Ddl_Arm_Bottom_Constant(3,1), ...
    'file', 'HFD_Coeffs_Ddl_Arm_Bottom_Constant.m', 'outputs', ...
    {...
    'A11',...
    'A21',...
    'A31',...
    }});
submit(job)
job.Tasks

l_Tau_Beta_Shoulder = X(2);

[coeffs_L_Tau_Beta_Shoulder, ~] = coeffs(l_Tau_Beta_Shoulder(1), [l_F_X, l_F_Y, l_F_Z]);

coeffs_L_Tau_Beta_Shoulder_Force = coeffs_L_Tau_Beta_Shoulder(1:3);

job = createJob(c);
createTask(job, @matlabFunction, 1,{...
    coeffs_L_Tau_Beta_Shoulder_Force(1,1), coeffs_L_Tau_Beta_Shoulder_Force(1,2), coeffs_L_Tau_Beta_Shoulder_Force(1,3), ...
    'file', 'HFD_Coeffs_L_Tau_Beta_Shoulder_Force.m', 'outputs', ...
    {...
    'A11','A12','A13',...
    }});
submit(job)
job.Tasks
%}

%% Full Reverse dynamics
%{
ddl_Arm_Bottom = diff(l_Arm_Bottom, t, t);
ddl_Arm_Bottom = subs(ddl_Arm_Bottom, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{ddl_Arm_Bottom, ...
    'file', 'FRD_Ddl_Arm_Bottom.m', 'outputs', ...
    {'ddl_Arm_Bottom'}});
submit(job)
job.Tasks

dl_Arm_Bottom = diff(l_Arm_Bottom, t);
dl_Arm_Bottom = subs(dl_Arm_Bottom, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{dl_Arm_Bottom, ...
    'file', 'FRD_Dl_Arm_Bottom.m', 'outputs', ...
    {'dl_Arm_Bottom'}});
submit(job)
job.Tasks

l_Arm_Bottom = subs(l_Arm_Bottom, syms_Replaced, syms_Replacing);

job = createJob(c);
createTask(job, @matlabFunction, 1,{l_Arm_Bottom, ...
    'file', 'FRD_L_Arm_Bottom.m', 'outputs', ...
    {'l_Arm_Bottom'}});
submit(job)
job.Tasks
%}





































