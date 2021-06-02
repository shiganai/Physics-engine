function dotq = ddt(~,q,r_P_Fixed,l_P_Fixed, g, length_Hand, m_Hand, m_Body, width_Body, height_Body, depth_Body)

r_X_Fixed = r_P_Fixed(1);
r_Y_Fixed = r_P_Fixed(2);
r_Z_Fixed = r_P_Fixed(3);

l_X_Fixed = l_P_Fixed(1);
l_Y_Fixed = l_P_Fixed(2);
l_Z_Fixed = l_P_Fixed(3);

r_Alpha_Hand = q(1);
dr_Alpha_Hand = q(2);

r_Beta_Hand = q(3);
dr_Beta_Hand = q(4);

l_Alpha_Hand = q(5);
dl_Alpha_Hand = q(6);

l_Beta_Hand = q(7);
dl_Beta_Hand = q(8);

offset_Body = q(9);
doffset_Body = q(10);

alpha_Body = q(11);
dalpha_Body = q(12);

beta_Body = q(13);
dbeta_Body = q(14);

gamma_Body = q(15);
dgamma_Body = q(16);

coeffs_Ddr_Arm_Bottom = find_Coeffs_Ddr_Arm_Bottom(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand);
coeffs_Ddr_Arm_Bottom = [coeffs_Ddr_Arm_Bottom(:, 1:3), zeros(3,3), coeffs_Ddr_Arm_Bottom(:, 4)];

coeffs_Ddl_Arm_Bottom = find_Coeffs_Ddl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,length_Hand,m_Hand);
coeffs_Ddl_Arm_Bottom = [zeros(3,3), coeffs_Ddl_Arm_Bottom(:, 1:3), coeffs_Ddl_Arm_Bottom(:, 4)];

coeffs_Ddl_Shoulder = find_Coeffs_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,m_Body,offset_Body,width_Body);

coeffs_Ddr_Shoulder = find_Coeffs_Ddr_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,m_Body,offset_Body,width_Body);

coeffs_Matrix = [
    coeffs_Ddl_Arm_Bottom(:, 1:6) - coeffs_Ddl_Shoulder(:, 1:6);
    coeffs_Ddr_Arm_Bottom(:, 1:6) - coeffs_Ddr_Shoulder(:, 1:6);
    ];

coeffs_Target = [
    coeffs_Ddl_Arm_Bottom(:, 7) - coeffs_Ddl_Shoulder(:, 7);
    coeffs_Ddr_Arm_Bottom(:, 7) - coeffs_Ddr_Shoulder(:, 7);
    ];

f_All = inv(coeffs_Matrix) * (-coeffs_Target);

r_F_X = f_All(1);
r_F_Y = f_All(2);
r_F_Z = f_All(3);

l_F_X = f_All(4);
l_F_Y = f_All(5);
l_F_Z = f_All(6);

[ddr_Alpha_Hand,ddr_Beta_Hand] = find_Dds_Arm_R(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z);
[ddl_Alpha_Hand,ddl_Beta_Hand] = find_Dds_Arm_L(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,length_Hand,m_Hand);
[ddoffset_Body,ddalpha_Body,ddbeta_Body,ddgamma_Body] = find_Dds_Body(alpha_Body,beta_Body,dalpha_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,offset_Body,r_F_X,r_F_Y,r_F_Z,width_Body);

% ddr_Arm_Bottom = find_Ddr_Arm_Bottom(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z);
% ddl_Arm_Bottom = find_Ddl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,length_Hand,m_Hand);
% ddl_Shoulder = find_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,offset_Body,r_F_X,r_F_Y,r_F_Z,width_Body);
% ddr_Shoulder = find_Ddr_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,offset_Body,r_F_X,r_F_Y,r_F_Z,width_Body);
% 
% [ddr_Arm_Bottom - ddr_Shoulder, ddl_Arm_Bottom - ddl_Shoulder]

dotq = [dr_Alpha_Hand, ddr_Alpha_Hand, dr_Beta_Hand, ddr_Beta_Hand, dl_Alpha_Hand, ddl_Alpha_Hand, dl_Beta_Hand, ddl_Beta_Hand, ...
    doffset_Body, ddoffset_Body, dalpha_Body, ddalpha_Body, dbeta_Body, ddbeta_Body, dgamma_Body, ddgamma_Body]';

end















































