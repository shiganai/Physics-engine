function dotq = ddt(t,q,r_P_Fixed,l_P_Fixed, g, length_Hand, m_Hand, m_Body, width_Body, height_Body, depth_Body)

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

alpha_Body = q(9);
dalpha_Body = q(10);

beta_Body = q(11);
dbeta_Body = q(12);

gamma_Body = q(13);
dgamma_Body = q(14);

% x_Head = q(15);
dx_Head = q(16);
% 
% y_Head = q(17);
dy_Head = q(18);
% 
% z_Head = q(19);
dz_Head = q(20);

dl_Arm_Bottom = find_Dl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,l_Alpha_Hand,l_Beta_Hand,length_Hand);
dr_Arm_Bottom = find_Dr_Arm_Bottom(dr_Beta_Hand,dr_Alpha_Hand,length_Hand,r_Alpha_Hand,r_Beta_Hand);

dx_Head = (dl_Arm_Bottom(1) + dr_Arm_Bottom(1))/2;
dy_Head = (dl_Arm_Bottom(2) + dr_Arm_Bottom(2))/2;
dz_Head = (dl_Arm_Bottom(3) + dr_Arm_Bottom(3))/2;

% [dbeta_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp] = find_Restrained_Velocity(alpha_Body,beta_Body,dalpha_Body,dgamma_Body,dl_Beta_Hand,dl_Alpha_Hand,dr_Beta_Hand,dr_Alpha_Hand,gamma_Body,l_Alpha_Hand,l_Beta_Hand,length_Hand,r_Alpha_Hand,r_Beta_Hand,width_Body);
% if any(isnan([dbeta_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp]))
%     1 + 1
% elseif any(isinf([dbeta_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp]))
%     1 + 1
% else
%     dbeta_Body = dbeta_Body_Tmp;
%     dx_Head = dx_Head_Tmp;
%     dy_Head = dy_Head_Tmp;
%     dz_Head = dz_Head_Tmp;
% end

% [dbeta_Body_Tmp,dgamma_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp] = find_Restrained_Velocity(alpha_Body,beta_Body,dalpha_Body,dl_Beta_Hand,dl_Alpha_Hand,dr_Beta_Hand,dr_Alpha_Hand,gamma_Body,l_Alpha_Hand,l_Beta_Hand,length_Hand,r_Alpha_Hand,r_Beta_Hand,width_Body);
% if any(isnan([dbeta_Body_Tmp,dgamma_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp]))
% elseif any(isinf([dbeta_Body_Tmp,dgamma_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp]))
% else
%     dbeta_Body = dbeta_Body_Tmp;
%     dgamma_Body = dgamma_Body_Tmp;
%     dx_Head = dx_Head_Tmp;
%     dy_Head = dy_Head_Tmp;
%     dz_Head = dz_Head_Tmp;
% end

% [dalpha_Body_Tmp,dbeta_Body_Tmp,dgamma_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp] = find_Restrained_Velocity(alpha_Body,beta_Body,dl_Beta_Hand,dl_Alpha_Hand,dr_Beta_Hand,dr_Alpha_Hand,gamma_Body,l_Alpha_Hand,l_Beta_Hand,length_Hand,r_Alpha_Hand,r_Beta_Hand,width_Body);
% if any(isnan([dalpha_Body_Tmp,dbeta_Body_Tmp,dgamma_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp]))
% elseif any(isinf([dalpha_Body_Tmp,dbeta_Body_Tmp,dgamma_Body_Tmp,dx_Head_Tmp,dy_Head_Tmp,dz_Head_Tmp]))
% else
%     dalpha_Body = dalpha_Body_Tmp;
%     dbeta_Body = dbeta_Body_Tmp;
%     dgamma_Body = dgamma_Body_Tmp;
%     dx_Head = dx_Head_Tmp;
%     dy_Head = dy_Head_Tmp;
%     dz_Head = dz_Head_Tmp;
% end

coeffs_Ddr_Arm_Bottom = find_Coeffs_Ddr_Arm_Bottom(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand);
coeffs_Ddr_Arm_Bottom = [coeffs_Ddr_Arm_Bottom(:, 1:3), zeros(3,3), coeffs_Ddr_Arm_Bottom(:, 4)];

coeffs_Ddl_Arm_Bottom = find_Coeffs_Ddl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,length_Hand,m_Hand);
coeffs_Ddl_Arm_Bottom = [zeros(3,3), coeffs_Ddl_Arm_Bottom(:, 1:3), coeffs_Ddl_Arm_Bottom(:, 4)];

coeffs_Ddr_Shoulder = find_Coeffs_Ddr_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,m_Body,width_Body);
coeffs_Ddl_Shoulder = find_Coeffs_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,m_Body,width_Body);

coeffs_Matrix = [
    coeffs_Ddr_Arm_Bottom(:, 1:6) - coeffs_Ddr_Shoulder(:, 1:6);
    coeffs_Ddl_Arm_Bottom(:, 1:6) - coeffs_Ddl_Shoulder(:, 1:6);
    ];

coeffs_Target = [
    coeffs_Ddr_Arm_Bottom(:, 7) - coeffs_Ddr_Shoulder(:, 7);
    coeffs_Ddl_Arm_Bottom(:, 7) - coeffs_Ddl_Shoulder(:, 7);
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
[ddalpha_Body,ddbeta_Body,ddgamma_Body,ddx_Head,ddy_Head,ddz_Head] = find_Dds_Body(alpha_Body,beta_Body,dalpha_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,r_F_X,r_F_Y,r_F_Z,width_Body);

% ddr_Arm_Bottom = find_Ddr_Arm_Bottom(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z);
% ddl_Arm_Bottom = find_Ddl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,length_Hand,m_Hand);
% ddr_Shoulder = find_Ddr_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,r_F_X,r_F_Y,r_F_Z,width_Body);
% ddl_Shoulder = find_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,r_F_X,r_F_Y,r_F_Z,width_Body);
% 
% [ddr_Arm_Bottom - ddr_Shoulder, ddl_Arm_Bottom - ddl_Shoulder]

dd_Threshold = 3;
if abs(ddbeta_Body) < dd_Threshold
    if (dbeta_Body == 0)
        ddbeta_Body = 0;
    end
end
if abs(ddgamma_Body) < dd_Threshold
    if (dgamma_Body == 0)
        ddgamma_Body = 0;
    end
end

if abs(ddr_Beta_Hand) < dd_Threshold
    if (dr_Beta_Hand == 0)
        ddr_Beta_Hand = 0;
    end
end
if abs(ddl_Beta_Hand) < dd_Threshold
    if (dl_Beta_Hand == 0)
        ddl_Beta_Hand = 0;
    end
end

dotq = [dr_Alpha_Hand, ddr_Alpha_Hand, dr_Beta_Hand, ddr_Beta_Hand, dl_Alpha_Hand, ddl_Alpha_Hand, dl_Beta_Hand, ddl_Beta_Hand, ...
    dalpha_Body, ddalpha_Body, dbeta_Body, ddbeta_Body, dgamma_Body, ddgamma_Body, ...
    dx_Head, ddx_Head, dy_Head, ddy_Head, dz_Head, ddz_Head]';

end















































