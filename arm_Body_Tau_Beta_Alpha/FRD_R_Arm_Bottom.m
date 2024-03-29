function r_Arm_Bottom = FRD_R_Arm_Bottom(length_Hand,r_Alpha_Hand,r_Beta_Hand,r_X_Fixed,r_Y_Fixed,r_Z_Fixed)
%FRD_R_ARM_BOTTOM
%    R_ARM_BOTTOM = FRD_R_ARM_BOTTOM(LENGTH_HAND,R_ALPHA_HAND,R_BETA_HAND,R_X_FIXED,R_Y_FIXED,R_Z_FIXED)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-Jun-2021 22:17:10

t2 = cos(r_Beta_Hand);
r_Arm_Bottom = [r_X_Fixed-length_Hand.*sin(r_Beta_Hand),r_Y_Fixed+length_Hand.*t2.*cos(r_Alpha_Hand),r_Z_Fixed+length_Hand.*t2.*sin(r_Alpha_Hand)];
