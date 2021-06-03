function r_Arm_Bottom = find_R_Arm_Bottom(length_Hand,r_Alpha_Hand,r_Beta_Hand,r_X_Fixed,r_Y_Fixed,r_Z_Fixed)
%FIND_R_ARM_BOTTOM
%    R_ARM_BOTTOM = FIND_R_ARM_BOTTOM(LENGTH_HAND,R_ALPHA_HAND,R_BETA_HAND,R_X_FIXED,R_Y_FIXED,R_Z_FIXED)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-Jun-2021 11:51:00

t2 = cos(r_Alpha_Hand);
r_Arm_Bottom = [r_X_Fixed-length_Hand.*t2.*sin(r_Beta_Hand),r_Y_Fixed+length_Hand.*t2.*cos(r_Beta_Hand),r_Z_Fixed+length_Hand.*sin(r_Alpha_Hand)];