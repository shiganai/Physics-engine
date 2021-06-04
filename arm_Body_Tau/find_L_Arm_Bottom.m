function l_Arm_Bottom = find_L_Arm_Bottom(l_Alpha_Hand,l_Beta_Hand,l_X_Fixed,l_Y_Fixed,l_Z_Fixed,length_Hand)
%FIND_L_ARM_BOTTOM
%    L_ARM_BOTTOM = FIND_L_ARM_BOTTOM(L_ALPHA_HAND,L_BETA_HAND,L_X_FIXED,L_Y_FIXED,L_Z_FIXED,LENGTH_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-Jun-2021 23:36:47

t2 = cos(l_Alpha_Hand);
l_Arm_Bottom = [l_X_Fixed+length_Hand.*t2.*sin(l_Beta_Hand),l_Y_Fixed+length_Hand.*t2.*cos(l_Beta_Hand),l_Z_Fixed+length_Hand.*sin(l_Alpha_Hand)];
