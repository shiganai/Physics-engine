function dr_Arm_Bottom = find_Dr_Arm_Bottom_FRD(dr_Beta_Hand,dr_Alpha_Hand,length_Hand,r_Alpha_Hand,r_Beta_Hand)
%FIND_DR_ARM_BOTTOM_FRD
%    DR_ARM_BOTTOM = FIND_DR_ARM_BOTTOM_FRD(DR_BETA_HAND,DR_ALPHA_HAND,LENGTH_HAND,R_ALPHA_HAND,R_BETA_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 15:28:23

t2 = cos(r_Alpha_Hand);
t3 = cos(r_Beta_Hand);
t4 = sin(r_Alpha_Hand);
t5 = sin(r_Beta_Hand);
dr_Arm_Bottom = [-dr_Beta_Hand.*length_Hand.*t2.*t3+dr_Alpha_Hand.*length_Hand.*t4.*t5,-dr_Beta_Hand.*length_Hand.*t2.*t5-dr_Alpha_Hand.*length_Hand.*t3.*t4,dr_Alpha_Hand.*length_Hand.*t2];
