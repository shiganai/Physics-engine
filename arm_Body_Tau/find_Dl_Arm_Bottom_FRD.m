function dl_Arm_Bottom = find_Dl_Arm_Bottom_FRD(dl_Beta_Hand,dl_Alpha_Hand,l_Alpha_Hand,l_Beta_Hand,length_Hand)
%FIND_DL_ARM_BOTTOM_FRD
%    DL_ARM_BOTTOM = FIND_DL_ARM_BOTTOM_FRD(DL_BETA_HAND,DL_ALPHA_HAND,L_ALPHA_HAND,L_BETA_HAND,LENGTH_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 15:28:09

t2 = cos(l_Alpha_Hand);
t3 = cos(l_Beta_Hand);
t4 = sin(l_Alpha_Hand);
t5 = sin(l_Beta_Hand);
dl_Arm_Bottom = [dl_Beta_Hand.*length_Hand.*t2.*t3-dl_Alpha_Hand.*length_Hand.*t4.*t5,-dl_Beta_Hand.*length_Hand.*t2.*t5-dl_Alpha_Hand.*length_Hand.*t3.*t4,dl_Alpha_Hand.*length_Hand.*t2];
