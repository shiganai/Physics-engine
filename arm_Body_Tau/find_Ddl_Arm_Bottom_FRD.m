function ddl_Arm_Bottom = find_Ddl_Arm_Bottom_FRD(ddl_Beta_Hand,ddl_Alpha_Hand,dl_Beta_Hand,dl_Alpha_Hand,l_Alpha_Hand,l_Beta_Hand,length_Hand)
%FIND_DDL_ARM_BOTTOM_FRD
%    DDL_ARM_BOTTOM = FIND_DDL_ARM_BOTTOM_FRD(DDL_BETA_HAND,DDL_ALPHA_HAND,DL_BETA_HAND,DL_ALPHA_HAND,L_ALPHA_HAND,L_BETA_HAND,LENGTH_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 15:28:08

t2 = cos(l_Alpha_Hand);
t3 = cos(l_Beta_Hand);
t4 = sin(l_Alpha_Hand);
t5 = sin(l_Beta_Hand);
t6 = dl_Beta_Hand.^2;
t7 = dl_Alpha_Hand.^2;
ddl_Arm_Bottom = [ddl_Beta_Hand.*length_Hand.*t2.*t3-ddl_Alpha_Hand.*length_Hand.*t4.*t5-length_Hand.*t2.*t5.*t6-length_Hand.*t2.*t5.*t7-dl_Beta_Hand.*dl_Alpha_Hand.*length_Hand.*t3.*t4.*2.0,-ddl_Beta_Hand.*length_Hand.*t2.*t5-ddl_Alpha_Hand.*length_Hand.*t3.*t4-length_Hand.*t2.*t3.*t6-length_Hand.*t2.*t3.*t7+dl_Beta_Hand.*dl_Alpha_Hand.*length_Hand.*t4.*t5.*2.0,ddl_Alpha_Hand.*length_Hand.*t2-length_Hand.*t4.*t7];