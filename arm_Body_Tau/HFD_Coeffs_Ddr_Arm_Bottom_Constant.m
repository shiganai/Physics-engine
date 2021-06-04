function [A11,A21,A31] = HFD_Coeffs_Ddr_Arm_Bottom_Constant(ddr_Beta_Hand,dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_Tau_Alpha_Shoulder)
%HFD_COEFFS_DDR_ARM_BOTTOM_CONSTANT
%    [A11,A21,A31] = HFD_COEFFS_DDR_ARM_BOTTOM_CONSTANT(DDR_BETA_HAND,DR_BETA_HAND,DR_ALPHA_HAND,G,LENGTH_HAND,M_HAND,R_ALPHA_HAND,R_BETA_HAND,R_TAU_ALPHA_SHOULDER)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 16:54:22

t2 = cos(r_Alpha_Hand);
t3 = cos(r_Beta_Hand);
t4 = sin(r_Alpha_Hand);
t5 = sin(r_Beta_Hand);
t6 = dr_Beta_Hand.^2;
t7 = dr_Alpha_Hand.^2;
t8 = length_Hand.^2;
t9 = t2.^2;
t10 = t3.^2;
t11 = t4.^2;
t12 = t5.^2;
t13 = m_Hand.*t8;
t16 = (g.*length_Hand.*m_Hand.*t2)./2.0;
t14 = r_Tau_Alpha_Shoulder.*t10;
t15 = r_Tau_Alpha_Shoulder.*t12;
t17 = t9.*t13.*3.0;
t18 = t10.*t11.*t13.*3.0;
t19 = t11.*t12.*t13.*3.0;
t20 = (t2.*t4.*t7.*t13)./4.0;
t22 = (t2.*t4.*t6.*t10.*t13)./4.0;
t24 = (t2.*t4.*t6.*t12.*t13)./4.0;
t21 = -t20;
t23 = t10.*t20;
t25 = t12.*t20;
t26 = t13+t17+t18+t19;
t27 = 1.0./t26;
t28 = t14+t15+t16+t21+t22+t23+t24+t25;
A11 = -ddr_Beta_Hand.*length_Hand.*t2.*t3+length_Hand.*t2.*t5.*t6+length_Hand.*t2.*t5.*t7+dr_Beta_Hand.*dr_Alpha_Hand.*length_Hand.*t3.*t4.*2.0-length_Hand.*t4.*t5.*t27.*t28.*1.2e+1;
if nargout > 1
    A21 = -ddr_Beta_Hand.*length_Hand.*t2.*t5-length_Hand.*t2.*t3.*t6-length_Hand.*t2.*t3.*t7+dr_Beta_Hand.*dr_Alpha_Hand.*length_Hand.*t4.*t5.*2.0+length_Hand.*t3.*t4.*t27.*t28.*1.2e+1;
end
if nargout > 2
    A31 = -length_Hand.*t4.*t7-length_Hand.*t2.*t27.*t28.*1.2e+1;
end
