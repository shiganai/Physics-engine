function [A11,A21,A31] = HFD_Coeffs_Ddr_Arm_Bottom_Constant(ddr_Beta_Hand,dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_Tau_Alpha_Shoulder)
%HFD_COEFFS_DDR_ARM_BOTTOM_CONSTANT
%    [A11,A21,A31] = HFD_COEFFS_DDR_ARM_BOTTOM_CONSTANT(DDR_BETA_HAND,DR_BETA_HAND,DR_ALPHA_HAND,G,LENGTH_HAND,M_HAND,R_ALPHA_HAND,R_BETA_HAND,R_TAU_ALPHA_SHOULDER)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 18:36:46

t2 = cos(r_Alpha_Hand);
t3 = cos(r_Beta_Hand);
t4 = sin(r_Alpha_Hand);
t5 = sin(r_Beta_Hand);
t6 = dr_Beta_Hand.^2;
A11 = -ddr_Beta_Hand.*length_Hand.*t3+length_Hand.*t5.*t6;
if nargout > 1
    t7 = dr_Alpha_Hand.^2;
    t8 = length_Hand.^2;
    t9 = t2.^2;
    t10 = t3.^2;
    t11 = t4.^2;
    t12 = m_Hand.*t8;
    t13 = (g.*length_Hand.*m_Hand.*t2.*t3)./2.0;
    t14 = t9.*t10.*t12.*3.0;
    t15 = t10.*t11.*t12.*3.0;
    t16 = (dr_Beta_Hand.*dr_Alpha_Hand.*t3.*t5.*t9.*t12)./2.0;
    t17 = (dr_Beta_Hand.*dr_Alpha_Hand.*t3.*t5.*t11.*t12)./2.0;
    t18 = -t16;
    t19 = -t17;
    t20 = t12+t14+t15;
    t21 = 1.0./t20;
    t22 = r_Tau_Alpha_Shoulder+t13+t18+t19;
    A21 = -ddr_Beta_Hand.*length_Hand.*t2.*t5-length_Hand.*t2.*t3.*t6-length_Hand.*t2.*t3.*t7+dr_Beta_Hand.*dr_Alpha_Hand.*length_Hand.*t4.*t5.*2.0+length_Hand.*t3.*t4.*t21.*t22.*1.2e+1;
end
if nargout > 2
    A31 = -ddr_Beta_Hand.*length_Hand.*t4.*t5-length_Hand.*t3.*t4.*t6-length_Hand.*t3.*t4.*t7-dr_Beta_Hand.*dr_Alpha_Hand.*length_Hand.*t2.*t5.*2.0-length_Hand.*t2.*t3.*t21.*t22.*1.2e+1;
end