function [ddr_Alpha_Hand,ddr_Beta_Hand] = FFD_Dds_Arm_R(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z,r_Tau_Beta_Shoulder,r_Tau_Alpha_Shoulder)
%FFD_DDS_ARM_R
%    [DDR_ALPHA_HAND,DDR_BETA_HAND] = FFD_DDS_ARM_R(DR_BETA_HAND,DR_ALPHA_HAND,G,LENGTH_HAND,M_HAND,R_ALPHA_HAND,R_BETA_HAND,R_F_X,R_F_Y,R_F_Z,R_TAU_BETA_SHOULDER,R_TAU_ALPHA_SHOULDER)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 18:11:21

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
ddr_Alpha_Hand = ((r_Tau_Alpha_Shoulder+length_Hand.*r_F_Z.*t2.*t3-length_Hand.*r_F_Y.*t3.*t4+(g.*length_Hand.*m_Hand.*t2.*t3)./2.0-(dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t9)./2.0-(dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t11)./2.0).*-1.2e+1)./(m_Hand.*t8.*(t9.*t10.*3.0+t10.*t11.*3.0+1.0));
if nargout > 1
    t12 = t5.^2;
    ddr_Beta_Hand = ((r_Tau_Beta_Shoulder.*t9+r_Tau_Beta_Shoulder.*t11-length_Hand.*r_F_X.*t3-length_Hand.*r_F_Y.*t2.*t5-length_Hand.*r_F_Z.*t4.*t5-(g.*length_Hand.*m_Hand.*t4.*t5)./2.0-(m_Hand.*t3.*t5.*t6.*t8)./4.0+(m_Hand.*t3.*t5.*t6.*t8.*t9)./4.0+(m_Hand.*t3.*t5.*t7.*t8.*t9)./4.0+(m_Hand.*t3.*t5.*t6.*t8.*t11)./4.0+(m_Hand.*t3.*t5.*t7.*t8.*t11)./4.0).*-1.2e+1)./(m_Hand.*t8+m_Hand.*t8.*t10.*3.0+m_Hand.*t8.*t9.*t12.*3.0+m_Hand.*t8.*t11.*t12.*3.0);
end
