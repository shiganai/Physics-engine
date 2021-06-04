function [ddl_Alpha_Hand,ddl_Beta_Hand] = FFD_Dds_Arm_L(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,l_Tau_Beta_Shoulder,l_Tau_Alpha_Shoulder,length_Hand,m_Hand)
%FFD_DDS_ARM_L
%    [DDL_ALPHA_HAND,DDL_BETA_HAND] = FFD_DDS_ARM_L(DL_BETA_HAND,DL_ALPHA_HAND,G,L_ALPHA_HAND,L_BETA_HAND,L_F_X,L_F_Y,L_F_Z,L_TAU_BETA_SHOULDER,L_TAU_ALPHA_SHOULDER,LENGTH_HAND,M_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 18:10:25

t2 = cos(l_Alpha_Hand);
t3 = cos(l_Beta_Hand);
t4 = sin(l_Alpha_Hand);
t5 = sin(l_Beta_Hand);
t6 = dl_Beta_Hand.^2;
t7 = dl_Alpha_Hand.^2;
t8 = length_Hand.^2;
t9 = t2.^2;
t10 = t3.^2;
t11 = t4.^2;
ddl_Alpha_Hand = ((l_Tau_Alpha_Shoulder+l_F_Z.*length_Hand.*t2.*t3-l_F_Y.*length_Hand.*t3.*t4+(g.*length_Hand.*m_Hand.*t2.*t3)./2.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t9)./2.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t11)./2.0).*-1.2e+1)./(m_Hand.*t8.*(t9.*t10.*3.0+t10.*t11.*3.0+1.0));
if nargout > 1
    t12 = t5.^2;
    ddl_Beta_Hand = ((l_Tau_Beta_Shoulder.*t9+l_Tau_Beta_Shoulder.*t11+l_F_X.*length_Hand.*t3-l_F_Y.*length_Hand.*t2.*t5-l_F_Z.*length_Hand.*t4.*t5-(g.*length_Hand.*m_Hand.*t4.*t5)./2.0-(m_Hand.*t3.*t5.*t6.*t8)./4.0+(m_Hand.*t3.*t5.*t6.*t8.*t9)./4.0+(m_Hand.*t3.*t5.*t7.*t8.*t9)./4.0+(m_Hand.*t3.*t5.*t6.*t8.*t11)./4.0+(m_Hand.*t3.*t5.*t7.*t8.*t11)./4.0).*-1.2e+1)./(m_Hand.*t8+m_Hand.*t8.*t10.*3.0+m_Hand.*t8.*t9.*t12.*3.0+m_Hand.*t8.*t11.*t12.*3.0);
end
