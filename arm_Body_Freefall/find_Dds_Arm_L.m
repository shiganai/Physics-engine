function [ddl_Alpha_Hand,ddl_Beta_Hand] = find_Dds_Arm_L(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,length_Hand,m_Hand)
%FIND_DDS_ARM_L
%    [DDL_ALPHA_HAND,DDL_BETA_HAND] = FIND_DDS_ARM_L(DL_BETA_HAND,DL_ALPHA_HAND,G,L_ALPHA_HAND,L_BETA_HAND,L_F_X,L_F_Y,L_F_Z,LENGTH_HAND,M_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-Jun-2021 11:50:36

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
t12 = t5.^2;
ddl_Alpha_Hand = ((l_F_Z.*length_Hand.*t2+(g.*length_Hand.*m_Hand.*t2)./2.0-l_F_Y.*length_Hand.*t3.*t4-l_F_X.*length_Hand.*t4.*t5-(m_Hand.*t2.*t4.*t7.*t8)./4.0+(m_Hand.*t2.*t4.*t6.*t8.*t10)./4.0+(m_Hand.*t2.*t4.*t7.*t8.*t10)./4.0+(m_Hand.*t2.*t4.*t6.*t8.*t12)./4.0+(m_Hand.*t2.*t4.*t7.*t8.*t12)./4.0).*-1.2e+1)./(m_Hand.*t8+m_Hand.*t8.*t9.*3.0+m_Hand.*t8.*t10.*t11.*3.0+m_Hand.*t8.*t11.*t12.*3.0);
if nargout > 1
    ddl_Beta_Hand = ((-l_F_X.*length_Hand.*t2.*t3+l_F_Y.*length_Hand.*t2.*t5+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t4.*t8.*t10)./2.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t4.*t8.*t12)./2.0).*1.2e+1)./(m_Hand.*t8.*(t9.*t10.*3.0+t9.*t12.*3.0+1.0));
end
