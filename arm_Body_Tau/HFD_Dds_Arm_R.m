function [ddr_Alpha_Hand,r_Tau_Beta_Shoulder] = HFD_Dds_Arm_R(ddr_Beta_Hand,dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z,r_Tau_Alpha_Shoulder)
%HFD_DDS_ARM_R
%    [DDR_ALPHA_HAND,R_TAU_BETA_SHOULDER] = HFD_DDS_ARM_R(DDR_BETA_HAND,DR_BETA_HAND,DR_ALPHA_HAND,G,LENGTH_HAND,M_HAND,R_ALPHA_HAND,R_BETA_HAND,R_F_X,R_F_Y,R_F_Z,R_TAU_ALPHA_SHOULDER)

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
ddr_Alpha_Hand = ((r_Tau_Alpha_Shoulder.*t10+r_Tau_Alpha_Shoulder.*t12+length_Hand.*r_F_Z.*t2+(g.*length_Hand.*m_Hand.*t2)./2.0-length_Hand.*r_F_Y.*t3.*t4+length_Hand.*r_F_X.*t4.*t5-(m_Hand.*t2.*t4.*t7.*t8)./4.0+(m_Hand.*t2.*t4.*t6.*t8.*t10)./4.0+(m_Hand.*t2.*t4.*t7.*t8.*t10)./4.0+(m_Hand.*t2.*t4.*t6.*t8.*t12)./4.0+(m_Hand.*t2.*t4.*t7.*t8.*t12)./4.0).*-1.2e+1)./(m_Hand.*t8+m_Hand.*t8.*t9.*3.0+m_Hand.*t8.*t10.*t11.*3.0+m_Hand.*t8.*t11.*t12.*3.0);
if nargout > 1
    r_Tau_Beta_Shoulder = ddr_Beta_Hand.*m_Hand.*t8.*(-1.0./1.2e+1)+length_Hand.*r_F_X.*t2.*t3+length_Hand.*r_F_Y.*t2.*t5-(ddr_Beta_Hand.*m_Hand.*t8.*t9.*t10)./4.0-(ddr_Beta_Hand.*m_Hand.*t8.*t9.*t12)./4.0+(dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t2.*t4.*t8.*t10)./2.0+(dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t2.*t4.*t8.*t12)./2.0;
end
