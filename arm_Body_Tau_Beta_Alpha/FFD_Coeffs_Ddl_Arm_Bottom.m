function [A11,A12,A13,A14,A21,A22,A23,A24,A31,A32,A33,A34] = FFD_Coeffs_Ddl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_Tau_Beta_Shoulder,l_Tau_Alpha_Shoulder,length_Hand,m_Hand)
%FFD_COEFFS_DDL_ARM_BOTTOM
%    [A11,A12,A13,A14,A21,A22,A23,A24,A31,A32,A33,A34] = FFD_COEFFS_DDL_ARM_BOTTOM(DL_BETA_HAND,DL_ALPHA_HAND,G,L_ALPHA_HAND,L_BETA_HAND,L_TAU_BETA_SHOULDER,L_TAU_ALPHA_SHOULDER,LENGTH_HAND,M_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 18:10:26

t2 = cos(l_Alpha_Hand);
t3 = cos(l_Beta_Hand);
t4 = sin(l_Alpha_Hand);
t5 = sin(l_Beta_Hand);
t6 = dl_Beta_Hand.^2;
t7 = dl_Alpha_Hand.^2;
t8 = length_Hand.^2;
t13 = 1.0./length_Hand;
t14 = 1.0./m_Hand;
t9 = t2.^2;
t10 = t3.^2;
t11 = t4.^2;
t12 = t5.^2;
t15 = m_Hand.*t8;
t21 = (g.*length_Hand.*m_Hand.*t2.*t3)./2.0;
t22 = (g.*length_Hand.*m_Hand.*t4.*t5)./2.0;
t16 = l_Tau_Beta_Shoulder.*t11;
t17 = l_Tau_Beta_Shoulder.*t9;
t18 = t9.*t10.*3.0;
t19 = t10.*t15.*3.0;
t20 = t10.*t11.*3.0;
t23 = -t22;
t24 = t9.*t12.*t15.*3.0;
t25 = t11.*t12.*t15.*3.0;
t26 = (t3.*t5.*t6.*t15)./4.0;
t28 = (dl_Beta_Hand.*dl_Alpha_Hand.*t3.*t5.*t9.*t15)./2.0;
t29 = (dl_Beta_Hand.*dl_Alpha_Hand.*t3.*t5.*t11.*t15)./2.0;
t31 = (t3.*t5.*t7.*t9.*t15)./4.0;
t34 = (t3.*t5.*t7.*t11.*t15)./4.0;
t27 = -t26;
t30 = t9.*t26;
t32 = -t28;
t33 = t11.*t26;
t35 = -t29;
t36 = t18+t20+1.0;
t39 = t15+t19+t24+t25;
t37 = 1.0./t36;
t40 = 1.0./t39;
A11 = t8.*t10.*t40.*-1.2e+1;
if nargout > 1
    t41 = l_Tau_Alpha_Shoulder+t21+t32+t35;
    t47 = t16+t17+t23+t27+t30+t31+t33+t34;
    t38 = t2.*t4.*t10.*t14.*t37.*1.2e+1;
    t42 = t2.*t3.*t5.*t8.*t40.*1.2e+1;
    A12 = t42;
end
if nargout > 2
    t43 = t3.*t4.*t5.*t8.*t40.*1.2e+1;
    A13 = t43;
end
if nargout > 3
    A14 = -length_Hand.*t5.*t6-length_Hand.*t3.*t40.*t47.*1.2e+1;
end
if nargout > 4
    A21 = t42;
end
if nargout > 5
    A22 = t8.*t9.*t12.*t40.*-1.2e+1-t10.*t11.*t14.*t37.*1.2e+1;
end
if nargout > 6
    t44 = t2.*t4.*t8.*t12.*t40.*1.2e+1;
    t45 = -t44;
    t46 = t38+t45;
    A23 = t46;
end
if nargout > 7
    A24 = -length_Hand.*t2.*t3.*t6-length_Hand.*t2.*t3.*t7+dl_Beta_Hand.*dl_Alpha_Hand.*length_Hand.*t4.*t5.*2.0+length_Hand.*t2.*t5.*t40.*t47.*1.2e+1+t3.*t4.*t13.*t14.*t37.*t41.*1.2e+1;
end
if nargout > 8
    A31 = t43;
end
if nargout > 9
    A32 = t46;
end
if nargout > 10
    A33 = t9.*t10.*t14.*t37.*-1.2e+1-t8.*t11.*t12.*t40.*1.2e+1;
end
if nargout > 11
    A34 = -length_Hand.*t3.*t4.*t6-length_Hand.*t3.*t4.*t7-dl_Beta_Hand.*dl_Alpha_Hand.*length_Hand.*t2.*t5.*2.0+length_Hand.*t4.*t5.*t40.*t47.*1.2e+1-t2.*t3.*t13.*t14.*t37.*t41.*1.2e+1;
end