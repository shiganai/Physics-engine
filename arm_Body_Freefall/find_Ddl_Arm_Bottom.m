function ddl_Arm_Bottom = find_Ddl_Arm_Bottom(dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,length_Hand,m_Hand)
%FIND_DDL_ARM_BOTTOM
%    DDL_ARM_BOTTOM = FIND_DDL_ARM_BOTTOM(DL_BETA_HAND,DL_ALPHA_HAND,G,L_ALPHA_HAND,L_BETA_HAND,L_F_X,L_F_Y,L_F_Z,LENGTH_HAND,M_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-Jun-2021 11:50:37

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
t15 = l_F_Z.*length_Hand.*t2;
t16 = m_Hand.*t8;
t17 = l_F_X.*length_Hand.*t2.*t3;
t18 = l_F_Y.*length_Hand.*t2.*t5;
t19 = l_F_Y.*length_Hand.*t3.*t4;
t20 = l_F_X.*length_Hand.*t4.*t5;
t21 = (g.*length_Hand.*m_Hand.*t2)./2.0;
t22 = -t17;
t23 = -t19;
t24 = -t20;
t25 = t9.*t10.*3.0;
t26 = t9.*t16.*3.0;
t27 = t9.*t12.*3.0;
t28 = t10.*t11.*t16.*3.0;
t29 = t11.*t12.*t16.*3.0;
t30 = (t2.*t4.*t7.*t16)./4.0;
t32 = (dl_Beta_Hand.*dl_Alpha_Hand.*t2.*t4.*t10.*t16)./2.0;
t33 = (dl_Beta_Hand.*dl_Alpha_Hand.*t2.*t4.*t12.*t16)./2.0;
t34 = (t2.*t4.*t6.*t10.*t16)./4.0;
t36 = (t2.*t4.*t6.*t12.*t16)./4.0;
t31 = -t30;
t35 = t10.*t30;
t37 = t12.*t30;
t38 = t25+t27+1.0;
t40 = t16+t26+t28+t29;
t42 = t18+t22+t32+t33;
t39 = 1.0./t38;
t41 = 1.0./t40;
t43 = t15+t21+t23+t24+t31+t34+t35+t36+t37;
ddl_Arm_Bottom = [-length_Hand.*t2.*t5.*t6-length_Hand.*t2.*t5.*t7-dl_Beta_Hand.*dl_Alpha_Hand.*length_Hand.*t3.*t4.*2.0+length_Hand.*t4.*t5.*t41.*t43.*1.2e+1+t2.*t3.*t13.*t14.*t39.*t42.*1.2e+1,-length_Hand.*t2.*t3.*t6-length_Hand.*t2.*t3.*t7+dl_Beta_Hand.*dl_Alpha_Hand.*length_Hand.*t4.*t5.*2.0+length_Hand.*t3.*t4.*t41.*t43.*1.2e+1-t2.*t5.*t13.*t14.*t39.*t42.*1.2e+1,-length_Hand.*t4.*t7-length_Hand.*t2.*t41.*t43.*1.2e+1];
