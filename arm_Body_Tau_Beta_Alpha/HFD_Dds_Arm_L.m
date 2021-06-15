function [ddl_Alpha_Hand,l_Tau_Beta_Shoulder] = HFD_Dds_Arm_L(ddl_Beta_Hand,dl_Beta_Hand,dl_Alpha_Hand,g,l_Alpha_Hand,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,l_Tau_Alpha_Shoulder,length_Hand,m_Hand)
%HFD_DDS_ARM_L
%    [DDL_ALPHA_HAND,L_TAU_BETA_SHOULDER] = HFD_DDS_ARM_L(DDL_BETA_HAND,DL_BETA_HAND,DL_ALPHA_HAND,G,L_ALPHA_HAND,L_BETA_HAND,L_F_X,L_F_Y,L_F_Z,L_TAU_ALPHA_SHOULDER,LENGTH_HAND,M_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    15-Jun-2021 15:05:36

t2 = cos(l_Alpha_Hand);
t3 = cos(l_Beta_Hand);
t4 = sin(l_Alpha_Hand);
t5 = sin(l_Beta_Hand);
t6 = dl_Beta_Hand.^2;
t7 = dl_Alpha_Hand.^2;
t8 = length_Hand.^2;
t25 = -l_Tau_Alpha_Shoulder;
t9 = t2.^2;
t10 = t2.^3;
t12 = t3.^2;
t13 = t3.^3;
t15 = t3.^5;
t17 = t4.^2;
t18 = t4.^3;
t20 = t5.^2;
t21 = t5.^3;
t23 = t5.^5;
t26 = l_F_Z.*length_Hand.*t2.*t3;
t27 = l_F_Y.*length_Hand.*t3.*t4;
t30 = (g.*length_Hand.*m_Hand.*t2.*t3)./2.0;
t11 = t9.^2;
t14 = t12.^2;
t16 = t12.^3;
t19 = t17.^2;
t22 = t20.^2;
t24 = t20.^3;
t28 = t9+t17;
t29 = -t26;
t32 = -t30;
t34 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t9)./2.0;
t35 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t17)./2.0;
t44 = (ddl_Beta_Hand.*m_Hand.*t2.*t3.*t5.*t8.*t18)./1.2e+1;
t46 = (m_Hand.*t2.*t7.*t8.*t12.*t18)./6.0;
t50 = (m_Hand.*t2.*t6.*t8.*t18.*t20)./6.0;
t56 = (ddl_Beta_Hand.*m_Hand.*t2.*t3.*t8.*t18.*t21)./1.2e+1;
t57 = (ddl_Beta_Hand.*m_Hand.*t2.*t5.*t8.*t13.*t18)./1.2e+1;
t58 = (ddl_Beta_Hand.*m_Hand.*t3.*t4.*t8.*t10.*t21)./1.2e+1;
t59 = (ddl_Beta_Hand.*m_Hand.*t4.*t5.*t8.*t10.*t13)./1.2e+1;
t60 = (ddl_Beta_Hand.*m_Hand.*t3.*t4.*t8.*t10.*t23)./1.2e+1;
t61 = (ddl_Beta_Hand.*m_Hand.*t4.*t5.*t8.*t10.*t15)./1.2e+1;
t63 = (m_Hand.*t3.*t6.*t8.*t9.*t17)./8.0;
t64 = (m_Hand.*t2.*t6.*t8.*t12.*t18)./1.2e+1;
t72 = (ddl_Beta_Hand.*m_Hand.*t8.*t9.*t17.*t21)./1.2e+1;
t75 = (m_Hand.*t4.*t6.*t8.*t10.*t20)./1.2e+1;
t83 = (ddl_Beta_Hand.*m_Hand.*t4.*t8.*t10.*t13.*t21)./6.0;
t90 = (m_Hand.*t6.*t8.*t9.*t13.*t17)./3.0;
t93 = (ddl_Beta_Hand.*m_Hand.*t5.*t8.*t9.*t12.*t17)./1.2e+1;
t94 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t17.*t21)./3.0;
t95 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t5.*t8.*t9.*t13.*t17)./3.0;
t96 = (m_Hand.*t6.*t8.*t9.*t15.*t17)./8.0;
t97 = (m_Hand.*t3.*t6.*t8.*t9.*t17.*t20)./3.0;
t99 = (m_Hand.*t4.*t6.*t8.*t10.*t12.*t20)./6.0;
t104 = (m_Hand.*t4.*t7.*t8.*t10.*t12.*t20)./6.0;
t112 = (m_Hand.*t6.*t8.*t9.*t13.*t17.*t20)./4.0;
t31 = 1.0./t28;
t33 = (ddl_Beta_Hand.*m_Hand.*t5.*t8.*t19)./2.4e+1;
t36 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t19)./6.0;
t37 = (ddl_Beta_Hand.*m_Hand.*t8.*t11.*t23)./2.4e+1;
t38 = (m_Hand.*t3.*t6.*t8.*t19)./1.2e+1;
t39 = (ddl_Beta_Hand.*m_Hand.*t5.*t8.*t11.*t14)./2.4e+1;
t40 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t11.*t23)./6.0;
t41 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t5.*t8.*t11.*t15)./6.0;
t42 = (m_Hand.*t6.*t8.*t11.*t15)./1.2e+1;
t43 = (m_Hand.*t6.*t8.*t11.*t13)./2.4e+1;
t45 = (m_Hand.*t6.*t8.*t13.*t19)./2.4e+1;
t47 = (m_Hand.*t2.*t7.*t8.*t14.*t18)./6.0;
t48 = (m_Hand.*t4.*t7.*t8.*t10.*t14)./6.0;
t49 = (m_Hand.*t4.*t7.*t8.*t10.*t16)./6.0;
t51 = (m_Hand.*t2.*t6.*t8.*t18.*t22)./4.0;
t52 = (m_Hand.*t4.*t6.*t8.*t10.*t22)./4.0;
t53 = (m_Hand.*t4.*t6.*t8.*t10.*t24)./6.0;
t62 = (ddl_Beta_Hand.*m_Hand.*t8.*t11.*t12.*t21)./1.2e+1;
t65 = (m_Hand.*t2.*t6.*t8.*t14.*t18)./1.2e+1;
t66 = (m_Hand.*t4.*t6.*t8.*t10.*t14)./1.2e+1;
t67 = (m_Hand.*t4.*t6.*t8.*t10.*t16)./1.2e+1;
t68 = -t46;
t70 = (m_Hand.*t3.*t6.*t8.*t11.*t22)./1.2e+1;
t71 = (m_Hand.*t3.*t6.*t8.*t11.*t20)./2.4e+1;
t73 = -t50;
t76 = (m_Hand.*t2.*t6.*t8.*t18.*t24)./1.2e+1;
t77 = (m_Hand.*t3.*t6.*t8.*t19.*t20)./2.4e+1;
t78 = (dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t11.*t13.*t21)./3.0;
t79 = -t56;
t80 = -t57;
t81 = -t60;
t82 = -t61;
t84 = -t63;
t91 = (m_Hand.*t6.*t8.*t11.*t13.*t20)./6.0;
t92 = -t83;
t98 = t12.*t50;
t101 = (m_Hand.*t2.*t6.*t8.*t12.*t18.*t22)./6.0;
t102 = (m_Hand.*t4.*t7.*t8.*t10.*t14.*t20)./3.0;
t103 = t20.*t46;
t105 = (m_Hand.*t4.*t7.*t8.*t10.*t12.*t22)./6.0;
t106 = -t96;
t107 = -t99;
t109 = t22.*t63;
t111 = -t104;
t113 = m_Hand.*t3.*t6.*t8.*t9.*t17.*t22.*(-1.0./8.0);
t114 = m_Hand.*t2.*t6.*t8.*t14.*t18.*t20.*(-1.0./1.2e+1);
t115 = -t112;
t54 = -t43;
t55 = -t45;
t69 = -t48;
t74 = -t52;
t85 = -t65;
t86 = -t67;
t87 = -t71;
t88 = -t76;
t89 = -t77;
t100 = t12.*t52;
t108 = -t101;
t110 = t20.*t65;
t116 = t25+t27+t29+t32+t33+t34+t35+t36+t37+t38+t39+t40+t41+t42+t44+t47+t49+t51+t53+t54+t55+t58+t59+t62+t64+t66+t68+t69+t70+t72+t73+t74+t75+t78+t79+t80+t81+t82+t84+t85+t86+t87+t88+t89+t90+t91+t92+t93+t94+t95+t97+t98+t100+t102+t103+t105+t106+t107+t108+t111+t113+t114+t115;
ddl_Alpha_Hand = (t116.*1.2e+1)./(m_Hand.*t8.*t9.*t12.*3.0+m_Hand.*t8.*t11.*t16+m_Hand.*t8.*t12.*t17.*3.0+m_Hand.*t8.*t12.*t19+m_Hand.*t8.*t9.*t14.*t17.*2.0+m_Hand.*t8.*t11.*t12.*t22+m_Hand.*t8.*t11.*t14.*t20.*2.0+m_Hand.*t8.*t9.*t12.*t17.*t20.*2.0);
if nargout > 1
    l_Tau_Beta_Shoulder = -t31.*(t5.*t46+t5.*t48+l_F_X.*length_Hand.*t3+(ddl_Beta_Hand.*m_Hand.*t8.*t12)./4.0-l_F_Y.*length_Hand.*t2.*t5-l_F_Z.*length_Hand.*t4.*t5-(g.*length_Hand.*m_Hand.*t4.*t5)./2.0+(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t20)./4.0+(ddl_Beta_Hand.*m_Hand.*t8.*t17.*t20)./4.0-(m_Hand.*t3.*t5.*t6.*t8)./4.0-(ddl_Beta_Hand.*m_Hand.*t2.*t3.*t8.*t18)./1.2e+1-(ddl_Beta_Hand.*m_Hand.*t4.*t8.*t10.*t13)./1.2e+1+(ddl_Beta_Hand.*m_Hand.*t4.*t8.*t10.*t15)./1.2e+1+(ddl_Beta_Hand.*m_Hand.*t2.*t8.*t13.*t18)./1.2e+1+(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t17.*t20)./1.2e+1-(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t17.*t22)./6.0+(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t17.*t24)./1.2e+1+(m_Hand.*t3.*t5.*t6.*t8.*t9)./4.0+(m_Hand.*t3.*t5.*t7.*t8.*t9)./4.0+(m_Hand.*t2.*t5.*t6.*t8.*t18)./2.4e+1+(m_Hand.*t3.*t5.*t6.*t8.*t17)./4.0-(m_Hand.*t2.*t5.*t7.*t8.*t18)./6.0+(m_Hand.*t3.*t5.*t7.*t8.*t17)./4.0+(m_Hand.*t3.*t5.*t7.*t8.*t19)./6.0-(m_Hand.*t5.*t7.*t8.*t11.*t13)./1.2e+1+(m_Hand.*t5.*t7.*t8.*t11.*t15)./6.0+(m_Hand.*t4.*t6.*t8.*t10.*t21)./2.4e+1-(m_Hand.*t3.*t7.*t8.*t11.*t21)./1.2e+1-(m_Hand.*t4.*t7.*t8.*t10.*t21)./6.0-(m_Hand.*t4.*t6.*t8.*t10.*t23)./2.4e+1+(m_Hand.*t3.*t7.*t8.*t11.*t23)./6.0+(m_Hand.*t4.*t7.*t8.*t10.*t23)./6.0-(m_Hand.*t5.*t7.*t8.*t13.*t19)./1.2e+1-(m_Hand.*t2.*t6.*t8.*t18.*t21)./2.4e+1+(m_Hand.*t2.*t7.*t8.*t18.*t21)./6.0-(m_Hand.*t3.*t7.*t8.*t19.*t21)./1.2e+1+(m_Hand.*t7.*t8.*t11.*t13.*t21)./3.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t19)./1.2e+1-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t11.*t13)./1.2e+1+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t11.*t15)./1.2e+1-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t13.*t19)./1.2e+1-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t17)./4.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t11.*t20)./1.2e+1+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t20)./6.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t11.*t22)./1.2e+1-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t22)./3.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t24)./6.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t9.*t13.*t17)./2.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t8.*t18.*t20)./6.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t9.*t15.*t17)./4.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t8.*t18.*t22)./3.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t19.*t20)./1.2e+1-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t8.*t18.*t24)./6.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t11.*t13.*t20)./6.0-(ddl_Beta_Hand.*m_Hand.*t3.*t4.*t8.*t10.*t20)./1.2e+1+(ddl_Beta_Hand.*m_Hand.*t3.*t4.*t8.*t10.*t22)./1.2e+1+(ddl_Beta_Hand.*m_Hand.*t2.*t3.*t8.*t18.*t20)./1.2e+1+(ddl_Beta_Hand.*m_Hand.*t4.*t8.*t10.*t13.*t20)./6.0-(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t12.*t17.*t20)./6.0+(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t12.*t17.*t22)./6.0+(ddl_Beta_Hand.*m_Hand.*t8.*t9.*t14.*t17.*t20)./1.2e+1+(m_Hand.*t4.*t5.*t6.*t8.*t10.*t12)./2.4e+1-(m_Hand.*t4.*t5.*t7.*t8.*t10.*t12)./6.0-(m_Hand.*t4.*t5.*t6.*t8.*t10.*t14)./2.4e+1+(m_Hand.*t3.*t5.*t6.*t8.*t9.*t17)./1.2e+1-(m_Hand.*t3.*t5.*t7.*t8.*t9.*t17)./4.0-(m_Hand.*t2.*t5.*t6.*t8.*t12.*t18)./2.4e+1-(m_Hand.*t5.*t6.*t8.*t9.*t13.*t17)./6.0+m_Hand.*t5.*t7.*t8.*t9.*t13.*t17.*(2.0./3.0)+(m_Hand.*t5.*t6.*t8.*t9.*t15.*t17)./1.2e+1-(m_Hand.*t4.*t6.*t8.*t10.*t12.*t21)./1.2e+1-(m_Hand.*t5.*t7.*t8.*t9.*t15.*t17)./4.0+(m_Hand.*t4.*t7.*t8.*t10.*t12.*t21)./3.0-(m_Hand.*t3.*t6.*t8.*t9.*t17.*t21)./6.0+m_Hand.*t3.*t7.*t8.*t9.*t17.*t21.*(2.0./3.0)+(m_Hand.*t3.*t6.*t8.*t9.*t17.*t23)./1.2e+1-(m_Hand.*t3.*t7.*t8.*t9.*t17.*t23)./4.0+(m_Hand.*t6.*t8.*t9.*t13.*t17.*t21)./6.0-(m_Hand.*t7.*t8.*t9.*t13.*t17.*t21)./2.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t12.*t20)./3.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t12.*t22)./3.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t14.*t20)./6.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t17.*t20)./2.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t17.*t22)./4.0+(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t8.*t12.*t18.*t20)./3.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t8.*t12.*t18.*t22)./3.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t2.*t8.*t14.*t18.*t20)./6.0-(dl_Beta_Hand.*dl_Alpha_Hand.*m_Hand.*t8.*t9.*t13.*t17.*t20)./2.0)+(t31.*t116.*(t5.*t19+t11.*t23+t5.*t11.*t14+t11.*t12.*t21.*2.0+t9.*t17.*t21.*2.0+t2.*t3.*t5.*t18.*2.0+t4.*t5.*t10.*t13.*2.0-t4.*t5.*t10.*t15.*2.0-t2.*t5.*t13.*t18.*2.0+t3.*t4.*t10.*t21.*2.0-t3.*t4.*t10.*t23.*2.0+t5.*t9.*t12.*t17.*2.0-t2.*t3.*t18.*t21.*2.0-t4.*t10.*t13.*t21.*4.0))./(t9.*t12.*6.0+t11.*t16.*2.0+t12.*t17.*6.0+t12.*t19.*2.0+t9.*t14.*t17.*4.0+t11.*t12.*t22.*2.0+t11.*t14.*t20.*4.0+t9.*t12.*t17.*t20.*4.0);
end
