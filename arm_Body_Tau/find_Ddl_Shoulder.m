function ddl_Shoulder = find_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,l_Tau_Alpha_Shoulder,m_Body,r_F_X,r_F_Y,r_F_Z,r_Tau_Alpha_Shoulder,width_Body)
%FIND_DDL_SHOULDER
%    DDL_SHOULDER = FIND_DDL_SHOULDER(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DBETA_BODY,DEPTH_BODY,DGAMMA_BODY,G,GAMMA_BODY,HEIGHT_BODY,L_F_X,L_F_Y,L_F_Z,L_TAU_ALPHA_SHOULDER,M_BODY,R_F_X,R_F_Y,R_F_Z,R_TAU_ALPHA_SHOULDER,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-Jun-2021 23:36:40

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = sin(alpha_Body);
t6 = sin(beta_Body);
t7 = sin(gamma_Body);
t8 = g.*m_Body;
t9 = dalpha_Body.^2;
t10 = dbeta_Body.^2;
t11 = depth_Body.^2;
t12 = dgamma_Body.^2;
t13 = height_Body.^2;
t15 = width_Body.^2;
t20 = 1.0./m_Body;
t14 = t13.^2;
t16 = t2.^2;
t17 = t4.^2;
t18 = t5.^2;
t19 = t7.^2;
t21 = t5.*t6;
t22 = -t8;
t23 = m_Body.*t11;
t24 = m_Body.*t15;
t25 = t2.*t3;
t26 = t2.*t6;
t27 = t3.*t5;
t32 = t11.*t13;
t33 = t11+t13;
t34 = t11.*t15;
t35 = t13.*t15;
t36 = t13+t15;
t37 = dalpha_Body.*dgamma_Body.*height_Body.*t5.*t7.*4.0;
t38 = dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*t2.*t7;
t41 = height_Body.*t2.*t4.*t9.*2.0;
t42 = height_Body.*t2.*t4.*t12.*2.0;
t46 = (l_F_X.*t3.*t7.*width_Body)./2.0;
t47 = (l_F_X.*t4.*t6.*width_Body)./2.0;
t48 = (r_F_X.*t3.*t7.*width_Body)./2.0;
t49 = (r_F_X.*t4.*t6.*width_Body)./2.0;
t50 = (height_Body.*t2.*t4.*t8)./2.0;
t51 = (height_Body.*t5.*t7.*t8)./2.0;
t53 = (height_Body.*m_Body.*t7.*t12)./2.0;
t67 = (height_Body.*m_Body.*t4.*t5.*t9)./2.0;
t68 = (height_Body.*m_Body.*t4.*t5.*t12)./2.0;
t78 = (m_Body.*t4.*t7.*t12.*t13)./4.0;
t28 = t7.*t25;
t29 = t7.*t26;
t30 = t7.*t27;
t31 = t7.*t21;
t43 = -t37;
t44 = 1.0./t33;
t45 = 1.0./t36;
t52 = t23+t24;
t54 = -t46;
t55 = -t47;
t56 = -t50;
t57 = -t51;
t58 = (l_F_Y.*t4.*t25.*width_Body)./2.0;
t59 = (r_F_Y.*t4.*t25.*width_Body)./2.0;
t60 = (l_F_Z.*t4.*t27.*width_Body)./2.0;
t61 = (r_F_Z.*t4.*t27.*width_Body)./2.0;
t62 = -t53;
t63 = t2.*t5.*t14.*t17;
t64 = t2.*t5.*t14.*t19;
t76 = t2.*t5.*t19.*t32;
t77 = t2.*t5.*t17.*t35;
t79 = -t78;
t81 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t7.*t13.*t16)./2.0;
t82 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t7.*t13.*t18)./2.0;
t83 = (m_Body.*t4.*t7.*t9.*t13.*t16)./4.0;
t84 = t16.*t78;
t85 = (m_Body.*t4.*t7.*t9.*t13.*t18)./4.0;
t86 = t18.*t78;
t99 = t14+t32+t34+t35;
t105 = l_F_Z+r_F_Z+t22+t38+t67+t68;
t39 = -t29;
t40 = -t30;
t65 = t21+t28;
t66 = t25+t31;
t69 = -t59;
t70 = -t61;
t71 = 1.0./t52;
t72 = -t64;
t73 = l_F_X+r_F_X+t62;
t80 = -t76;
t100 = 1.0./t99;
t101 = t41+t42+t43;
t109 = height_Body.*t2.*t4.*t20.*t44.*t105.*6.0;
t110 = height_Body.*t5.*t7.*t20.*t45.*t105.*6.0;
t74 = t26+t40;
t75 = t27+t39;
t87 = (l_F_Z.*t65.*width_Body)./2.0;
t88 = (l_F_Z.*t66.*width_Body)./2.0;
t89 = (r_F_Z.*t65.*width_Body)./2.0;
t90 = (r_F_Z.*t66.*width_Body)./2.0;
t102 = height_Body.*t4.*t20.*t45.*t73.*6.0;
t103 = (m_Body.*t101)./4.0;
t106 = t63+t72+t77+t80;
t111 = -t109;
t115 = t48+t54+t57+t58+t60+t69+t70+t79+t83+t84+t85+t86;
t91 = -t87;
t92 = -t88;
t93 = (l_F_Y.*t74.*width_Body)./2.0;
t94 = (l_F_Y.*t75.*width_Body)./2.0;
t95 = (r_F_Y.*t74.*width_Body)./2.0;
t96 = (r_F_Y.*t75.*width_Body)./2.0;
t104 = l_F_Y+r_F_Y+t103;
t116 = t20.*t45.*t115.*1.2e+1;
t97 = -t93;
t98 = -t96;
t107 = height_Body.*t4.*t5.*t20.*t44.*t104.*6.0;
t108 = height_Body.*t2.*t7.*t20.*t45.*t104.*6.0;
t117 = -t116;
t112 = t49+t55+t90+t92+t94+t98;
t113 = l_Tau_Alpha_Shoulder+r_Tau_Alpha_Shoulder+t56+t81+t82+t89+t91+t95+t97;
t119 = t102+t108+t110+t117;
t114 = t20.*t44.*t113.*1.2e+1;
t118 = t107+t111+t114;
ddl_Shoulder = [(t3.*t4.*t10.*width_Body)./2.0+(t3.*t4.*t12.*width_Body)./2.0+(t3.*t7.*t119.*width_Body)./2.0+t20.*t45.*t73.*(t36+t13.*t17.*3.0)+t4.*t6.*t71.*width_Body.*(t47-t49+t88-t90-t94+t96).*6.0-dbeta_Body.*dgamma_Body.*t6.*t7.*width_Body-height_Body.*t4.*t20.*t45.*t115.*6.0+t2.*t4.*t7.*t13.*t20.*t45.*t104.*3.0+t4.*t5.*t7.*t13.*t20.*t45.*t105.*3.0,(width_Body.*(t9.*t21+t10.*t21+t9.*t28+t10.*t28+t12.*t28-t26.*t118+t30.*t118-t27.*t71.*(t47-t49+t88-t90-t94+t96).*1.2e+1+t29.*t71.*(t47-t49+t88-t90-t94+t96).*1.2e+1-dalpha_Body.*dbeta_Body.*t25.*2.0-dalpha_Body.*dbeta_Body.*t31.*2.0-t4.*t25.*t119+dalpha_Body.*dgamma_Body.*t4.*t27.*2.0+dbeta_Body.*dgamma_Body.*t4.*t26.*2.0))./2.0-t20.*t100.*t105.*t106.*3.0+t20.*t100.*t104.*(t99+t14.*t16.*t19.*3.0+t14.*t17.*t18.*3.0+t16.*t19.*t32.*3.0+t17.*t18.*t35.*3.0)+height_Body.*t4.*t5.*t20.*t44.*t113.*6.0-height_Body.*t2.*t7.*t20.*t45.*t115.*6.0+t2.*t4.*t7.*t13.*t20.*t45.*t73.*3.0,width_Body.*(t9.*t26+t10.*t26+t9.*t40+t10.*t40+t12.*t40+t21.*t118+t28.*t118-t25.*t71.*(t47-t49+t88-t90-t94+t96).*1.2e+1-t31.*t71.*(t47-t49+t88-t90-t94+t96).*1.2e+1+dalpha_Body.*dbeta_Body.*t27.*2.0-dalpha_Body.*dbeta_Body.*t29.*2.0+t4.*t27.*t119+dalpha_Body.*dgamma_Body.*t4.*t25.*2.0-dbeta_Body.*dgamma_Body.*t4.*t21.*2.0).*(-1.0./2.0)-t20.*t100.*t104.*t106.*3.0+t20.*t100.*t105.*(t99+t14.*t16.*t17.*3.0+t14.*t18.*t19.*3.0+t16.*t17.*t35.*3.0+t18.*t19.*t32.*3.0)-height_Body.*t2.*t4.*t20.*t44.*t113.*6.0-height_Body.*t5.*t7.*t20.*t45.*t115.*6.0+t4.*t5.*t7.*t13.*t20.*t45.*t73.*3.0];
