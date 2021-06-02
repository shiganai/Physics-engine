function ddr_Shoulder = find_Ddr_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,offset_Body,r_F_X,r_F_Y,r_F_Z,width_Body)
%FIND_DDR_SHOULDER
%    DDR_SHOULDER = FIND_DDR_SHOULDER(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DBETA_BODY,DEPTH_BODY,DGAMMA_BODY,DOFFSET_BODY,G,GAMMA_BODY,HEIGHT_BODY,L_F_X,L_F_Y,L_F_Z,M_BODY,OFFSET_BODY,R_F_X,R_F_Y,R_F_Z,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    02-Jun-2021 19:37:01

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = sin(alpha_Body);
t6 = sin(beta_Body);
t7 = sin(gamma_Body);
t8 = dalpha_Body.^2;
t9 = dbeta_Body.^2;
t10 = depth_Body.^2;
t11 = dgamma_Body.^2;
t12 = height_Body.^2;
t13 = offset_Body.^2;
t14 = width_Body.^2;
t25 = 1.0./m_Body;
t15 = t2.^2;
t16 = t4.^2;
t17 = t5.^2;
t18 = t7.^2;
t19 = l_F_Y.*t2;
t20 = offset_Body.*t2;
t21 = r_F_Y.*t2;
t22 = l_F_Z.*t5;
t23 = offset_Body.*t5;
t24 = r_F_Z.*t5;
t26 = g.*m_Body.*t5;
t27 = t5.*t6;
t28 = m_Body.*t10;
t29 = m_Body.*t12;
t30 = m_Body.*t14;
t31 = t2.*t3;
t32 = t2.*t6;
t33 = t3.*t5;
t51 = (l_F_X.*t3.*t7.*width_Body)./2.0;
t52 = (l_F_X.*t4.*t6.*width_Body)./2.0;
t53 = (r_F_X.*t3.*t7.*width_Body)./2.0;
t54 = (r_F_X.*t4.*t6.*width_Body)./2.0;
t55 = (g.*height_Body.*m_Body.*t2.*t4)./2.0;
t34 = m_Body.*t17;
t35 = g.*m_Body.*t20;
t36 = -t20;
t37 = -t23;
t38 = m_Body.*t15;
t39 = -t26;
t40 = t7.*t31;
t41 = t7.*t32;
t42 = t7.*t33;
t43 = t7.*t27;
t44 = t12.*t16.*3.0;
t56 = (height_Body.*t7.*t26)./2.0;
t57 = t28+t30;
t60 = -t51;
t61 = -t52;
t68 = (t3.*t4.*t19.*width_Body)./2.0;
t69 = (t3.*t4.*t21.*width_Body)./2.0;
t71 = (t3.*t4.*t22.*width_Body)./2.0;
t72 = (t3.*t4.*t24.*width_Body)./2.0;
t76 = t12.*t15.*t18.*3.0;
t77 = t12.*t17.*t18.*3.0;
t82 = t15.*t16.*t29.*3.0;
t84 = t16.*t17.*t29.*3.0;
t91 = (t4.*t7.*t11.*t29)./4.0;
t98 = (dalpha_Body.*dgamma_Body.*t4.*t7.*t15.*t29)./2.0;
t99 = (dalpha_Body.*dgamma_Body.*t4.*t7.*t17.*t29)./2.0;
t100 = (t4.*t7.*t8.*t15.*t29)./4.0;
t103 = (t4.*t7.*t8.*t17.*t29)./4.0;
t45 = -t41;
t46 = -t42;
t47 = offset_Body.*t8.*t38;
t48 = dalpha_Body.*doffset_Body.*offset_Body.*t38.*2.0;
t49 = offset_Body.*t8.*t34;
t50 = dalpha_Body.*doffset_Body.*offset_Body.*t34.*2.0;
t58 = t13.*t38.*1.2e+1;
t59 = t13.*t34.*1.2e+1;
t62 = dalpha_Body.*doffset_Body.*height_Body.*t4.*t38;
t63 = dalpha_Body.*doffset_Body.*height_Body.*t4.*t34;
t64 = -t56;
t65 = dalpha_Body.*dgamma_Body.*height_Body.*offset_Body.*t7.*t38;
t66 = dalpha_Body.*dgamma_Body.*height_Body.*offset_Body.*t7.*t34;
t67 = height_Body.*offset_Body.*t4.*t38.*1.2e+1;
t70 = height_Body.*offset_Body.*t4.*t34.*1.2e+1;
t73 = t27+t40;
t74 = t31+t43;
t75 = t34+t38;
t78 = -t69;
t79 = -t72;
t83 = 1.0./t57;
t87 = (height_Body.*t4.*t8.*t38)./2.0;
t88 = (height_Body.*t4.*t11.*t38)./2.0;
t89 = (height_Body.*t4.*t8.*t34)./2.0;
t90 = (height_Body.*t4.*t11.*t34)./2.0;
t92 = t12+t14+t44;
t96 = -t91;
t101 = t15.*t91;
t102 = -t98;
t104 = t17.*t91;
t105 = -t99;
t80 = -t65;
t81 = -t66;
t85 = t32+t46;
t86 = t33+t45;
t93 = (height_Body.*t7.*t47)./2.0;
t94 = (height_Body.*t7.*t49)./2.0;
t95 = 1.0./t75;
t97 = (t73.*width_Body)./2.0;
t106 = (l_F_Z.*t74.*width_Body)./2.0;
t107 = (r_F_Z.*t74.*width_Body)./2.0;
t108 = 1.0./t92;
t124 = t76+t77+t92;
t125 = t28+t29+t58+t59+t67+t70+t82+t84;
t127 = t19+t21+t22+t24+t39+t47+t49+t87+t88+t89+t90;
t109 = (t85.*width_Body)./2.0;
t110 = -t106;
t111 = (l_F_Y.*t86.*width_Body)./2.0;
t112 = (r_F_Y.*t86.*width_Body)./2.0;
t114 = t20+t97;
t117 = t36+t97;
t119 = -l_F_Z.*(t20-t97);
t126 = 1.0./t125;
t129 = height_Body.*t7.*t25.*t108.*t127.*6.0;
t131 = t95.*t108.*t124.*t127;
t132 = t53+t60+t64+t68+t71+t78+t79+t93+t94+t96+t100+t101+t103+t104;
t113 = -t112;
t115 = t23+t109;
t116 = r_F_Z.*t114;
t118 = t37+t109;
t122 = -r_F_Y.*(t23-t109);
t123 = r_F_Y.*(t23-t109);
t130 = -t129;
t133 = t25.*t108.*t132.*1.2e+1;
t134 = height_Body.*t7.*t25.*t108.*t132.*6.0;
t120 = -t116;
t121 = l_F_Y.*t115;
t128 = t54+t61+t107+t110+t111+t113;
t135 = -t134;
t137 = t130+t133;
t136 = t35+t48+t50+t55+t62+t63+t80+t81+t102+t105+t119+t120+t121+t123;
t138 = t131+t135;
ddr_Shoulder = [t3.*t4.*t9.*width_Body.*(-1.0./2.0)-(t3.*t4.*t11.*width_Body)./2.0-(t3.*t7.*width_Body.*(t129-t133))./2.0-t4.*t6.*t83.*width_Body.*(t52-t54+t106-t107-t111+t112).*6.0+dbeta_Body.*dgamma_Body.*t6.*t7.*width_Body,t8.*t36+t2.*t138-(width_Body.*(t8.*t27+t9.*t27+t8.*t40+t9.*t40+t11.*t40-t33.*t83.*(t52-t54+t106-t107-t111+t112).*1.2e+1+t41.*t83.*(t52-t54+t106-t107-t111+t112).*1.2e+1-dalpha_Body.*dbeta_Body.*t31.*2.0-dalpha_Body.*dbeta_Body.*t43.*2.0+t32.*t126.*t136.*1.2e+1-t42.*t126.*t136.*1.2e+1-t4.*t31.*(t129-t133)+dalpha_Body.*dgamma_Body.*t4.*t33.*2.0+dbeta_Body.*dgamma_Body.*t4.*t32.*2.0))./2.0-dalpha_Body.*doffset_Body.*t5.*2.0+t23.*t126.*t136.*1.2e+1,t8.*t37+t5.*t138-(width_Body.*(-t8.*t32-t9.*t32+t8.*t42+t9.*t42+t11.*t42+t31.*t83.*(t52-t54+t106-t107-t111+t112).*1.2e+1+t43.*t83.*(t52-t54+t106-t107-t111+t112).*1.2e+1-dalpha_Body.*dbeta_Body.*t33.*2.0+dalpha_Body.*dbeta_Body.*t41.*2.0+t27.*t126.*t136.*1.2e+1+t40.*t126.*t136.*1.2e+1-t4.*t33.*(t129-t133)-dalpha_Body.*dgamma_Body.*t4.*t31.*2.0+dbeta_Body.*dgamma_Body.*t4.*t27.*2.0))./2.0+dalpha_Body.*doffset_Body.*t2.*2.0-t20.*t126.*t136.*1.2e+1];
