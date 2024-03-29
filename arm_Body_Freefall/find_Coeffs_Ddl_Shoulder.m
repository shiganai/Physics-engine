function coeffs_Ddl_Shoulder = find_Coeffs_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,m_Body,width_Body)
%FIND_COEFFS_DDL_SHOULDER
%    COEFFS_DDL_SHOULDER = FIND_COEFFS_DDL_SHOULDER(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DBETA_BODY,DEPTH_BODY,DGAMMA_BODY,G,GAMMA_BODY,HEIGHT_BODY,M_BODY,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    02-Jun-2021 20:43:59

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
t14 = height_Body.^3;
t16 = width_Body.^2;
t22 = 1.0./m_Body;
t15 = t13.^2;
t17 = t2.^2;
t18 = t4.^2;
t19 = t5.^2;
t20 = t6.^2;
t21 = t7.^2;
t23 = t5.*t6;
t24 = -t8;
t25 = m_Body.*t11;
t26 = m_Body.*t16;
t27 = t2.*t3;
t28 = t2.*t6;
t29 = t3.*t5;
t34 = t11.*t13;
t35 = t11+t13;
t36 = t11.*t16;
t37 = t13.*t16;
t38 = t13+t16;
t39 = dalpha_Body.*dgamma_Body.*height_Body.*t5.*t7;
t42 = dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*t2.*t7;
t45 = height_Body.*t2.*t4.*t9.*2.0;
t46 = height_Body.*t2.*t4.*t12.*2.0;
t51 = (height_Body.*t2.*t4.*t8)./2.0;
t52 = (height_Body.*t5.*t7.*t8)./2.0;
t54 = (height_Body.*t2.*t4.*t9)./2.0;
t55 = (height_Body.*t2.*t4.*t12)./2.0;
t66 = (height_Body.*m_Body.*t4.*t5.*t9)./2.0;
t67 = (height_Body.*m_Body.*t4.*t5.*t12)./2.0;
t74 = (m_Body.*t4.*t7.*t12.*t13)./4.0;
t30 = t7.*t27;
t31 = t7.*t28;
t32 = t7.*t29;
t33 = t7.*t23;
t40 = t39.*4.0;
t41 = t13.*t18.*3.0;
t47 = -t39;
t49 = 1.0./t35;
t50 = 1.0./t38;
t53 = t25+t26;
t56 = -t51;
t57 = -t52;
t58 = t2.*t5.*t15.*t18;
t59 = t2.*t5.*t15.*t21;
t62 = t15.*t17.*t18.*3.0;
t63 = t15.*t17.*t21.*3.0;
t64 = t15.*t18.*t19.*3.0;
t65 = t15.*t19.*t21.*3.0;
t72 = t2.*t5.*t21.*t34;
t73 = t2.*t5.*t18.*t37;
t75 = t17.*t21.*t34.*3.0;
t76 = t17.*t18.*t37.*3.0;
t78 = t19.*t21.*t34.*3.0;
t79 = t18.*t19.*t37.*3.0;
t80 = -t74;
t82 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t7.*t13.*t17)./2.0;
t83 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t7.*t13.*t19)./2.0;
t84 = (m_Body.*t4.*t7.*t9.*t13.*t17)./4.0;
t85 = t17.*t74;
t86 = (m_Body.*t4.*t7.*t9.*t13.*t19)./4.0;
t87 = t19.*t74;
t104 = t15+t34+t36+t37;
t139 = t24+t42+t66+t67;
t43 = -t31;
t44 = -t32;
t48 = -t40;
t60 = t23+t30;
t61 = t27+t33;
t68 = 1.0./t53;
t69 = -t59;
t77 = t38+t41;
t81 = -t72;
t88 = height_Body.*t4.*t22.*t50.*6.0;
t89 = height_Body.*t2.*t4.*t22.*t49.*6.0;
t90 = height_Body.*t4.*t5.*t22.*t49.*6.0;
t91 = height_Body.*t2.*t7.*t22.*t50.*6.0;
t92 = height_Body.*t5.*t7.*t22.*t50.*6.0;
t93 = t3.*t7.*t22.*t50.*width_Body.*6.0;
t97 = t4.*t7.*t12.*t13.*t50.*3.0;
t98 = t4.*t22.*t27.*t50.*width_Body.*6.0;
t99 = t4.*t22.*t29.*t50.*width_Body.*6.0;
t102 = height_Body.*t3.*t4.*t7.*t22.*t50.*width_Body.*3.0;
t105 = t2.*t4.*t7.*t13.*t22.*t50.*3.0;
t106 = t4.*t5.*t7.*t13.*t22.*t50.*3.0;
t108 = height_Body.*t18.*t22.*t27.*t50.*width_Body.*3.0;
t109 = height_Body.*t21.*t22.*t27.*t50.*width_Body.*3.0;
t110 = height_Body.*t18.*t22.*t29.*t50.*width_Body.*3.0;
t111 = height_Body.*t21.*t22.*t29.*t50.*width_Body.*3.0;
t116 = 1.0./t104;
t117 = height_Body.*t4.*t5.*t22.*t30.*t50.*width_Body.*3.0;
t128 = t47+t54+t55;
t146 = t56+t82+t83;
t157 = height_Body.*t5.*t7.*t22.*t50.*t139.*-6.0;
t158 = height_Body.*t2.*t4.*t22.*t49.*t139.*-6.0;
t160 = t62+t65+t76+t78+t104;
t161 = t63+t64+t75+t79+t104;
t164 = t57+t80+t84+t85+t86+t87;
t70 = t28+t44;
t71 = t29+t43;
t94 = -t89;
t95 = -t90;
t96 = -t93;
t100 = -t98;
t101 = -t99;
t103 = t3.*t4.*t23.*t68.*width_Body.*6.0;
t107 = t4.*t6.*t27.*t68.*width_Body.*6.0;
t112 = t2.*t4.*t7.*t20.*t68.*width_Body.*6.0;
t113 = t4.*t5.*t7.*t20.*t68.*width_Body.*6.0;
t114 = t16.*t18.*t20.*t68.*3.0;
t115 = t19.*t102;
t118 = t17.*t102;
t119 = -t117;
t120 = t22.*t49.*t60.*width_Body.*6.0;
t121 = t45+t46+t48;
t122 = t22.*t50.*t77;
t124 = t27.*t61.*t68.*width_Body.*6.0;
t125 = t29.*t61.*t68.*width_Body.*6.0;
t126 = height_Body.*t2.*t4.*t22.*t49.*t60.*width_Body.*3.0;
t127 = height_Body.*t4.*t5.*t22.*t49.*t60.*width_Body.*3.0;
t129 = t31.*t61.*t68.*width_Body.*6.0;
t130 = t33.*t61.*t68.*width_Body.*6.0;
t133 = t4.*t6.*t16.*t61.*t68.*3.0;
t140 = t88+t93;
t142 = t91+t98;
t143 = t92+t99;
t154 = t58+t69+t73+t81;
t155 = t92.*t139;
t156 = t89.*t139;
t159 = t22.*t49.*t146.*1.2e+1;
t165 = t22.*t116.*t160;
t166 = t22.*t116.*t161;
t167 = t22.*t50.*t164.*1.2e+1;
t123 = t22.*t49.*t70.*width_Body.*6.0;
t131 = t27.*t68.*t71.*width_Body.*6.0;
t132 = t29.*t68.*t71.*width_Body.*6.0;
t134 = height_Body.*t2.*t4.*t22.*t49.*t70.*width_Body.*3.0;
t135 = height_Body.*t4.*t5.*t22.*t49.*t70.*width_Body.*3.0;
t136 = t31.*t68.*t71.*width_Body.*6.0;
t137 = t33.*t68.*t71.*width_Body.*6.0;
t138 = t4.*t6.*t16.*t68.*t71.*3.0;
t141 = t88+t96;
t144 = t91+t100;
t145 = t92+t101;
t147 = height_Body.*t4.*t5.*t49.*t121.*(3.0./2.0);
t148 = height_Body.*t2.*t7.*t50.*t121.*(3.0./2.0);
t150 = t89+t120;
t151 = t94+t120;
t162 = t22.*t116.*t154.*3.0;
t149 = -t148;
t152 = t90+t123;
t153 = t95+t123;
t163 = -t162;
t168 = t147+t158+t159;
t169 = t97+t149+t157+t167;
coeffs_Ddl_Shoulder = reshape([-t102-t114+t122+(t3.*t7.*t141.*width_Body)./2.0,t105-t109-(width_Body.*(-t103+t112+t4.*t27.*t141))./2.0,t106-t111-(width_Body.*(t107+t113+t4.*t29.*t141))./2.0,t105+t108+t138+(t3.*t7.*t142.*width_Body)./2.0,t118+t135+t166-(width_Body.*(t132-t136+t28.*t152+t44.*t152+t4.*t27.*t142))./2.0,t117-t134+t163-(width_Body.*(-t131-t137+t23.*t152+t30.*t152+t4.*t29.*t142))./2.0,t106+t110-t133+(t3.*t7.*t143.*width_Body)./2.0,t117+t127+t163+(width_Body.*(t125-t129+t28.*(t89-t120)+t44.*(t89-t120)-t4.*t27.*t143))./2.0,t115-t126+t165-(width_Body.*(t124+t130-t23.*(t89-t120)-t30.*(t89-t120)+t4.*t29.*t143))./2.0,t102+t114+t122+(t3.*t7.*t140.*width_Body)./2.0,t105+t109-(width_Body.*(t103-t112+t4.*t27.*t140))./2.0,t106+t111+(width_Body.*(t107+t113-t4.*t29.*t140))./2.0,t105-t108-t138+(t3.*t7.*t144.*width_Body)./2.0,-t135+t166-(width_Body.*(-t132+t136+t28.*(t90-t123)+t44.*(t90-t123)+t4.*t27.*t144))./2.0-height_Body.*t3.*t4.*t7.*t17.*t22.*t50.*width_Body.*3.0,t119+t134+t163-(width_Body.*(t131+t137+t23.*(t90-t123)+t30.*(t90-t123)+t4.*t29.*t144))./2.0,t106-t110+t133+(t3.*t7.*t145.*width_Body)./2.0,t119-t127+t163+(width_Body.*(-t125+t129+t28.*t150+t44.*t150-t4.*t27.*t145))./2.0,t126+t165+(width_Body.*(t124+t130+t23.*t150+t30.*t150-t4.*t29.*t145))./2.0-height_Body.*t3.*t4.*t7.*t19.*t22.*t50.*width_Body.*3.0,t106.*t139+(t3.*t4.*t10.*width_Body)./2.0+(t3.*t4.*t12.*width_Body)./2.0-(t3.*t7.*t169.*width_Body)./2.0-dbeta_Body.*dgamma_Body.*t6.*t7.*width_Body-(height_Body.*t7.*t12.*t50.*t77)./2.0-height_Body.*t4.*t22.*t50.*t164.*6.0+t2.*t4.*t7.*t13.*t50.*t121.*(3.0./4.0),(width_Body.*(t9.*t23+t10.*t23+t9.*t30+t10.*t30+t12.*t30-t28.*t168+t32.*t168-dalpha_Body.*dbeta_Body.*t27.*2.0-dalpha_Body.*dbeta_Body.*t33.*2.0+t4.*t27.*t169+dalpha_Body.*dgamma_Body.*t4.*t29.*2.0+dbeta_Body.*dgamma_Body.*t4.*t28.*2.0))./2.0+t90.*t146+t116.*t128.*t161-t22.*t116.*t139.*t154.*3.0-height_Body.*t2.*t7.*t22.*t50.*t164.*6.0-t2.*t4.*t12.*t14.*t21.*t50.*(3.0./2.0),width_Body.*(t9.*t28+t10.*t28+t9.*t44+t10.*t44+t12.*t44+t23.*t168+t30.*t168+dalpha_Body.*dbeta_Body.*t29.*2.0-dalpha_Body.*dbeta_Body.*t31.*2.0-t4.*t29.*t169+dalpha_Body.*dgamma_Body.*t4.*t27.*2.0-dbeta_Body.*dgamma_Body.*t4.*t23.*2.0).*(-1.0./2.0)+t139.*t165-t116.*t128.*t154.*3.0-height_Body.*t2.*t4.*t22.*t49.*t146.*6.0-height_Body.*t5.*t7.*t22.*t50.*t164.*6.0-t4.*t5.*t12.*t14.*t21.*t50.*(3.0./2.0)],[3,7]);
