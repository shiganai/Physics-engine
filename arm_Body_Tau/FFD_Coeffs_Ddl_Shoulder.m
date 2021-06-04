function [A11,A12,A13,A14,A15,A16,A17,A21,A22,A23,A24,A25,A26,A27,A31,A32,A33,A34,A35,A36,A37] = FFD_Coeffs_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_Alpha_Hand,l_Tau_Beta_Shoulder,l_Tau_Alpha_Shoulder,m_Body,r_Alpha_Hand,r_Tau_Beta_Shoulder,r_Tau_Alpha_Shoulder,width_Body)
%FFD_COEFFS_DDL_SHOULDER
%    [A11,A12,A13,A14,A15,A16,A17,A21,A22,A23,A24,A25,A26,A27,A31,A32,A33,A34,A35,A36,A37] = FFD_COEFFS_DDL_SHOULDER(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DBETA_BODY,DEPTH_BODY,DGAMMA_BODY,G,GAMMA_BODY,HEIGHT_BODY,L_ALPHA_HAND,L_TAU_BETA_SHOULDER,L_TAU_ALPHA_SHOULDER,M_BODY,R_ALPHA_HAND,R_TAU_BETA_SHOULDER,R_TAU_ALPHA_SHOULDER,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 18:12:57

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = cos(l_Alpha_Hand);
t6 = sin(alpha_Body);
t7 = cos(r_Alpha_Hand);
t8 = sin(beta_Body);
t9 = sin(gamma_Body);
t10 = sin(l_Alpha_Hand);
t11 = sin(r_Alpha_Hand);
t12 = g.*m_Body;
t13 = l_Tau_Alpha_Shoulder+r_Tau_Alpha_Shoulder;
t14 = dalpha_Body.^2;
t15 = dbeta_Body.^2;
t16 = depth_Body.^2;
t17 = dgamma_Body.^2;
t18 = height_Body.^2;
t19 = height_Body.^3;
t21 = width_Body.^2;
t31 = 1.0./m_Body;
t20 = t18.^2;
t22 = t2.^2;
t23 = t4.^2;
t24 = t6.^2;
t25 = t8.^2;
t26 = t9.^2;
t27 = l_Tau_Beta_Shoulder.*t5;
t28 = r_Tau_Beta_Shoulder.*t7;
t29 = l_Tau_Beta_Shoulder.*t10;
t30 = r_Tau_Beta_Shoulder.*t11;
t32 = -t12;
t33 = m_Body.*t16;
t34 = m_Body.*t21;
t35 = t9.*t13;
t36 = t16.*t18;
t37 = t16+t18;
t38 = t16.*t21;
t39 = t18.*t21;
t40 = t18+t21;
t41 = dalpha_Body.*dgamma_Body.*height_Body.*t6.*t9;
t44 = (t2.*t3.*width_Body)./2.0;
t45 = (t2.*t8.*width_Body)./2.0;
t46 = (t3.*t6.*width_Body)./2.0;
t47 = (t6.*t8.*width_Body)./2.0;
t48 = dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*t2.*t9;
t51 = height_Body.*t2.*t4.*t14.*2.0;
t52 = height_Body.*t2.*t4.*t17.*2.0;
t57 = (height_Body.*t2.*t4.*t12)./2.0;
t59 = (height_Body.*t6.*t9.*t12)./2.0;
t64 = (height_Body.*t2.*t4.*t14)./2.0;
t65 = (height_Body.*t2.*t4.*t17)./2.0;
t67 = t2.*t8.*t9.*width_Body.*(-1.0./2.0);
t68 = t3.*t6.*t9.*width_Body.*(-1.0./2.0);
t77 = (height_Body.*m_Body.*t4.*t6.*t14)./2.0;
t78 = (height_Body.*m_Body.*t4.*t6.*t17)./2.0;
t85 = (m_Body.*t4.*t9.*t17.*t18)./4.0;
t42 = t41.*4.0;
t43 = t18.*t23.*3.0;
t49 = t27+t28;
t50 = t29+t30;
t53 = -t41;
t55 = 1.0./t37;
t56 = 1.0./t40;
t58 = t9.*t44;
t60 = t9.*t45;
t61 = t9.*t46;
t62 = t9.*t47;
t63 = t33+t34;
t66 = -t57;
t70 = t2.*t6.*t20.*t23;
t72 = t2.*t6.*t20.*t26;
t73 = t20.*t22.*t23.*3.0;
t74 = t20.*t22.*t26.*3.0;
t75 = t20.*t23.*t24.*3.0;
t76 = t20.*t24.*t26.*3.0;
t83 = t2.*t6.*t26.*t36;
t84 = t2.*t6.*t23.*t39;
t86 = t22.*t26.*t36.*3.0;
t87 = t22.*t23.*t39.*3.0;
t89 = t24.*t26.*t36.*3.0;
t90 = t23.*t24.*t39.*3.0;
t93 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t9.*t18.*t22)./2.0;
t94 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t9.*t18.*t24)./2.0;
t95 = (m_Body.*t4.*t9.*t14.*t18.*t22)./4.0;
t96 = t22.*t85;
t97 = (m_Body.*t4.*t9.*t14.*t18.*t24)./4.0;
t98 = t24.*t85;
t101 = m_Body.*t4.*t9.*t17.*t18.*t22.*(-1.0./4.0);
t103 = m_Body.*t4.*t9.*t17.*t18.*t24.*(-1.0./4.0);
t111 = t45+t68;
t112 = t46+t67;
t123 = t20+t36+t38+t39;
t146 = t32+t48+t77+t78;
t54 = -t42;
t69 = t2.*t49;
t71 = t6.*t50;
t79 = t4.*t6.*t49;
t80 = t2.*t4.*t50;
t81 = 1.0./t63;
t82 = -t72;
t88 = t40+t43;
t92 = -t83;
t99 = height_Body.*t4.*t31.*t56.*6.0;
t100 = -t95;
t102 = -t97;
t104 = t44+t62;
t105 = t47+t58;
t106 = height_Body.*t2.*t4.*t31.*t55.*6.0;
t107 = height_Body.*t4.*t6.*t31.*t55.*6.0;
t108 = height_Body.*t2.*t9.*t31.*t56.*6.0;
t109 = height_Body.*t6.*t9.*t31.*t56.*6.0;
t110 = t3.*t9.*t31.*t56.*width_Body.*6.0;
t116 = t4.*t9.*t17.*t18.*t56.*3.0;
t117 = t2.*t3.*t4.*t31.*t56.*width_Body.*6.0;
t118 = t3.*t4.*t6.*t31.*t56.*width_Body.*6.0;
t122 = height_Body.*t3.*t4.*t9.*t31.*t56.*width_Body.*3.0;
t124 = t2.*t4.*t9.*t18.*t31.*t56.*3.0;
t125 = t4.*t6.*t9.*t18.*t31.*t56.*3.0;
t126 = height_Body.*t2.*t3.*t23.*t31.*t56.*width_Body.*3.0;
t127 = height_Body.*t2.*t3.*t26.*t31.*t56.*width_Body.*3.0;
t128 = height_Body.*t3.*t6.*t23.*t31.*t56.*width_Body.*3.0;
t129 = height_Body.*t3.*t6.*t26.*t31.*t56.*width_Body.*3.0;
t134 = 1.0./t123;
t137 = height_Body.*t2.*t3.*t4.*t6.*t9.*t31.*t56.*width_Body.*-3.0;
t142 = t53+t64+t65;
t145 = t31.*t55.*t111.*1.2e+1;
t169 = t13+t66+t93+t94;
t177 = height_Body.*t2.*t4.*t31.*t55.*t146.*-6.0;
t179 = t73+t76+t87+t89+t123;
t180 = t74+t75+t86+t90+t123;
t91 = -t79;
t113 = -t106;
t114 = -t107;
t115 = -t110;
t119 = -t116;
t120 = -t117;
t121 = -t118;
t130 = t2.*t3.*t4.*t8.*t21.*t81.*3.0;
t131 = t3.*t4.*t6.*t8.*t21.*t81.*3.0;
t132 = t21.*t23.*t25.*t81.*3.0;
t133 = t24.*t122;
t135 = t2.*t6.*t122;
t136 = t22.*t122;
t138 = t2.*t4.*t9.*t21.*t25.*t81.*3.0;
t139 = t4.*t6.*t9.*t21.*t25.*t81.*3.0;
t140 = t51+t52+t54;
t141 = t31.*t56.*t88;
t144 = t31.*t55.*t105.*1.2e+1;
t147 = t99+t110;
t148 = t2.*t3.*t81.*t104.*width_Body.*6.0;
t149 = t3.*t6.*t81.*t104.*width_Body.*6.0;
t150 = t4.*t8.*t81.*t104.*width_Body.*6.0;
t151 = t105.*t106;
t152 = t105.*t107;
t154 = t2.*t3.*t81.*t112.*width_Body.*6.0;
t155 = t3.*t6.*t81.*t112.*width_Body.*6.0;
t156 = t4.*t8.*t81.*t112.*width_Body.*6.0;
t157 = t106.*t111;
t158 = t107.*t111;
t159 = t2.*t8.*t9.*t81.*t104.*width_Body.*6.0;
t160 = t6.*t8.*t9.*t81.*t104.*width_Body.*6.0;
t161 = t108+t117;
t162 = t2.*t8.*t9.*t81.*t112.*width_Body.*6.0;
t163 = t109+t118;
t164 = t6.*t8.*t9.*t81.*t112.*width_Body.*6.0;
t170 = t70+t82+t84+t92;
t173 = t107+t145;
t174 = t109.*t146;
t176 = t106.*t146;
t178 = t31.*t55.*t169.*1.2e+1;
t183 = t31.*t134.*t179;
t184 = t31.*t134.*t180;
t185 = t59+t69+t71+t85+t100+t101+t102+t103;
t143 = t35+t80+t91;
t153 = t99+t115;
A11 = -t122-t132+t141+(t3.*t9.*t153.*width_Body)./2.0;
if nargout > 1
    A12 = t124+t126+t156+(t3.*t9.*t161.*width_Body)./2.0;
end
if nargout > 2
    A13 = t125+t128-t150+(t3.*t9.*t163.*width_Body)./2.0;
end
if nargout > 3
    A14 = t122+t132+t141+(t3.*t9.*t147.*width_Body)./2.0;
end
if nargout > 4
    t165 = t108+t120;
    A15 = t124-t126-t156+(t3.*t9.*t165.*width_Body)./2.0;
end
if nargout > 5
    t166 = t109+t121;
    A16 = t125-t128+t150+(t3.*t9.*t166.*width_Body)./2.0;
end
if nargout > 6
    t167 = height_Body.*t4.*t6.*t55.*t140.*(3.0./2.0);
    t168 = height_Body.*t2.*t9.*t56.*t140.*(3.0./2.0);
    t171 = t106+t144;
    t172 = t113+t144;
    t175 = t114+t145;
    t181 = t31.*t134.*t170.*3.0;
    t186 = t31.*t56.*t185.*1.2e+1;
    t182 = -t181;
    t187 = t167+t177+t178;
    t188 = t119+t168+t174+t186;
    A17 = t125.*t146+t99.*t185+(t3.*t4.*t15.*width_Body)./2.0+(t3.*t4.*t17.*width_Body)./2.0+(t3.*t9.*t188.*width_Body)./2.0-dbeta_Body.*dgamma_Body.*t8.*t9.*width_Body-(height_Body.*t9.*t17.*t56.*t88)./2.0-t4.*t8.*t81.*t143.*width_Body.*6.0+t2.*t4.*t9.*t18.*t56.*t140.*(3.0./4.0);
end
if nargout > 7
    A21 = t124-t127+t131-t138-(t2.*t3.*t4.*t153.*width_Body)./2.0;
end
if nargout > 8
    A22 = t136-t155+t158+t162+t184+t61.*t173-(t2.*t8.*t173.*width_Body)./2.0-(t2.*t3.*t4.*t161.*width_Body)./2.0;
end
if nargout > 9
    A23 = t135+t149+t152-t159+t182-t61.*(t106-t144)+(t2.*t8.*width_Body.*(t106-t144))./2.0-(t2.*t3.*t4.*t163.*width_Body)./2.0;
end
if nargout > 10
    A24 = t124+t127-t131+t138-(t2.*t3.*t4.*t147.*width_Body)./2.0;
end
if nargout > 11
    A25 = t155-t162+t184-t45.*(t107-t145)-t68.*(t107-t145)-(t2.*t3.*t4.*t165.*width_Body)./2.0-height_Body.*t4.*t6.*t31.*t55.*t111.*6.0-height_Body.*t3.*t4.*t9.*t22.*t31.*t56.*width_Body.*3.0;
end
if nargout > 12
    A26 = t137-t149+t159+t182+t45.*t171+t68.*t171-(t2.*t3.*t4.*t166.*width_Body)./2.0-height_Body.*t4.*t6.*t31.*t55.*t105.*6.0;
end
if nargout > 13
    A27 = t14.*t47+t15.*t47+t14.*t58+t15.*t58+t17.*t58+t61.*t187+t107.*t169+t108.*t185+t134.*t142.*t180-t31.*t134.*t146.*t170.*3.0-(t2.*t8.*t187.*width_Body)./2.0-dalpha_Body.*dbeta_Body.*t2.*t3.*width_Body-(t2.*t3.*t4.*t188.*width_Body)./2.0+t3.*t6.*t81.*t143.*width_Body.*6.0-t2.*t4.*t17.*t19.*t26.*t56.*(3.0./2.0)-t2.*t8.*t9.*t81.*t143.*width_Body.*6.0-dalpha_Body.*dbeta_Body.*t6.*t8.*t9.*width_Body+dalpha_Body.*dgamma_Body.*t3.*t4.*t6.*width_Body+dbeta_Body.*dgamma_Body.*t2.*t4.*t8.*width_Body;
end
if nargout > 14
    A31 = t125-t129-t130-t139-(t3.*t4.*t6.*t153.*width_Body)./2.0;
end
if nargout > 15
    A32 = t135+t154+t164+t182-(t6.*t8.*t173.*width_Body)./2.0-(t3.*t4.*t6.*t161.*width_Body)./2.0-(t2.*t3.*t9.*t173.*width_Body)./2.0-height_Body.*t2.*t4.*t31.*t55.*t111.*6.0;
end
if nargout > 16
    A33 = t133-t148-t160+t183+(t6.*t8.*width_Body.*(t106-t144))./2.0-(t3.*t4.*t6.*t163.*width_Body)./2.0+(t2.*t3.*t9.*width_Body.*(t106-t144))./2.0-height_Body.*t2.*t4.*t31.*t55.*t105.*6.0;
end
if nargout > 17
    A34 = t125+t129+t130+t139-(t3.*t4.*t6.*t147.*width_Body)./2.0;
end
if nargout > 18
    A35 = t137-t154+t157-t164+t182-t47.*(t107-t145)-t58.*(t107-t145)-(t3.*t4.*t6.*t165.*width_Body)./2.0;
end
if nargout > 19
    A36 = t148+t151+t160+t183+t47.*t171+t58.*t171-(t3.*t4.*t6.*t166.*width_Body)./2.0-height_Body.*t3.*t4.*t9.*t24.*t31.*t56.*width_Body.*3.0;
end
if nargout > 20
    A37 = t14.*t61+t15.*t61+t17.*t61+t109.*t185+t146.*t183-t134.*t142.*t170.*3.0-(t2.*t8.*t14.*width_Body)./2.0-(t2.*t8.*t15.*width_Body)./2.0-(t6.*t8.*t187.*width_Body)./2.0-dalpha_Body.*dbeta_Body.*t3.*t6.*width_Body-(t2.*t3.*t9.*t187.*width_Body)./2.0-(t3.*t4.*t6.*t188.*width_Body)./2.0-t2.*t3.*t81.*t143.*width_Body.*6.0-height_Body.*t2.*t4.*t31.*t55.*t169.*6.0-t4.*t6.*t17.*t19.*t26.*t56.*(3.0./2.0)-t6.*t8.*t9.*t81.*t143.*width_Body.*6.0+dalpha_Body.*dbeta_Body.*t2.*t8.*t9.*width_Body-dalpha_Body.*dgamma_Body.*t2.*t3.*t4.*width_Body+dbeta_Body.*dgamma_Body.*t4.*t6.*t8.*width_Body;
end
