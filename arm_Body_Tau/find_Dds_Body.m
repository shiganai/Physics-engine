function [ddalpha_Body,ddx_Head,ddy_Head,ddz_Head,r_Tau_Beta_Shoulder,l_Tau_Beta_Shoulder] = find_Dds_Body(alpha_Body,beta_Body,dalpha_Body,ddbeta_Body,ddgamma_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_Beta_Hand,l_F_X,l_F_Y,l_F_Z,l_Tau_Alpha_Shoulder,m_Body,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z,r_Tau_Alpha_Shoulder,width_Body)
%FIND_DDS_BODY
%    [DDALPHA_BODY,DDX_HEAD,DDY_HEAD,DDZ_HEAD,R_TAU_BETA_SHOULDER,L_TAU_BETA_SHOULDER] = FIND_DDS_BODY(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DDBETA_BODY,DDGAMMA_BODY,DEPTH_BODY,DGAMMA_BODY,G,GAMMA_BODY,HEIGHT_BODY,L_BETA_HAND,L_F_X,L_F_Y,L_F_Z,L_TAU_ALPHA_SHOULDER,M_BODY,R_BETA_HAND,R_F_X,R_F_Y,R_F_Z,R_TAU_ALPHA_SHOULDER,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 15:17:58

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = conj(l_Beta_Hand);
t6 = sin(alpha_Body);
t7 = conj(r_Beta_Hand);
t8 = sin(beta_Body);
t9 = sin(gamma_Body);
t11 = g.*m_Body;
t15 = dalpha_Body.^2;
t16 = depth_Body.^2;
t17 = dgamma_Body.^2;
t18 = height_Body.^2;
t19 = width_Body.^2;
t10 = cos(t5);
t12 = cos(t7);
t13 = sin(t5);
t14 = sin(t7);
t20 = t2.^2;
t21 = t4.^2;
t22 = t6.^2;
t23 = t9.^2;
t27 = -t11;
t29 = ddgamma_Body.*height_Body.*t4.*2.0;
t32 = height_Body.*t9.*t17.*2.0;
t33 = ddgamma_Body.*height_Body.*t2.*t9.*2.0;
t34 = dalpha_Body.*dgamma_Body.*height_Body.*t6.*t9.*4.0;
t35 = (t2.*t3.*width_Body)./2.0;
t36 = (t2.*t8.*width_Body)./2.0;
t37 = (t3.*t6.*width_Body)./2.0;
t38 = (t6.*t8.*width_Body)./2.0;
t39 = dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*t2.*t9;
t40 = (ddbeta_Body.*m_Body.*t19)./1.2e+1;
t41 = (ddgamma_Body.*m_Body.*t19)./1.2e+1;
t42 = height_Body.*t2.*t4.*t15.*2.0;
t43 = height_Body.*t2.*t4.*t17.*2.0;
t45 = (l_F_X.*t3.*t9.*width_Body)./2.0;
t46 = (l_F_X.*t4.*t8.*width_Body)./2.0;
t47 = (r_F_X.*t3.*t9.*width_Body)./2.0;
t48 = (r_F_X.*t4.*t8.*width_Body)./2.0;
t49 = (ddbeta_Body.*m_Body.*t16)./1.2e+1;
t50 = (ddgamma_Body.*m_Body.*t18)./1.2e+1;
t51 = (height_Body.*t2.*t4.*t11)./2.0;
t53 = (ddgamma_Body.*height_Body.*m_Body.*t6.*t9)./2.0;
t54 = (height_Body.*t6.*t9.*t11)./2.0;
t63 = t2.*t8.*t9.*width_Body.*(-1.0./2.0);
t64 = t3.*t6.*t9.*width_Body.*(-1.0./2.0);
t71 = (height_Body.*m_Body.*t4.*t6.*t15)./2.0;
t72 = (height_Body.*m_Body.*t4.*t6.*t17)./2.0;
t73 = r_F_Y.*t2.*t3.*t4.*width_Body.*(-1.0./2.0);
t74 = r_F_Z.*t3.*t4.*t6.*width_Body.*(-1.0./2.0);
t76 = (m_Body.*t4.*t9.*t17.*t18)./4.0;
t24 = r_Tau_Alpha_Shoulder.*t12;
t25 = l_Tau_Alpha_Shoulder.*t13;
t26 = r_Tau_Alpha_Shoulder.*t14;
t28 = l_Tau_Alpha_Shoulder.*t10;
t31 = -t29;
t44 = -t34;
t52 = t9.*t35;
t55 = t9.*t36;
t56 = t9.*t37;
t57 = t9.*t38;
t59 = -t45;
t60 = -t46;
t61 = -t51;
t62 = -t54;
t65 = l_F_Y.*t4.*t35;
t66 = r_F_Y.*t4.*t35;
t67 = l_F_Z.*t4.*t37;
t68 = r_F_Z.*t4.*t37;
t69 = (ddgamma_Body.*m_Body.*t18.*t21)./4.0;
t79 = -t76;
t80 = (ddgamma_Body.*m_Body.*t18.*t20.*t23)./4.0;
t81 = (ddgamma_Body.*m_Body.*t18.*t22.*t23)./4.0;
t82 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t9.*t18.*t20)./2.0;
t83 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t9.*t18.*t22)./2.0;
t86 = (m_Body.*t4.*t9.*t15.*t18.*t20)./4.0;
t87 = t20.*t76;
t88 = (m_Body.*t4.*t9.*t15.*t18.*t22)./4.0;
t89 = t22.*t76;
t92 = m_Body.*(t29-t32).*(-1.0./4.0);
t93 = t36+t64;
t94 = t37+t63;
t95 = (m_Body.*(t29-t32))./4.0;
t112 = l_F_Z+r_F_Z+t27+t39+t53+t71+t72;
t30 = -t26;
t58 = t24+t28;
t77 = t31+t32;
t90 = t35+t57;
t91 = t38+t52;
t100 = l_F_Y.*t93;
t101 = l_F_Y.*t94;
t102 = r_F_Y.*t93;
t103 = r_F_Y.*t94;
t104 = l_F_X+r_F_X+t95;
t109 = t33+t42+t43+t44;
t70 = t25+t30;
t75 = t9.*t58;
t96 = l_F_Z.*t90;
t97 = l_F_Z.*t91;
t98 = r_F_Z.*t90;
t99 = r_F_Z.*t91;
t107 = -t100;
t108 = -t103;
t110 = (m_Body.*t109)./4.0;
t78 = t6.*t70;
t84 = t2.*t4.*t70;
t105 = -t96;
t106 = -t97;
t111 = l_F_Y+r_F_Y+t110;
t85 = -t78;
t113 = t58+t61+t82+t83+t99+t102+t106+t107;
t114 = t40+t48+t49+t60+t75+t84+t98+t101+t105+t108;
t115 = sign(t114);
t117 = t41+t47+t50+t59+t62+t65+t67+t69+t73+t74+t79+t80+t81+t85+t86+t87+t88+t89;
t116 = -t115;
t118 = sign(t117);
t119 = -t118;
t120 = t104.*Inf+t111.*Inf+t112.*Inf+t113.*Inf+t116.*Inf+t119.*Inf;
ddalpha_Body = t120;
if nargout > 1
    ddx_Head = t120;
end
if nargout > 2
    ddy_Head = t120;
end
if nargout > 3
    ddz_Head = t120;
end
if nargout > 4
    r_Tau_Beta_Shoulder = t120;
end
if nargout > 5
    l_Tau_Beta_Shoulder = t120;
end
