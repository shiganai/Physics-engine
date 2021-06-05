function [ddalpha_Body,ddbeta_Body,ddgamma_Body,ddx_Head,ddy_Head,ddz_Head] = FFD_Dds_Body(alpha_Body,beta_Body,dalpha_Body,depth_Body,dgamma_Body,g,gamma_Body,height_Body,l_Alpha_Hand,l_F_X,l_F_Y,l_F_Z,l_Tau_Beta_Shoulder,l_Tau_Alpha_Shoulder,m_Body,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z,r_Tau_Beta_Shoulder,r_Tau_Alpha_Shoulder,width_Body)
%FFD_DDS_BODY
%    [DDALPHA_BODY,DDBETA_BODY,DDGAMMA_BODY,DDX_HEAD,DDY_HEAD,DDZ_HEAD] = FFD_DDS_BODY(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DEPTH_BODY,DGAMMA_BODY,G,GAMMA_BODY,HEIGHT_BODY,L_ALPHA_HAND,L_F_X,L_F_Y,L_F_Z,L_TAU_BETA_SHOULDER,L_TAU_ALPHA_SHOULDER,M_BODY,R_BETA_HAND,R_F_X,R_F_Y,R_F_Z,R_TAU_BETA_SHOULDER,R_TAU_ALPHA_SHOULDER,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    05-Jun-2021 09:07:45

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = cos(l_Alpha_Hand);
t6 = sin(alpha_Body);
t7 = cos(r_Beta_Hand);
t8 = sin(beta_Body);
t9 = sin(gamma_Body);
t10 = sin(l_Alpha_Hand);
t11 = sin(r_Beta_Hand);
t12 = g.*m_Body;
t13 = dalpha_Body.^2;
t14 = depth_Body.^2;
t15 = dgamma_Body.^2;
t16 = height_Body.^2;
t18 = width_Body.^2;
t27 = 1.0./m_Body;
t17 = t16.^2;
t19 = t2.^2;
t20 = t4.^2;
t21 = t6.^2;
t22 = t9.^2;
t23 = l_Tau_Beta_Shoulder.*t5;
t24 = r_Tau_Alpha_Shoulder.*t7;
t25 = l_Tau_Beta_Shoulder.*t10;
t26 = r_Tau_Alpha_Shoulder.*t11;
t28 = -t12;
t31 = t14.*t16;
t32 = t14+t16;
t33 = t14.*t18;
t34 = t16.*t18;
t35 = t16+t18;
t36 = dalpha_Body.*dgamma_Body.*height_Body.*t6.*t9.*4.0;
t37 = (t2.*t3.*width_Body)./2.0;
t38 = (t2.*t8.*width_Body)./2.0;
t39 = (t3.*t6.*width_Body)./2.0;
t40 = (t6.*t8.*width_Body)./2.0;
t42 = dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*t2.*t9;
t43 = height_Body.*t2.*t4.*t13.*2.0;
t44 = height_Body.*t2.*t4.*t15.*2.0;
t48 = (l_F_X.*t3.*t9.*width_Body)./2.0;
t49 = (r_F_X.*t3.*t9.*width_Body)./2.0;
t50 = (height_Body.*t2.*t4.*t12)./2.0;
t52 = (height_Body.*t6.*t9.*t12)./2.0;
t57 = (height_Body.*m_Body.*t9.*t15)./2.0;
t60 = t2.*t8.*t9.*width_Body.*(-1.0./2.0);
t61 = t3.*t6.*t9.*width_Body.*(-1.0./2.0);
t69 = (height_Body.*m_Body.*t4.*t6.*t13)./2.0;
t70 = (height_Body.*m_Body.*t4.*t6.*t15)./2.0;
t71 = l_F_Y.*t2.*t3.*t4.*width_Body.*(-1.0./2.0);
t72 = l_F_Z.*t3.*t4.*t6.*width_Body.*(-1.0./2.0);
t78 = (m_Body.*t4.*t9.*t15.*t16)./4.0;
t29 = r_Tau_Beta_Shoulder+t23;
t30 = -t26;
t45 = -t36;
t46 = 1.0./t32;
t47 = 1.0./t35;
t51 = t9.*t37;
t53 = t9.*t38;
t54 = t9.*t39;
t55 = t9.*t40;
t58 = -t49;
t59 = -t50;
t62 = l_F_Y.*t4.*t37;
t63 = r_F_Y.*t4.*t37;
t64 = l_F_Z.*t4.*t39;
t65 = r_F_Z.*t4.*t39;
t66 = -t57;
t67 = t2.*t6.*t17.*t20;
t68 = t2.*t6.*t17.*t22;
t76 = t2.*t6.*t22.*t31;
t77 = t2.*t6.*t20.*t34;
t80 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t9.*t16.*t19)./2.0;
t81 = (dalpha_Body.*dgamma_Body.*m_Body.*t4.*t9.*t16.*t21)./2.0;
t82 = (m_Body.*t4.*t9.*t13.*t16.*t19)./4.0;
t83 = t19.*t78;
t84 = (m_Body.*t4.*t9.*t13.*t16.*t21)./4.0;
t85 = t21.*t78;
t87 = m_Body.*t4.*t9.*t15.*t16.*t19.*(-1.0./4.0);
t89 = m_Body.*t4.*t9.*t15.*t16.*t21.*(-1.0./4.0);
t92 = t38+t61;
t93 = t39+t60;
t98 = t17+t31+t33+t34;
t105 = l_F_Z+r_F_Z+t28+t42+t69+t70;
t41 = t2.*t29;
t56 = t25+t30;
t74 = -t68;
t75 = l_F_X+r_F_X+t66;
t79 = -t76;
t86 = -t82;
t88 = -t84;
t90 = t37+t55;
t91 = t40+t51;
t96 = l_F_Y.*t92;
t97 = r_F_Y.*t92;
t101 = 1.0./t98;
t102 = t43+t44+t45;
t73 = t6.*t56;
t94 = l_F_Z.*t91;
t95 = r_F_Z.*t91;
t100 = -t96;
t103 = (m_Body.*t102)./4.0;
t106 = t67+t74+t77+t79;
t99 = -t94;
t104 = l_F_Y+r_F_Y+t103;
t108 = t41+t48+t52+t58+t63+t65+t71+t72+t73+t78+t86+t87+t88+t89;
t107 = l_Tau_Alpha_Shoulder+t24+t59+t80+t81+t95+t97+t99+t100;
ddalpha_Body = t27.*t46.*t107.*1.2e+1-height_Body.*t2.*t4.*t27.*t46.*t105.*6.0+height_Body.*t4.*t6.*t27.*t46.*t104.*6.0;
if nargout > 1
    ddbeta_Body = ((t9.*(l_Tau_Alpha_Shoulder+t24)-l_F_Z.*t90+l_F_Y.*t93+r_F_Z.*t90-r_F_Y.*t93-t4.*t6.*t29+t2.*t4.*t56-(l_F_X.*t4.*t8.*width_Body)./2.0+(r_F_X.*t4.*t8.*width_Body)./2.0).*-1.2e+1)./(m_Body.*t14+m_Body.*t18);
end
if nargout > 2
    ddgamma_Body = t27.*t47.*t108.*1.2e+1+height_Body.*t4.*t27.*t47.*t75.*6.0+height_Body.*t2.*t9.*t27.*t47.*t104.*6.0+height_Body.*t6.*t9.*t27.*t47.*t105.*6.0;
end
if nargout > 3
    ddx_Head = t27.*t47.*t75.*(t35+t16.*t20.*3.0)+height_Body.*t4.*t27.*t47.*t108.*6.0+t2.*t4.*t9.*t16.*t27.*t47.*t104.*3.0+t4.*t6.*t9.*t16.*t27.*t47.*t105.*3.0;
end
if nargout > 4
    ddy_Head = t27.*t101.*t105.*t106.*-3.0+t27.*t101.*t104.*(t98+t17.*t19.*t22.*3.0+t17.*t20.*t21.*3.0+t19.*t22.*t31.*3.0+t20.*t21.*t34.*3.0)+height_Body.*t4.*t6.*t27.*t46.*t107.*6.0+height_Body.*t2.*t9.*t27.*t47.*t108.*6.0+t2.*t4.*t9.*t16.*t27.*t47.*t75.*3.0;
end
if nargout > 5
    ddz_Head = t27.*t101.*t104.*t106.*-3.0+t27.*t101.*t105.*(t98+t17.*t19.*t20.*3.0+t17.*t21.*t22.*3.0+t19.*t20.*t34.*3.0+t21.*t22.*t31.*3.0)-height_Body.*t2.*t4.*t27.*t46.*t107.*6.0+height_Body.*t6.*t9.*t27.*t47.*t108.*6.0+t4.*t6.*t9.*t16.*t27.*t47.*t75.*3.0;
end