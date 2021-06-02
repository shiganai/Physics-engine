function [ddoffset_Body,ddalpha_Body,ddbeta_Body,ddgamma_Body] = find_Dds_Body(alpha_Body,beta_Body,dalpha_Body,depth_Body,dgamma_Body,doffset_Body,g,gamma_Body,height_Body,l_F_X,l_F_Y,l_F_Z,m_Body,offset_Body,r_F_X,r_F_Y,r_F_Z,width_Body)
%FIND_DDS_BODY
%    [DDOFFSET_BODY,DDALPHA_BODY,DDBETA_BODY,DDGAMMA_BODY] = FIND_DDS_BODY(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DEPTH_BODY,DGAMMA_BODY,DOFFSET_BODY,G,GAMMA_BODY,HEIGHT_BODY,L_F_X,L_F_Y,L_F_Z,M_BODY,OFFSET_BODY,R_F_X,R_F_Y,R_F_Z,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    02-Jun-2021 19:11:27

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = sin(alpha_Body);
t6 = sin(beta_Body);
t7 = sin(gamma_Body);
t8 = dalpha_Body.^2;
t9 = depth_Body.^2;
t10 = dgamma_Body.^2;
t11 = height_Body.^2;
t12 = offset_Body.^2;
t13 = width_Body.^2;
t24 = 1.0./m_Body;
t14 = t2.^2;
t15 = t4.^2;
t16 = t5.^2;
t17 = t7.^2;
t18 = l_F_Y.*t2;
t19 = offset_Body.*t2;
t20 = r_F_Y.*t2;
t21 = l_F_Z.*t5;
t22 = offset_Body.*t5;
t23 = r_F_Z.*t5;
t25 = g.*m_Body.*t5;
t26 = t5.*t6;
t27 = m_Body.*t9;
t28 = t2.*t3;
t29 = t2.*t6;
t30 = t3.*t5;
t41 = (l_F_X.*t3.*t7.*width_Body)./2.0;
t42 = (r_F_X.*t3.*t7.*width_Body)./2.0;
t60 = (m_Body.*t4.*t7.*t10.*t11)./4.0;
t31 = -t25;
t32 = t7.*t28;
t33 = t7.*t29;
t34 = t7.*t30;
t35 = t7.*t26;
t36 = t11.*t15.*3.0;
t39 = m_Body.*offset_Body.*t8.*t14;
t40 = m_Body.*offset_Body.*t8.*t16;
t43 = (height_Body.*t7.*t25)./2.0;
t44 = -t41;
t46 = (t3.*t4.*t18.*width_Body)./2.0;
t47 = (t3.*t4.*t20.*width_Body)./2.0;
t48 = (t3.*t4.*t21.*width_Body)./2.0;
t49 = (t3.*t4.*t23.*width_Body)./2.0;
t56 = (height_Body.*m_Body.*t4.*t8.*t14)./2.0;
t57 = (height_Body.*m_Body.*t4.*t10.*t14)./2.0;
t58 = (height_Body.*m_Body.*t4.*t8.*t16)./2.0;
t59 = (height_Body.*m_Body.*t4.*t10.*t16)./2.0;
t64 = -t60;
t66 = (m_Body.*t4.*t7.*t8.*t11.*t14)./4.0;
t67 = t14.*t60;
t68 = (m_Body.*t4.*t7.*t8.*t11.*t16)./4.0;
t69 = t16.*t60;
t37 = -t33;
t38 = -t34;
t45 = -t43;
t50 = t26+t32;
t51 = t28+t35;
t52 = -t47;
t53 = -t49;
t61 = t11+t13+t36;
t62 = (height_Body.*t7.*t39)./2.0;
t63 = (height_Body.*t7.*t40)./2.0;
t72 = t18+t20+t21+t23+t31+t39+t40+t56+t57+t58+t59;
t54 = t29+t38;
t55 = t30+t37;
t65 = (t50.*width_Body)./2.0;
t70 = 1.0./t61;
t73 = t42+t44+t45+t46+t48+t52+t53+t62+t63+t64+t66+t67+t68+t69;
ddoffset_Body = (t70.*t72.*(t61+t11.*t14.*t17.*3.0+t11.*t16.*t17.*3.0))./(m_Body.*t14+m_Body.*t16)-height_Body.*t7.*t24.*t70.*t73.*6.0;
if nargout > 1
    t71 = (t54.*width_Body)./2.0;
    ddalpha_Body = ((l_F_Y.*(t22+t71)-r_F_Z.*(t19+t65)-l_F_Z.*(t19-t65)+r_F_Y.*(t22-t71)+g.*m_Body.*t19+dalpha_Body.*doffset_Body.*m_Body.*offset_Body.*t14.*2.0+dalpha_Body.*doffset_Body.*m_Body.*offset_Body.*t16.*2.0+(g.*height_Body.*m_Body.*t2.*t4)./2.0+dalpha_Body.*doffset_Body.*height_Body.*m_Body.*t4.*t14+dalpha_Body.*doffset_Body.*height_Body.*m_Body.*t4.*t16-dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*offset_Body.*t7.*t14-dalpha_Body.*dgamma_Body.*height_Body.*m_Body.*offset_Body.*t7.*t16-(dalpha_Body.*dgamma_Body.*m_Body.*t4.*t7.*t11.*t14)./2.0-(dalpha_Body.*dgamma_Body.*m_Body.*t4.*t7.*t11.*t16)./2.0).*-1.2e+1)./(t27+m_Body.*t11+m_Body.*t12.*t14.*1.2e+1+m_Body.*t12.*t16.*1.2e+1+m_Body.*t14.*t36+m_Body.*t16.*t36+height_Body.*m_Body.*offset_Body.*t4.*t14.*1.2e+1+height_Body.*m_Body.*offset_Body.*t4.*t16.*1.2e+1);
end
if nargout > 2
    ddbeta_Body = (((l_F_Z.*t51.*width_Body)./2.0-(l_F_Y.*t55.*width_Body)./2.0-(r_F_Z.*t51.*width_Body)./2.0+(r_F_Y.*t55.*width_Body)./2.0+(l_F_X.*t4.*t6.*width_Body)./2.0-(r_F_X.*t4.*t6.*width_Body)./2.0).*1.2e+1)./(t27+m_Body.*t13);
end
if nargout > 3
    ddgamma_Body = t24.*t70.*t73.*-1.2e+1+height_Body.*t7.*t24.*t70.*t72.*6.0;
end
