function dl_Shoulder = FRD_Dl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,dgamma_Body,dx_Head,dy_Head,dz_Head,gamma_Body,width_Body)
%FRD_DL_SHOULDER
%    DL_SHOULDER = FRD_DL_SHOULDER(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DBETA_BODY,DGAMMA_BODY,DX_HEAD,DY_HEAD,DZ_HEAD,GAMMA_BODY,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    05-Jun-2021 09:07:49

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = sin(alpha_Body);
t6 = sin(beta_Body);
t7 = sin(gamma_Body);
dl_Shoulder = [dx_Head+(dbeta_Body.*t4.*t6.*width_Body)./2.0+(dgamma_Body.*t3.*t7.*width_Body)./2.0,dy_Head-(dalpha_Body.*t2.*t6.*width_Body)./2.0-(dbeta_Body.*t3.*t5.*width_Body)./2.0+(dalpha_Body.*t3.*t5.*t7.*width_Body)./2.0+(dbeta_Body.*t2.*t6.*t7.*width_Body)./2.0-(dgamma_Body.*t2.*t3.*t4.*width_Body)./2.0,dz_Head-(dalpha_Body.*t5.*t6.*width_Body)./2.0+(dbeta_Body.*t2.*t3.*width_Body)./2.0-(dalpha_Body.*t2.*t3.*t7.*width_Body)./2.0+(dbeta_Body.*t5.*t6.*t7.*width_Body)./2.0-(dgamma_Body.*t3.*t4.*t5.*width_Body)./2.0];
