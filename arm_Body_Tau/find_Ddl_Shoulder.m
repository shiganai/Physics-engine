function ddl_Shoulder = find_Ddl_Shoulder(alpha_Body,beta_Body,dalpha_Body,dbeta_Body,ddalpha_Body,ddbeta_Body,ddgamma_Body,ddx_Head,ddy_Head,ddz_Head,dgamma_Body,gamma_Body,width_Body)
%FIND_DDL_SHOULDER
%    DDL_SHOULDER = FIND_DDL_SHOULDER(ALPHA_BODY,BETA_BODY,DALPHA_BODY,DBETA_BODY,DDALPHA_BODY,DDBETA_BODY,DDGAMMA_BODY,DDX_HEAD,DDY_HEAD,DDZ_HEAD,DGAMMA_BODY,GAMMA_BODY,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 15:17:31

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = sin(alpha_Body);
t6 = sin(beta_Body);
t7 = sin(gamma_Body);
t8 = dalpha_Body.^2;
t9 = dbeta_Body.^2;
t10 = dgamma_Body.^2;
ddl_Shoulder = [ddx_Head+(ddbeta_Body.*t4.*t6.*width_Body)./2.0+(ddgamma_Body.*t3.*t7.*width_Body)./2.0+(t3.*t4.*t9.*width_Body)./2.0+(t3.*t4.*t10.*width_Body)./2.0-dbeta_Body.*dgamma_Body.*t6.*t7.*width_Body,ddy_Head-(ddalpha_Body.*t2.*t6.*width_Body)./2.0-(ddbeta_Body.*t3.*t5.*width_Body)./2.0+(t5.*t6.*t8.*width_Body)./2.0+(t5.*t6.*t9.*width_Body)./2.0-dalpha_Body.*dbeta_Body.*t2.*t3.*width_Body+(ddalpha_Body.*t3.*t5.*t7.*width_Body)./2.0+(ddbeta_Body.*t2.*t6.*t7.*width_Body)./2.0-(ddgamma_Body.*t2.*t3.*t4.*width_Body)./2.0+(t2.*t3.*t7.*t8.*width_Body)./2.0+(t2.*t3.*t7.*t9.*width_Body)./2.0+(t2.*t3.*t7.*t10.*width_Body)./2.0-dalpha_Body.*dbeta_Body.*t5.*t6.*t7.*width_Body+dalpha_Body.*dgamma_Body.*t3.*t4.*t5.*width_Body+dbeta_Body.*dgamma_Body.*t2.*t4.*t6.*width_Body,ddz_Head-(ddalpha_Body.*t5.*t6.*width_Body)./2.0+(ddbeta_Body.*t2.*t3.*width_Body)./2.0-(t2.*t6.*t8.*width_Body)./2.0-(t2.*t6.*t9.*width_Body)./2.0-dalpha_Body.*dbeta_Body.*t3.*t5.*width_Body-(ddalpha_Body.*t2.*t3.*t7.*width_Body)./2.0+(ddbeta_Body.*t5.*t6.*t7.*width_Body)./2.0-(ddgamma_Body.*t3.*t4.*t5.*width_Body)./2.0+(t3.*t5.*t7.*t8.*width_Body)./2.0+(t3.*t5.*t7.*t9.*width_Body)./2.0+(t3.*t5.*t7.*t10.*width_Body)./2.0+dalpha_Body.*dbeta_Body.*t2.*t6.*t7.*width_Body-dalpha_Body.*dgamma_Body.*t2.*t3.*t4.*width_Body+dbeta_Body.*dgamma_Body.*t4.*t5.*t6.*width_Body];
