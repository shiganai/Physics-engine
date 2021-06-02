function l_Hip = find_L_Hip(alpha_Body,beta_Body,gamma_Body,height_Body,width_Body,x_Head,y_Head,z_Head)
%FIND_L_HIP
%    L_HIP = FIND_L_HIP(ALPHA_BODY,BETA_BODY,GAMMA_BODY,HEIGHT_BODY,WIDTH_BODY,X_HEAD,Y_HEAD,Z_HEAD)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    02-Jun-2021 20:44:01

t2 = cos(alpha_Body);
t3 = cos(beta_Body);
t4 = cos(gamma_Body);
t5 = sin(alpha_Body);
t6 = sin(beta_Body);
t7 = sin(gamma_Body);
l_Hip = [x_Head-height_Body.*t7-(t3.*t4.*width_Body)./2.0,y_Head-(width_Body.*(t5.*t6+t2.*t3.*t7))./2.0+height_Body.*t2.*t4,z_Head+(width_Body.*(t2.*t6-t3.*t5.*t7))./2.0+height_Body.*t4.*t5];
