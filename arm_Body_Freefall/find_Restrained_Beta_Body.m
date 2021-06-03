function beta_Body = find_Restrained_Beta_Body(alpha_Body,gamma_Body,l_Alpha_Hand,l_Beta_Hand,l_X_Fixed,l_Y_Fixed,length_Hand,r_Alpha_Hand,r_Beta_Hand,r_X_Fixed,r_Y_Fixed,width_Body)
%FIND_RESTRAINED_BETA_BODY
%    BETA_BODY = FIND_RESTRAINED_BETA_BODY(ALPHA_BODY,GAMMA_BODY,L_ALPHA_HAND,L_BETA_HAND,L_X_FIXED,L_Y_FIXED,LENGTH_HAND,R_ALPHA_HAND,R_BETA_HAND,R_X_FIXED,R_Y_FIXED,WIDTH_BODY)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    03-Jun-2021 15:16:24

t2 = cos(alpha_Body);
t3 = cos(gamma_Body);
t4 = cos(l_Alpha_Hand);
t5 = cos(r_Alpha_Hand);
t6 = sin(gamma_Body);
t7 = sin(l_Beta_Hand);
t8 = sin(r_Beta_Hand);
t9 = 1.0./width_Body;
t10 = 1.0./t3;
beta_Body = atan2((t9.*t10.*(-l_Y_Fixed.*t3+r_Y_Fixed.*t3+l_X_Fixed.*t2.*t6-r_X_Fixed.*t2.*t6-length_Hand.*t3.*t4.*cos(l_Beta_Hand)+length_Hand.*t3.*t5.*cos(r_Beta_Hand)+length_Hand.*t2.*t4.*t6.*t7+length_Hand.*t2.*t5.*t6.*t8))./sin(alpha_Body),-t9.*t10.*(l_X_Fixed-r_X_Fixed+length_Hand.*t4.*t7+length_Hand.*t5.*t8));