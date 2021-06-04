function [A11,A12,A13] = HFD_Coeffs_L_Tau_Beta_Shoulder_Force(l_Alpha_Hand,l_Beta_Hand,length_Hand)
%HFD_COEFFS_L_TAU_BETA_SHOULDER_FORCE
%    [A11,A12,A13] = HFD_COEFFS_L_TAU_BETA_SHOULDER_FORCE(L_ALPHA_HAND,L_BETA_HAND,LENGTH_HAND)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    04-Jun-2021 18:36:58

t2 = cos(l_Alpha_Hand);
t3 = sin(l_Alpha_Hand);
t4 = sin(l_Beta_Hand);
t5 = t2.^2;
t6 = t3.^2;
t7 = t5+t6;
t8 = 1.0./t7;
A11 = -length_Hand.*t8.*cos(l_Beta_Hand);
if nargout > 1
    A12 = length_Hand.*t2.*t4.*t8;
end
if nargout > 2
    A13 = length_Hand.*t3.*t4.*t8;
end
