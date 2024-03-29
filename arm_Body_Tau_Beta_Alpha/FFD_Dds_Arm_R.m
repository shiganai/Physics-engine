function [ddr_Alpha_Hand,ddr_Beta_Hand] = FFD_Dds_Arm_R(dr_Beta_Hand,dr_Alpha_Hand,g,length_Hand,m_Hand,r_Alpha_Hand,r_Beta_Hand,r_F_X,r_F_Y,r_F_Z,r_Tau_Beta_Shoulder,r_Tau_Alpha_Shoulder)
%FFD_DDS_ARM_R
%    [DDR_ALPHA_HAND,DDR_BETA_HAND] = FFD_DDS_ARM_R(DR_BETA_HAND,DR_ALPHA_HAND,G,LENGTH_HAND,M_HAND,R_ALPHA_HAND,R_BETA_HAND,R_F_X,R_F_Y,R_F_Z,R_TAU_BETA_SHOULDER,R_TAU_ALPHA_SHOULDER)

%    This function was generated by the Symbolic Math Toolbox version 8.6.
%    12-Jun-2021 22:17:11

t2 = cos(r_Alpha_Hand);
t3 = cos(r_Beta_Hand);
t4 = sin(r_Alpha_Hand);
t5 = sin(r_Beta_Hand);
t6 = dr_Beta_Hand.^2;
t7 = dr_Alpha_Hand.^2;
t8 = length_Hand.^2;
t32 = -r_Tau_Alpha_Shoulder;
t9 = t2.^2;
t10 = t2.^3;
t12 = t2.^5;
t14 = t2.^7;
t15 = t3.^2;
t16 = t3.^3;
t18 = t3.^5;
t20 = t3.^7;
t21 = t4.^2;
t22 = t4.^3;
t24 = t4.^5;
t26 = t4.^7;
t27 = t5.^2;
t28 = t5.^3;
t30 = t5.^5;
t33 = length_Hand.*r_F_X.*t3;
t37 = t3.^11;
t43 = length_Hand.*r_F_Z.*t4.*t5;
t45 = length_Hand.*r_F_Z.*t2.*t3;
t46 = length_Hand.*r_F_Y.*t2.*t5;
t47 = length_Hand.*r_F_Y.*t3.*t4;
t52 = (g.*length_Hand.*m_Hand.*t2.*t3)./2.0;
t53 = (g.*length_Hand.*m_Hand.*t4.*t5)./2.0;
t70 = (m_Hand.*t3.*t5.*t6.*t8)./4.0;
t11 = t9.^2;
t13 = t9.^3;
t17 = t15.^2;
t19 = t15.^3;
t23 = t21.^2;
t25 = t21.^3;
t29 = t27.^2;
t31 = t27.^3;
t36 = t16.^3;
t40 = t27.^5;
t41 = r_Tau_Beta_Shoulder.*t9;
t42 = r_Tau_Beta_Shoulder.*t21;
t51 = -t45;
t55 = t2.*t3.*t5.*t22.*2.0;
t56 = -t52;
t58 = t9.*t21.*t28.*2.0;
t60 = t2.*t3.*t22.*t28.*2.0;
t61 = t2.*t5.*t16.*t22.*2.0;
t62 = t3.*t4.*t10.*t28.*2.0;
t63 = t4.*t5.*t10.*t16.*2.0;
t64 = t3.*t4.*t10.*t30.*2.0;
t65 = t4.*t5.*t10.*t18.*2.0;
t77 = m_Hand.*t2.*t8.*t16.*t26.*4.0;
t78 = m_Hand.*t2.*t8.*t18.*t26.*4.0;
t79 = t5.*t9.*t15.*t21.*2.0;
t80 = t4.*t10.*t16.*t28.*4.0;
t84 = m_Hand.*t2.*t8.*t16.*t24.*1.2e+1;
t86 = m_Hand.*t2.*t8.*t18.*t24.*1.2e+1;
t87 = m_Hand.*t4.*t8.*t12.*t18.*1.2e+1;
t88 = m_Hand.*t4.*t8.*t14.*t37.*4.0;
t89 = m_Hand.*t4.*t8.*t12.*t20.*1.2e+1;
t91 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t9)./2.0;
t92 = t9.*t15.*t21.*t27.*2.0;
t95 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t5.*t8.*t21)./2.0;
t99 = m_Hand.*t2.*t3.*t8.*t26.*t27.*4.0;
t104 = m_Hand.*t8.*t10.*t16.*t22.*1.2e+1;
t105 = m_Hand.*t8.*t10.*t18.*t24.*1.2e+1;
t106 = m_Hand.*t8.*t10.*t20.*t22.*1.2e+1;
t107 = m_Hand.*t8.*t10.*t20.*t24.*1.2e+1;
t108 = m_Hand.*t8.*t12.*t20.*t22.*1.2e+1;
t119 = t9.*t70;
t120 = (m_Hand.*t3.*t5.*t7.*t8.*t9)./4.0;
t126 = (m_Hand.*t2.*t5.*t7.*t8.*t22)./6.0;
t127 = t21.*t70;
t128 = (m_Hand.*t3.*t5.*t7.*t8.*t21)./4.0;
t143 = m_Hand.*t3.*t5.*t6.*t8.*t9.*(-1.0./4.0);
t146 = (m_Hand.*t2.*t5.*t6.*t8.*t22)./2.4e+1;
t147 = m_Hand.*t3.*t5.*t6.*t8.*t21.*(-1.0./4.0);
t150 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t21)./4.0;
t155 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t2.*t8.*t22.*t27)./6.0;
t156 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t27)./6.0;
t164 = m_Hand.*t2.*t8.*t16.*t26.*t27.*8.0;
t166 = m_Hand.*t4.*t8.*t12.*t16.*t27.*1.2e+1;
t172 = m_Hand.*t4.*t8.*t14.*t20.*t27.*1.6e+1;
t173 = m_Hand.*t4.*t8.*t12.*t18.*t27.*2.4e+1;
t180 = (m_Hand.*t2.*t7.*t8.*t15.*t22)./6.0;
t188 = (m_Hand.*t2.*t6.*t8.*t22.*t27)./6.0;
t192 = (m_Hand.*t2.*t7.*t8.*t22.*t28)./6.0;
t193 = (m_Hand.*t4.*t7.*t8.*t10.*t28)./6.0;
t194 = (m_Hand.*t4.*t7.*t8.*t10.*t30)./6.0;
t201 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t9.*t16.*t21)./2.0;
t202 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t9.*t18.*t21)./4.0;
t217 = (m_Hand.*t3.*t6.*t8.*t9.*t21)./8.0;
t218 = (m_Hand.*t2.*t6.*t8.*t15.*t22)./1.2e+1;
t233 = (m_Hand.*t4.*t6.*t8.*t10.*t27)./1.2e+1;
t235 = (m_Hand.*t2.*t6.*t8.*t22.*t28)./2.4e+1;
t236 = (m_Hand.*t4.*t6.*t8.*t10.*t28)./2.4e+1;
t237 = (m_Hand.*t4.*t6.*t8.*t10.*t30)./2.4e+1;
t243 = (m_Hand.*t6.*t8.*t9.*t16.*t21)./3.0;
t260 = m_Hand.*t8.*t10.*t16.*t24.*t27.*2.4e+1;
t261 = m_Hand.*t8.*t10.*t18.*t22.*t27.*2.4e+1;
t266 = m_Hand.*t8.*t10.*t18.*t24.*t27.*3.6e+1;
t268 = m_Hand.*t8.*t12.*t18.*t22.*t27.*3.6e+1;
t271 = m_Hand.*t8.*t12.*t20.*t22.*t27.*4.8e+1;
t274 = (m_Hand.*t4.*t5.*t7.*t8.*t10.*t15)./6.0;
t283 = (m_Hand.*t6.*t8.*t9.*t18.*t21)./8.0;
t291 = m_Hand.*t8.*t9.*t15.*t21.*t27.*7.2e+1;
t293 = (m_Hand.*t3.*t5.*t6.*t8.*t9.*t21)./1.2e+1;
t295 = (m_Hand.*t4.*t5.*t6.*t8.*t10.*t15)./2.4e+1;
t297 = m_Hand.*t2.*t5.*t7.*t8.*t15.*t22.*(-1.0./6.0);
t299 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t21.*t27)./2.0;
t300 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t2.*t8.*t15.*t22.*t27)./3.0;
t301 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t9.*t21.*t28)./3.0;
t302 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t5.*t8.*t9.*t16.*t21)./3.0;
t303 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t15.*t27)./3.0;
t312 = (m_Hand.*t3.*t6.*t8.*t9.*t21.*t27)./3.0;
t314 = (m_Hand.*t3.*t6.*t8.*t9.*t21.*t28)./6.0;
t315 = (m_Hand.*t5.*t6.*t8.*t9.*t16.*t21)./6.0;
t316 = (m_Hand.*t4.*t6.*t8.*t10.*t15.*t27)./6.0;
t319 = m_Hand.*t3.*t7.*t8.*t9.*t21.*t28.*(2.0./3.0);
t320 = m_Hand.*t5.*t7.*t8.*t9.*t16.*t21.*(2.0./3.0);
t321 = (m_Hand.*t4.*t7.*t8.*t10.*t15.*t28)./3.0;
t324 = (m_Hand.*t3.*t7.*t8.*t9.*t21.*t30)./4.0;
t325 = (m_Hand.*t5.*t7.*t8.*t9.*t18.*t21)./4.0;
t326 = (m_Hand.*t4.*t7.*t8.*t10.*t15.*t27)./6.0;
t336 = (m_Hand.*t4.*t6.*t8.*t10.*t15.*t28)./1.2e+1;
t338 = (m_Hand.*t3.*t6.*t8.*t9.*t21.*t30)./1.2e+1;
t339 = (m_Hand.*t5.*t6.*t8.*t9.*t18.*t21)./1.2e+1;
t344 = (m_Hand.*t6.*t8.*t9.*t16.*t21.*t27)./4.0;
t345 = (m_Hand.*t6.*t8.*t9.*t16.*t21.*t28)./6.0;
t346 = (m_Hand.*t7.*t8.*t9.*t16.*t21.*t28)./2.0;
t34 = t11.^2;
t35 = t17.^2;
t38 = t23.^2;
t39 = t29.^2;
t44 = t5.*t23;
t48 = -t41;
t49 = -t42;
t50 = t11.*t30;
t54 = t5.*t11.*t17;
t57 = t11.*t15.*t28.*2.0;
t66 = m_Hand.*t8.*t9.*t17.*3.6e+1;
t68 = m_Hand.*t8.*t17.*t23.*1.2e+1;
t69 = m_Hand.*t8.*t17.*t21.*3.6e+1;
t71 = -t60;
t72 = -t61;
t73 = -t64;
t74 = -t65;
t82 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t23)./1.2e+1;
t83 = -t77;
t85 = m_Hand.*t4.*t8.*t14.*t36.*4.0;
t90 = -t80;
t93 = m_Hand.*t8.*t9.*t25.*t29.*4.0;
t94 = m_Hand.*t8.*t11.*t23.*t31.*6.0;
t97 = (m_Hand.*t3.*t6.*t8.*t23)./1.2e+1;
t100 = m_Hand.*t2.*t3.*t8.*t26.*t29.*4.0;
t101 = -t84;
t103 = -t87;
t109 = m_Hand.*t8.*t9.*t19.*t21.*2.4e+1;
t112 = m_Hand.*t8.*t11.*t17.*t29.*1.2e+1;
t114 = m_Hand.*t8.*t13.*t15.*t31.*1.2e+1;
t115 = m_Hand.*t8.*t13.*t19.*t27.*1.2e+1;
t116 = m_Hand.*t8.*t11.*t19.*t27.*2.4e+1;
t117 = m_Hand.*t8.*t13.*t17.*t29.*2.4e+1;
t118 = m_Hand.*t8.*t11.*t15.*t27.*3.6e+1;
t124 = m_Hand.*t8.*t15.*t25.*t27.*1.2e+1;
t125 = m_Hand.*t8.*t15.*t23.*t27.*3.6e+1;
t130 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t11.*t16)./1.2e+1;
t131 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t11.*t18)./1.2e+1;
t132 = -t99;
t134 = m_Hand.*t3.*t4.*t8.*t14.*t40.*4.0;
t135 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t16.*t23)./1.2e+1;
t136 = -t104;
t137 = -t105;
t138 = -t108;
t139 = m_Hand.*t8.*t12.*t22.*t36.*1.2e+1;
t144 = -t120;
t148 = -t128;
t152 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t5.*t8.*t11.*t18)./6.0;
t153 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t2.*t8.*t22.*t29)./3.0;
t154 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t29)./3.0;
t157 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t2.*t8.*t22.*t31)./6.0;
t158 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t31)./6.0;
t159 = (m_Hand.*t6.*t8.*t11.*t18)./1.2e+1;
t160 = (m_Hand.*t6.*t8.*t11.*t16)./2.4e+1;
t163 = (m_Hand.*t6.*t8.*t16.*t23)./2.4e+1;
t165 = t27.*t84;
t167 = m_Hand.*t3.*t8.*t10.*t24.*t29.*1.2e+1;
t168 = m_Hand.*t4.*t8.*t12.*t16.*t29.*1.2e+1;
t169 = m_Hand.*t3.*t8.*t10.*t24.*t31.*1.2e+1;
t170 = m_Hand.*t3.*t8.*t12.*t22.*t31.*1.2e+1;
t171 = m_Hand.*t4.*t8.*t14.*t16.*t31.*1.6e+1;
t174 = m_Hand.*t4.*t8.*t14.*t18.*t29.*2.4e+1;
t175 = m_Hand.*t4.*t8.*t14.*t18.*t31.*4.0e+1;
t176 = m_Hand.*t4.*t8.*t14.*t20.*t29.*4.0e+1;
t179 = -t146;
t181 = (m_Hand.*t2.*t7.*t8.*t17.*t22)./6.0;
t182 = (m_Hand.*t4.*t7.*t8.*t10.*t17)./6.0;
t183 = (m_Hand.*t4.*t7.*t8.*t10.*t19)./6.0;
t185 = (m_Hand.*t5.*t7.*t8.*t11.*t18)./6.0;
t186 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t11.*t27)./1.2e+1;
t187 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t11.*t29)./1.2e+1;
t189 = (m_Hand.*t2.*t6.*t8.*t22.*t29)./4.0;
t190 = (m_Hand.*t4.*t6.*t8.*t10.*t29)./4.0;
t191 = (m_Hand.*t4.*t6.*t8.*t10.*t31)./6.0;
t196 = -t156;
t203 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t11.*t16.*t28)./3.0;
t204 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t8.*t11.*t16.*t27)./6.0;
t205 = -t166;
t210 = -t172;
t212 = m_Hand.*t4.*t8.*t14.*t27.*t36.*2.0e+1;
t214 = m_Hand.*t8.*t9.*t15.*t25.*t27.*4.0;
t215 = m_Hand.*t8.*t11.*t17.*t23.*t27.*6.0;
t216 = m_Hand.*t8.*t13.*t19.*t21.*t27.*4.0;
t219 = (m_Hand.*t2.*t6.*t8.*t17.*t22)./1.2e+1;
t220 = (m_Hand.*t4.*t6.*t8.*t10.*t17)./1.2e+1;
t221 = (m_Hand.*t4.*t6.*t8.*t10.*t19)./1.2e+1;
t222 = -t180;
t224 = (m_Hand.*t3.*t6.*t8.*t11.*t29)./1.2e+1;
t225 = (m_Hand.*t3.*t6.*t8.*t11.*t27)./2.4e+1;
t228 = (m_Hand.*t3.*t7.*t8.*t11.*t28)./1.2e+1;
t229 = (m_Hand.*t5.*t7.*t8.*t11.*t16)./1.2e+1;
t231 = -t188;
t234 = (m_Hand.*t2.*t6.*t8.*t22.*t31)./1.2e+1;
t238 = -t192;
t239 = -t194;
t240 = (m_Hand.*t3.*t6.*t8.*t23.*t27)./2.4e+1;
t241 = (m_Hand.*t3.*t7.*t8.*t23.*t28)./1.2e+1;
t244 = (m_Hand.*t6.*t8.*t11.*t16.*t27)./6.0;
t245 = (m_Hand.*t7.*t8.*t11.*t16.*t28)./3.0;
t246 = -t201;
t251 = m_Hand.*t8.*t11.*t15.*t21.*t27.*1.2e+1;
t252 = t29.*t104;
t253 = m_Hand.*t8.*t9.*t15.*t23.*t31.*1.2e+1;
t254 = m_Hand.*t8.*t9.*t19.*t23.*t27.*1.2e+1;
t255 = m_Hand.*t8.*t11.*t15.*t23.*t29.*1.2e+1;
t257 = m_Hand.*t8.*t13.*t17.*t21.*t29.*1.2e+1;
t258 = m_Hand.*t8.*t9.*t15.*t23.*t27.*2.4e+1;
t259 = m_Hand.*t8.*t9.*t17.*t21.*t27.*2.4e+1;
t262 = m_Hand.*t8.*t9.*t17.*t23.*t29.*2.4e+1;
t263 = m_Hand.*t8.*t11.*t15.*t21.*t31.*2.4e+1;
t265 = m_Hand.*t8.*t10.*t16.*t24.*t29.*3.6e+1;
t267 = m_Hand.*t8.*t12.*t16.*t22.*t29.*3.6e+1;
t269 = m_Hand.*t8.*t11.*t17.*t21.*t29.*4.8e+1;
t270 = m_Hand.*t8.*t12.*t16.*t22.*t31.*4.8e+1;
t272 = t21.*t120;
t273 = t15.*t126;
t276 = -t217;
t281 = -t236;
t286 = m_Hand.*t8.*t13.*t15.*t21.*t31.*-1.2e+1;
t288 = -t260;
t290 = -t268;
t292 = m_Hand.*t8.*t12.*t18.*t22.*t29.*7.2e+1;
t294 = t15.*t146;
t296 = (m_Hand.*t4.*t5.*t6.*t8.*t10.*t17)./2.4e+1;
t298 = m_Hand.*t4.*t5.*t7.*t8.*t10.*t17.*(-1.0./6.0);
t305 = t29.*t150;
t307 = t17.*t155;
t308 = t17.*t156;
t309 = -t283;
t310 = -t293;
t311 = -t295;
t313 = t15.*t188;
t318 = (m_Hand.*t2.*t6.*t8.*t15.*t22.*t29)./6.0;
t322 = (m_Hand.*t4.*t7.*t8.*t10.*t17.*t27)./3.0;
t323 = t27.*t180;
t327 = (m_Hand.*t4.*t7.*t8.*t10.*t15.*t29)./6.0;
t328 = -t299;
t329 = -t300;
t330 = dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t15.*t29.*(-1.0./3.0);
t331 = dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t4.*t8.*t10.*t17.*t27.*(-1.0./6.0);
t332 = t27.*t201;
t333 = -t316;
t335 = t29.*t217;
t340 = -t319;
t341 = -t320;
t342 = -t321;
t343 = -t326;
t347 = m_Hand.*t3.*t6.*t8.*t9.*t21.*t29.*(-1.0./8.0);
t348 = m_Hand.*t2.*t6.*t8.*t17.*t22.*t27.*(-1.0./1.2e+1);
t349 = -t338;
t350 = -t339;
t351 = -t344;
t352 = -t345;
t59 = m_Hand.*t8.*t27.*t38;
t67 = m_Hand.*t8.*t34.*t40;
t75 = m_Hand.*t8.*t11.*t35.*1.2e+1;
t96 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t44)./6.0;
t98 = -t82;
t102 = -t85;
t110 = m_Hand.*t8.*t27.*t34.*t35;
t111 = m_Hand.*t8.*t19.*t29.*t34.*4.0;
t113 = m_Hand.*t8.*t17.*t31.*t34.*6.0;
t121 = -t93;
t122 = -t94;
t123 = m_Hand.*t8.*t13.*t21.*t39.*4.0;
t129 = (m_Hand.*t3.*t7.*t8.*t44)./6.0;
t133 = m_Hand.*t3.*t4.*t8.*t14.*t39.*4.0;
t140 = m_Hand.*t8.*t15.*t34.*t39.*4.0;
t151 = (dr_Beta_Hand.*dr_Alpha_Hand.*m_Hand.*t3.*t8.*t50)./6.0;
t161 = -t131;
t184 = (m_Hand.*t3.*t7.*t8.*t50)./6.0;
t195 = -t153;
t197 = -t158;
t198 = t27.*t82;
t199 = -t160;
t200 = -t163;
t206 = -t167;
t207 = -t170;
t208 = m_Hand.*t3.*t8.*t12.*t22.*t39.*1.2e+1;
t209 = -t171;
t211 = m_Hand.*t4.*t8.*t14.*t16.*t39.*2.0e+1;
t213 = -t174;
t223 = -t182;
t227 = -t185;
t230 = -t187;
t232 = -t190;
t242 = (m_Hand.*t7.*t8.*t16.*t44)./1.2e+1;
t247 = -t204;
t248 = -t214;
t249 = -t215;
t250 = -t216;
t256 = t21.*t114;
t264 = t21.*t116;
t275 = t5.*t182;
t277 = -t219;
t278 = -t221;
t279 = -t225;
t280 = -t234;
t282 = -t240;
t284 = -t245;
t285 = -t255;
t287 = -t257;
t289 = -t267;
t304 = t15.*t153;
t306 = t15.*t154;
t317 = t15.*t190;
t334 = -t318;
t337 = t27.*t219;
t353 = t44+t50+t54+t55+t57+t58+t62+t63+t71+t72+t73+t74+t79+t90;
t76 = -t59;
t81 = -t67;
t141 = -t111;
t142 = -t113;
t145 = -t123;
t149 = -t129;
t162 = -t133;
t177 = -t110;
t178 = -t140;
t226 = -t184;
t354 = t32+t47+t51+t56+t91+t95+t96+t97+t151+t152+t159+t181+t183+t189+t191+t199+t200+t203+t218+t220+t222+t223+t224+t231+t232+t233+t243+t244+t276+t277+t278+t279+t280+t282+t301+t302+t309+t312+t313+t317+t322+t323+t327+t333+t334+t343+t347+t348+t351;
t355 = t33+t43+t46+t48+t49+t53+t70+t98+t126+t130+t135+t143+t144+t147+t148+t149+t150+t154+t155+t157+t161+t179+t186+t193+t195+t196+t197+t198+t202+t226+t227+t228+t229+t230+t235+t237+t238+t239+t241+t242+t246+t247+t272+t274+t281+t284+t294+t296+t297+t298+t303+t304+t305+t307+t310+t311+t314+t315+t324+t325+t328+t329+t330+t331+t332+t336+t340+t341+t342+t346+t349+t350+t352;
t356 = t66+t68+t69+t75+t76+t78+t81+t83+t86+t88+t89+t100+t101+t102+t103+t106+t107+t109+t112+t114+t115+t116+t117+t118+t121+t122+t124+t125+t132+t134+t136+t137+t138+t139+t141+t142+t145+t162+t164+t165+t168+t169+t173+t175+t176+t177+t178+t205+t206+t207+t208+t209+t210+t211+t212+t213+t248+t249+t250+t251+t252+t253+t254+t258+t259+t261+t262+t263+t264+t265+t266+t269+t270+t271+t285+t286+t287+t288+t289+t290+t291+t292;
t357 = 1.0./t356;
ddr_Alpha_Hand = t354.*t357.*(t15.*3.0-t92+t9.*t27.*3.0+t21.*t27.*3.0-t2.*t3.*t22-t4.*t10.*t16+t4.*t10.*t18+t2.*t16.*t22+t9.*t21.*t27-t9.*t21.*t29.*2.0+t9.*t21.*t31-t3.*t4.*t10.*t27+t3.*t4.*t10.*t29+t2.*t3.*t22.*t27+t4.*t10.*t16.*t27.*2.0+t9.*t15.*t21.*t29.*2.0+t9.*t17.*t21.*t27).*4.8e+1+t353.*t355.*t357.*2.4e+1;
if nargout > 1
    ddr_Beta_Hand = t355.*t357.*(t92+t9.*t15.*3.0+t11.*t19+t15.*t21.*3.0+t15.*t23+t9.*t17.*t21.*2.0+t11.*t15.*t29+t11.*t17.*t27.*2.0).*4.8e+1+t353.*t354.*t357.*2.4e+1;
end
