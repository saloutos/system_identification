function Jparam = singleBodyParamJacobian(in1,in2)
%SINGLEBODYPARAMJACOBIAN
%    JPARAM = SINGLEBODYPARAMJACOBIAN(IN1,IN2)

%    This function was generated by the Symbolic Math Toolbox version 8.4.
%    06-Sep-2021 00:30:12

a = in1(10,:);
d1 = in1(1,:);
d2 = in1(2,:);
d3 = in1(3,:);
pi01 = in2(1,:);
pi02 = in2(2,:);
pi03 = in2(3,:);
pi04 = in2(4,:);
pi05 = in2(5,:);
pi06 = in2(6,:);
pi07 = in2(7,:);
pi08 = in2(8,:);
pi09 = in2(9,:);
pi010 = in2(10,:);
s12 = in1(4,:);
s13 = in1(5,:);
s23 = in1(6,:);
t1 = in1(7,:);
t2 = in1(8,:);
t3 = in1(9,:);
t5 = exp(a);
t6 = exp(d1);
t7 = exp(d2);
t8 = exp(d3);
t9 = a.*2.0;
t10 = d1.*2.0;
t11 = d2.*2.0;
t12 = d3.*2.0;
t13 = s12.^2;
t14 = s13.^2;
t15 = s23.^2;
t16 = t1.^2;
t32 = pi05./2.0;
t33 = pi06./2.0;
t34 = pi07./2.0;
t17 = t5.^2;
t18 = pi04.*s23.*t5;
t19 = pi08.*s23.*t5;
t20 = pi09.*s23.*t5;
t21 = pi01.*t2.*t5;
t22 = pi01.*t3.*t5;
t23 = pi02.*t2.*t5;
t24 = pi02.*t3.*t5;
t25 = pi03.*t2.*t5;
t26 = pi03.*t3.*t5;
t27 = pi04.*t2.*t5;
t28 = pi04.*t3.*t5;
t29 = d1+t9;
t30 = d2+t9;
t31 = d3+t9;
t35 = pi03.*t5.*t7;
t36 = pi04.*t5.*t8;
t37 = pi08.*t5.*t7;
t38 = pi08.*t5.*t8;
t39 = pi09.*t5.*t8;
t40 = pi010.*t5.*t7;
t41 = -t33;
t42 = -t34;
t61 = t9+t10;
t62 = t9+t11;
t63 = t9+t12;
t43 = exp(t29);
t44 = exp(t30);
t45 = exp(t31);
t46 = pi01.*t17;
t47 = pi04.*t17;
t48 = -t19;
t49 = -t23;
t50 = pi03.*s12.*t17;
t52 = pi05.*s12.*t17;
t53 = pi05.*s13.*t17;
t54 = pi06.*s12.*t17;
t55 = pi06.*s13.*t17;
t56 = pi07.*s12.*t17;
t57 = pi07.*s13.*t17;
t64 = exp(t61);
t65 = exp(t62);
t66 = exp(t63);
t67 = -t37;
t68 = -t38;
t69 = -t39;
t72 = pi08.*s12.*t17.*2.0;
t73 = pi08.*s13.*t17.*2.0;
t75 = pi03.*t1.*t17.*2.0;
t77 = pi02.*t6.*t17;
t78 = pi03.*t7.*t17;
t90 = pi05.*t13.*t17;
t91 = pi05.*t14.*t17;
t92 = pi06.*t13.*t17;
t93 = pi06.*t14.*t17;
t94 = pi07.*t13.*t17;
t95 = pi07.*t14.*t17;
t96 = pi08.*s12.*s13.*t17.*4.0;
t128 = t22+t36;
t132 = t32+t34+t41;
t133 = t32+t33+t42;
t134 = t18+t21+t35;
t51 = s13.*t47;
t58 = s23.*t47;
t59 = t1.*t46;
t60 = t2.*t46;
t70 = t50.*2.0;
t76 = t1.*t47.*2.0;
t79 = t8.*t47;
t80 = -t50;
t82 = -t54;
t83 = -t57;
t84 = -t72;
t85 = -t73;
t87 = pi02.*t43.*2.0;
t88 = pi09.*t43.*2.0;
t89 = pi010.*t43.*2.0;
t97 = t1.*t50.*4.0;
t99 = pi05.*t64;
t100 = pi05.*t65;
t101 = pi06.*t64;
t102 = pi05.*t66;
t103 = pi06.*t65;
t104 = pi07.*t64;
t105 = pi06.*t66;
t106 = pi07.*t65;
t107 = pi07.*t66;
t110 = pi09.*s13.*t43.*4.0;
t111 = pi010.*s12.*t43.*4.0;
t113 = pi02.*t1.*t43.*4.0;
t116 = -t77;
t117 = t16.*t46.*2.0;
t118 = -t96;
t119 = pi09.*s13.*t43.*-2.0;
t120 = pi010.*s12.*t43.*-2.0;
t123 = -t92;
t124 = -t95;
t129 = t24+t69;
t130 = t26+t68;
t131 = t5.*t128;
t136 = t20+t40+t49;
t137 = t5.*t134;
t138 = s23.*t5.*t133;
t139 = t5.*t7.*t132;
t140 = t5.*t8.*t133;
t71 = t51.*2.0;
t74 = t59.*2.0;
t81 = -t51;
t86 = -t59;
t98 = t1.*t51.*4.0;
t108 = s13.*t88;
t109 = s12.*t89;
t112 = t1.*t87;
t114 = -t88;
t115 = -t89;
t121 = -t110;
t122 = -t111;
t125 = -t99;
t126 = -t103;
t127 = -t107;
t135 = -t131;
t141 = t28+t140;
t145 = t25+t48+t139;
t146 = t27+t67+t138;
t142 = t5.*t141;
t144 = t70+t71+t74+t87;
t147 = t5.*t146;
t148 = t80+t81+t86+t116;
t149 = t52+t56+t75+t82+t85+t115;
t150 = t53+t55+t76+t83+t84+t114;
t151 = t101+t104+t112+t119+t120+t125;
t143 = -t142;
Jparam = reshape([0.0,t77,0.0,0.0,0.0,t151,t151,0.0,-t5.*t6.*t129,t5.*t6.*t136,0.0,0.0,t78,0.0,t2.*t78+t7.^2.*t17.*t132+t5.*t7.*t145-pi08.*s23.*t7.*t17,0.0,t100+t106+t126-pi08.*s23.*t44.*2.0+pi03.*t2.*t44.*2.0,-t5.*t7.*t130,0.0,-t1.*t78+pi08.*s13.*t7.*t17+pi010.*t6.*t7.*t17-s12.*t7.*t17.*t132,0.0,0.0,0.0,t79,t3.*t79+t8.*t142+t8.^2.*t17.*t133,t102+t105+t127+pi04.*t3.*t45.*2.0,0.0,-t2.*t79+pi08.*t7.*t8.*t17-s23.*t8.*t17.*t133,-t1.*t79+pi08.*s12.*t8.*t17+pi09.*t6.*t8.*t17-s13.*t8.*t17.*t133,0.0,0.0,pi03.*t17,0.0,0.0,0.0,t149,t149,0.0,-t5.*t130,-t5.*t145,0.0,t47,0.0,0.0,0.0,t150,t150,0.0,t143,-t147,0.0,0.0,t47,0.0,t147+t2.*t47-pi08.*t7.*t17+s23.*t17.*t133,0.0,pi08.*t44.*-2.0+t2.*t47.*2.0+pi05.*s23.*t17+pi06.*s23.*t17-pi07.*s23.*t17,t143,0.0,-t1.*t47+pi08.*s12.*t17+pi09.*t6.*t17-s13.*t17.*t133,0.0,t46,0.0,0.0,0.0,t144,t144,0.0,t135,-t137,0.0,0.0,t46,0.0,t58+t60+t78+t137,0.0,t58.*2.0+t60.*2.0+pi03.*t44.*2.0,t135,0.0,t148,0.0,0.0,0.0,t46,t79+t131+t3.*t46,pi04.*t45.*2.0+t3.*t46.*2.0,0.0,-t58-t60-t78,t148,0.0,t46.*2.0,t5.*(pi03.*s12.*t5+pi04.*s13.*t5+pi01.*t1.*t5+pi02.*t5.*t6).*2.0,t137.*2.0,t131.*2.0,s23.*t147.*2.0+t3.*t131.*2.0+t2.*t137.*2.0+t8.*t142.*2.0+t5.*t7.*t145.*2.0,t90+t91+t93+t94+t97+t98+t101+t102+t104+t105+t113+t117+t118+t121+t122+t123+t124+t125+t127+t3.^2.*t46.*2.0+pi04.*t3.*t45.*4.0,t90+t91+t93+t94+t97+t98+t100+t101+t104+t106+t113+t117+t118+t121+t122+t123+t124+t125+t126+t2.*t58.*4.0+t2.*t60.*2.0-pi08.*s23.*t44.*4.0+pi05.*t15.*t17+pi06.*t15.*t17-pi07.*t15.*t17+pi03.*t2.*t44.*4.0,s23.*t142.*-2.0-t2.*t131.*2.0-t5.*t7.*t130.*2.0,s13.*t142.*-2.0-t1.*t131.*2.0-s12.*t5.*t130.*2.0-t5.*t6.*t129.*2.0,s13.*t147.*-2.0-t1.*t137.*2.0-s12.*t5.*t145.*2.0+t5.*t6.*t136.*2.0],[10,10]);
