function dotq = ddt_Double(~, q,g,l1,l2,m1,m2,tau2)
%DDT_SINGLE この関数の概要をここに記述
%   詳細説明をここに記述

th1 = q(1);
dth1 = q(2);
th2 = q(3);
dth2 = q(4);
x2 = q(5);
dx2 = q(6);
y2 = q(7);
dy2 = q(8);

f2_X = 0;
f2_Y = 0;

v1 = l1 * dth1 * [-sin(th1), cos(th1)];
dx2 = v1(1);
dy2 = v1(2);

[~,p1ddx_0,p1ddy_0] = find_Dds1(dth1,f2_X,f2_Y,g,l1,m1,tau2,th1);
[ddx2_0,ddy2_0,~] = find_Dds2(dth2,f2_X,f2_Y,g,l2,m2,tau2,th2);

f2_X = 1;
f2_Y = 0;

[~,p1ddx_1,p1ddy_1] = find_Dds1(dth1,f2_X,f2_Y,g,l1,m1,tau2,th1);
[ddx2_1,ddy2_1,~] = find_Dds2(dth2,f2_X,f2_Y,g,l2,m2,tau2,th2);

xfx = (p1ddx_1 - p1ddx_0) - (ddx2_1 - ddx2_0);
yfx = (p1ddy_1 - p1ddy_0) - (ddy2_1 - ddy2_0);

f2_X = 0;
f2_Y = 1;

[~,p1ddx_1,p1ddy_1] = find_Dds1(dth1,f2_X,f2_Y,g,l1,m1,tau2,th1);
[ddx2_1,ddy2_1,~] = find_Dds2(dth2,f2_X,f2_Y,g,l2,m2,tau2,th2);

xfy = (p1ddx_1 - p1ddx_0) - (ddx2_1 - ddx2_0);
yfy = (p1ddy_1 - p1ddy_0) - (ddy2_1 - ddy2_0);

F = inv([xfx, xfy; yfx, yfy]) * (-[p1ddx_0 - ddx2_0, p1ddy_0 - ddy2_0]');

f2_X = F(1);
f2_Y = F(2);

[ddth1,~,~] = find_Dds1(dth1,f2_X,f2_Y,g,l1,m1,tau2,th1);
[ddx2,ddy2,ddth2] = find_Dds2(dth2,f2_X,f2_Y,g,l2,m2,tau2,th2);

dotq = [dth1, ddth1, dth2, ddth2, dx2, ddx2, dy2, ddy2]';
end






































