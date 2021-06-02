function dotq = ddt_Single(~, q, f2_X,f2_Y,g,l1,m1,tau2)
%DDT_SINGLE この関数の概要をここに記述
%   詳細説明をここに記述

th1 = q(1);
dth1 = q(2);
x1 = q(3);
dx1 = q(4);
y1 = q(5);
dy1 = q(6);

[ddth1,p1ddx,p1ddy] = find_Dds1(dth1,f2_X,f2_Y,g,l1,m1,tau2,th1);

dotq = [dth1, ddth1, dx1, p1ddx, dy1, p1ddy]';
end

