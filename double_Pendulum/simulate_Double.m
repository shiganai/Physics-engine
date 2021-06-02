
clear all

g = 1;

m1 = 1;
l1 = 1;

m2 = 1;
l2 = 1;

th1 = 0;
dth1 = 0;

th2 = 0;
dth2 = 0;

p1 = l1 * [cos(th1), sin(th1)];
v1 = l1 * dth1 * [-sin(th1), cos(th1)];

x2 = p1(1);
y2 = p1(2);

dx2 = v1(1);
dy2 = v1(2);

f2_X = 0;
f2_Y = 0;
tau2 = 0.1;

q = [th1, dth1, th2, dth2, x2, dx2, y2, dy2]';
time = 0:1e-2:100;

[time,q] = ode45(@(t,q) ddt_Double(t, q,g,l1,l2,m1,m2,tau2), time, q);

th1 = q(:, 1);
dth1 = q(:, 2);
th2 = q(:, 3);
dth2 = q(:, 4);
x2 = q(:, 5);
dx2 = q(:, 6);
y2 = q(:, 7);
dy2 = q(:, 8);

p1 = l1 * [cos(th1), sin(th1)];
p20 = [x2, y2];
p2 = p20 + l2 * [cos(th2), sin(th2)];

zero_Array = zeros(size(time));
nan_Array = nan(size(time));
x_Array = [zero_Array, p1(:,1), nan_Array, p20(:,1), p2(:,1)];
z_Array = [zero_Array, p1(:,2), nan_Array, p20(:,2), p2(:,2)];
y_Array = [zero_Array, zero_Array, nan_Array, zero_Array, zero_Array];

anime = AnimeAndData(time, x_Array, y_Array, z_Array);
plot_Lim = (l1 + l2) * [-1, 1];
xlim(anime.axAnime, plot_Lim)
ylim(anime.axAnime, plot_Lim)
zlim(anime.axAnime, plot_Lim)
view(anime.axAnime, [0,-1,0])





















































