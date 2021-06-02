
m1 = 1;
g = 1;
l1 = 1;

th1 = 0;
dth1 = 0;

p1 = l1 * [cos(th1), sin(th1)];
v1 = l1 * dth1 * [-sin(th1), cos(th1)];

x1 = p1(1);
y1 = p1(2);

dx1 = v1(1);
dy1 = v1(2);

f2_X = 0;
f2_Y = 0;
tau2 = 0;

q = [th1, dth1, x1, dx1, y1, dy1]';
time = 0:1e-2:10;

[time,q] = ode45(@(t,q) ddt_Single(t, q, f2_X,f2_Y,g,l1,m1,tau2), time, q);

th1 = q(:, 1);
dth1 = q(:, 2);
x1 = q(:, 3);
dx1 = q(:, 4);
y1 = q(:, 5);
dy1 = q(:, 6);

p1 = l1 * [cos(th1), sin(th1)];

zero_Array = zeros(size(time));
x_Array = [zero_Array, p1(:,1)];
z_Array = [zero_Array, p1(:,2)];
y_Array = [zero_Array, zero_Array];

anime = AnimeAndData(time, x_Array, y_Array, z_Array);
plot_Lim = l1 * [-1, 1];
xlim(anime.axAnime, plot_Lim)
ylim(anime.axAnime, plot_Lim)
zlim(anime.axAnime, plot_Lim)
view(anime.axAnime, [0,-1,0])

x_Array = [zero_Array, x1];
z_Array = [zero_Array, y1];
y_Array = [zero_Array, zero_Array];

anime = AnimeAndData(time, x_Array, y_Array, z_Array);
plot_Lim = l1 * [-1, 1];
xlim(anime.axAnime, plot_Lim)
ylim(anime.axAnime, plot_Lim)
zlim(anime.axAnime, plot_Lim)
view(anime.axAnime, [0,-1,0])












































