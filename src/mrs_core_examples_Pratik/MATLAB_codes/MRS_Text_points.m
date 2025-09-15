clc
clear all
close all

pt1_M = [0;0;2];
pt2_M = [0;0;6];
pt3_M = [0;-2;4];
pt4_M = [0;-4;6];
pt5_M = [0;-4;2];

plot3(pt1_M(1),pt1_M(2),pt1_M(3),"b.","MarkerSize",10); hold on;
plot3(pt2_M(1),pt2_M(2),pt2_M(3),"b.","MarkerSize",10); hold on;
plot3(pt3_M(1),pt3_M(2),pt3_M(3),"b.","MarkerSize",10); hold on;
plot3(pt4_M(1),pt4_M(2),pt4_M(3),"b.","MarkerSize",10); hold on;
plot3(pt5_M(1),pt5_M(2),pt5_M(3),"b.","MarkerSize",10); hold on;
grid on
ylabel("y")
view(-90,0)
axis equal
daspect([1 1 1])
ylim([-15  1])
zlim([0 7])


pt1_R = [0;-6;2];
pt2_R = [0;-6;6];
pt3_R = [0;-9;6];
pt4_R = [0;-9;4];
pt5_R = [0;-6;4];
pt6_R = [0;-9;2];

plot3(pt1_R(1),pt1_R(2),pt1_R(3),"b.","MarkerSize",10); hold on;
plot3(pt2_R(1),pt2_R(2),pt2_R(3),"b.","MarkerSize",10); hold on;
plot3(pt3_R(1),pt3_R(2),pt3_R(3),"b.","MarkerSize",10); hold on;
plot3(pt4_R(1),pt4_R(2),pt4_R(3),"b.","MarkerSize",10); hold on;
plot3(pt5_R(1),pt5_R(2),pt5_R(3),"b.","MarkerSize",10); hold on;
plot3(pt6_R(1),pt6_R(2),pt6_R(3),"b.","MarkerSize",10); hold on;

pt1_S = [0;-11;2];
pt2_S = [0;-14;2];
pt3_S = [0;-14;4];
pt4_S = [0;-11;4];
pt5_S = [0;-11;6];
pt6_S = [0;-14;6];

plot3(pt1_S(1),pt1_S(2),pt1_S(3),"b.","MarkerSize",10); hold on;
plot3(pt2_S(1),pt2_S(2),pt2_S(3),"b.","MarkerSize",10); hold on;
plot3(pt3_S(1),pt3_S(2),pt3_S(3),"b.","MarkerSize",10); hold on;
plot3(pt4_S(1),pt4_S(2),pt4_S(3),"b.","MarkerSize",10); hold on;
plot3(pt5_S(1),pt5_S(2),pt5_S(3),"b.","MarkerSize",10); hold on;
plot3(pt6_S(1),pt6_S(2),pt6_S(3),"b.","MarkerSize",10); hold on;
