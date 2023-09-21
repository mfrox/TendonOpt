% clear all; close all;
% C = [0 0; pi/2 0; pi 0; 3*pi/2 0];
tau = [0 0 0 0 0];
% D = [0.0055; 0.0055; 0.0055; 0.0055];
L = 0.070;
L_f = 0.075;
% s_end = [0.015 0.075 0.024 0.075];
robot = fkin_terminating(tau,L,C,D,s_end);
L1 = robot.L1;
L2 = robot.L2;
L3 = robot.L3;
L4 = robot.L4;
p = robot.p;
R = robot.R;
r = robot.r;

%%

for i=1:length(L2)
    scatter3(r(1,2,i),r(2,2,i),s(i)); hold on;
end
%%
for i=1:length(L2)
    pt_2(i,1:3)= p(i,:)'+R(:,:,i)*r(:,2,i);
end
for i=1:length(L3)
    pt_3(i,1:3)= p(i,:)'+R(:,:,i)*r(:,3,i);
end
for i=1:length(L4)
    pt_4(i,1:3)= p(i,:)'+R(:,:,i)*r(:,4,i);
end
%%
for i=1:length(L1)
    pt_1(i,1:3)= p(i,:)'+R(:,:,i)*r(:,1,i);
end
for i=1:length(L2)
    pt_2(i,1:3)= p(i,:)'+R(:,:,i)*r(:,2,i);
%     draw_coordinate_system2(0.003,R(:,:,i),p(i,:),'bbb')
end
for i=1:length(L3)
    pt_3(i,1:3)= p(i,:)'+R(:,:,i)*r(:,3,i);
end
for i=1:length(L4)
    pt_4(i,1:3)= p(i,:)'+R(:,:,i)*r(:,4,i);
end
 plot3(pt_1(:,1),pt_1(:,2),pt_1(:,3),'k','LineWidth',3)
hold on;
plot3(pt_2(:,1),pt_2(:,2),pt_2(:,3),'r','LineWidth',3)

plot3(pt_3(:,1),pt_3(:,2),pt_3(:,3),'g','LineWidth',3)
plot3(pt_4(:,1),pt_4(:,2),pt_4(:,3),'c','LineWidth',3)
axis equal;

plot3(p(:,1),p(:,2),p(:,3),'b','LineWidth',3)