C_star = [params_star(1:3)';
params_star(4:6)';
params_star(7:9)';
params_star(10:12)'];
D_star = [params_star(13:15)';
params_star(16:18)';
params_star(19:21)';
params_star(22:24)'];
% L_star = params_star(25);

robot = fkin_opt_four_mex(zeros(1,4),L,C_star,D_star);
    s=robot.s;
    for j = 1:length(s)
        D1(j)=D_star(1,1)+D_star(1,2)*s(j)+D_star(1,3)*s(j)^2 
        D2(j)=D_star(2,1)+D_star(2,2)*s(j)+D_star(2,3)*s(j)^2
        D3(j)=D_star(3,1)+D_star(3,2)*s(j)+D_star(3,3)*s(j)^2 
        D4(j)=D_star(4,1)+D_star(4,2)*s(j)+D_star(4,3)*s(j)^2

    end
%%
robot = fkin_opt_four_mex(zeros(1,4),0.07,C,D);
    s=robot.s;
 for j = 1:length(s)
  [r(:,:,j),~,~]=get_r_info_t(s(j),N,C,D)
 end
%%
 for k=1:N.N_t

    pt=zeros(length(s),3);
%     length(s)
    for i=1:length(s)
        p = [robot.p(i,1);robot.p(i,2);robot.p(i,3)];
        R = robot.R(:,:,i);
        [r,~,~]=get_r_info_t(s(i),N,C,D)
        pt(i,1:3)= p+R*r(:,k)
    end

plot3(pt(:,1),pt(:,2),pt(:,3),'k','LineWidth',3)
hold on
axis equal
plot3(robot.p(:,1),robot.p(:,2),robot.p(:,3),'b','LineWidth',3)
 end

 %%
     for j = 1:length(s)
        d1(j) = sqrt(r(1,1,j)^2+r(2,1,j)^2) ;
       d2(j) =  sqrt(r(1,2,j)^2+r(2,2,j)^2);
       d3(j) =  sqrt(r(1,3,j)^2+r(2,3,j)^2);
       d4(j) =  sqrt(r(1,4,j)^2+r(2,4,j)^2);

    end