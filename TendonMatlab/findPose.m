     function [E,e,W] = findPose(q,C_star,D_star,L,s_end,y_d,R_d)
        robot = fkin_terminating_mex(q,L,C_star,D_star,s_end);
        y = robot.p(end,:);
        R = robot.R(:,:,end);          
%         R_diff = R*R_d';
        sigma_p = 0.002;
        sigma_R = pi/16;
        z_R = R(:,3);
        z_Rd = R_d(:,3);
        ang = atan2(norm(cross(z_R,z_Rd)),dot(z_R,z_Rd));
        W = diag([(1/sigma_p^2); (1/sigma_R^2)]);

        e_p = norm(y_d - y);
%         e_R = rotm2eul(R_diff);
        e_R = ang;
        e = [e_p e_R]'; % absolute error
        E = W*e;    % weighted error
     end