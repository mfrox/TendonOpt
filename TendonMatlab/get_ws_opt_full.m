function [p3d,R_log,X_log] = get_ws_opt_full(x1,x2,x3,x4,C,D)

p3d=[];
R_log = [];
X_log = [];
L = 0.070;
    parfor ii = 1:length(x1)
        for jj = 1:length(x2)
            for kk = 1:length(x3)
                for mm = 1:length(x4)
                        X = [x1(ii) x2(jj) x3(kk) x4(mm)];
                    %     % Create psis and betas from grid
                    %     lims{ii} = grid(1:4, ii);
                    %     lims{ii};
                        
                        robot = fkin_opt_four_mex(X,L,C,D);
                        
                        % x, y, and z position of tip
                        p_tip = robot.p(end,:)
                        R_tip = robot.R(:,:,end)
%                         y = robot.p(end,2);
%                         z = robot.p(end,3);
        
                     % create matrix of x, y, z
                    
                     p3d = [p3d; p_tip];
                     R_log = [R_log; R_tip];

                     X_log = [X_log; X]
%                      robot_log = [robot_log; robot];
                end
            end
        end    
     end
end



