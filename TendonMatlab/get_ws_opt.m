function [R3d,X_log] = get_ws_opt(x1,x2,x3,x4,C,D)
% C=[-pi/2,0;
%    -4*pi/14+pi/2, 0; 
%    pi/2, 0;
%    -pi/2, -2*pi/0.2];
% D=[.0055;
%    .0055;
%    .0055;
%    .0055];
% robot=struct('p',zeros(86,3),'R',zeros(3,3,86));
% coder.extrinsic('fkin_opt_four_mex')
% robot.p=zeros(86,3)
% robot.R=zeros(3,3,86)
% robot.r=zeros(3,4,86)
% R3d = zeros(83521,3);
% X_log = zeros(83521,4);
% robot_log = [];
R3d=[];
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
%                         y = robot.p(end,2);
%                         z = robot.p(end,3);
        
                     % create matrix of x, y, z
                    
                     R3d = [R3d; p_tip];
                     X_log = [X_log; X]
%                      robot_log = [robot_log; robot];
                end
            end
        end    
     end
end



