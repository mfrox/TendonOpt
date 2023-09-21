function [Err] = err_cost(params,y,R)

%      R_d = eye(3,3);
%      D = [0.0035;
%           0.0035;
%           0.0035;
%           0.0035];
    C = [params(1:4)'; 
       params(5:8)';
       params(9:12)';
       params(13:16)'];
%     D = [params(10:12)'; 
%        params(13:15)';
%        params(16:18)'];
D = [0.0035; 0.0035; 0.0035; 0.0035];
    s_end = zeros(1,4);
    s_end = params(17:20)';
    s_end = round(s_end*1000)/1000;
%     L = params(16);
    L = 0.07;
%      s_end(4) = L;
     options = optimoptions('lsqnonlin','Display','Off','FunctionTolerance',...
        1e-4,'StepTolerance',1e-4,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);
     options.Algorithm = 'levenberg-marquardt';
     tic
%      if s_end(3) >= s_end(4)-0.002
%          s_end(3) = s_end(4)-0.003;
%      end
%      if s_end(2) >= s_end(3)-0.002
%             s_end(2) = s_end(3)-0.003;
%      end
%      if s_end(1) >= s_end(2)-0.002
%          s_end(1) = s_end(2)-0.003;
%      end

%      s_end
% 
%      if abs(s_end(1)-s_end(2)) > 0.002 && abs(s_end(2)-s_end(3)) > 0.002 && ...
%              abs(s_end(1)-s_end(3)) > 0.002 && abs(s_end(3)-s_end(4)) > 0.002 ...
%              && abs(s_end(1)-s_end(4)) > 0.002 &&  abs(s_end(2)-s_end(4)) > 0.002
    for i = 1:10
        problem = createOptimProblem('lsqnonlin','objective',@(q) findPose(q,C,D,L,s_end,y(i,:),R(:,:,i)),...
            'x0',[0 0 0 0 0],'lb',[0 0 0 0 -0.03],'ub',[6 6 6 6 0.03],'options',options);
        ms = MultiStart;
        ms.UseParallel = true;
        ms.Display = 'off';
        [q_star(i,:),f(i)] =run(ms,problem,20);
        [E,e,W] = findPose(q_star(i,:),C,D,L,s_end,y(i,:),R(:,:,i));
        Err_met(i) = e'*W*e;
    end
     
    Err = sum(Err_met);
%      if s_end(1) > L || s_end(2) > L || s_end(3) > L || s_end(4) > L
%          Err = 1e6;
%      end
 
