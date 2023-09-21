% LSQ opt
% 


%%
C = [0 0 0;
pi/2 0 0;
pi 0 0
3*pi/2 0 0];

D = [0.0035 0 0;
0.0035 0 0;
0.0035 0 0
0.0035 0 0];

q0 = [0 0 0];

% q_i =  3.*rand(10,4);

x1 = 0:0.25:1;
x2 = 0:0.25:1;
x3 = 0:0.25:1;
x4 = 0:0.25:1;
I = randi(625,[10,1]);

[R3d,X_log] = get_ws_opt(x1,x2,x3,x4,C,D);
% q0 = [2.3882    2.7276    2.7277    0.9473];

%%

for m = 1:10
    y_start(m,:) = R3d(I(m),:);
    q_start(m,:) = X_log(I(m),:);
    R_id(:,:,m) = eye(3,3);
end



    
%%
L = 0.07;
for j = 1:10
 robot = fkin_opt_four_mex(q_start(j,:),L,C,D);
 y_d(j,:) = robot.p(end,:);
 R_d(:,:,j) = robot.R(:,:,end);
 p_d{j} = robot.p;
%  Rd{j} = robot.R;
end


%%

params_0 = [zeros(16,1);0.070];
[fstar, xstar, grad, hessian, state] = asamin('minimize','err_cost',params_0,[-Inf*ones(12,1);zeros(4,1);0.05],...
    [Inf*ones(12,1);0.1*ones(4,1);0.1],-ones(25,1));

%%
options = optimset('Display','Iter','UseParallel',true,'MaxIter',600,'TolFun',1e-2)
%  options = optimoptions('lsqnonlin','Display','Iter','StepTolerance',1e-4,'MaxFunctionEvaluations',...
% 1e4,'MaxIterations',1e4,'UseParallel',true);
%      options.Algorithm = 'levenberg-marquardt';

tic
[C_star,~] = fminsearch(@(C) err_cost(C,y_d,R_d), C0, options);
% problem = createOptimProblem('lsqnonlin','objective',@(C) err_cost(C,y_d,R_d),'x0',C0,'lb',[],'ub',[],'options',options);
% ms = MultiStart;
% [C_star,~] = run(ms,problem,20);
% [C_star,~] = lsqnonlin(@(C) err_cost(C,y_d,R_d), C0,[],[], options);
toc
% [E] = err_cost(C);
%%
% C0 = zeros(4,3);
% params_0 = [zeros(18,1)];
params_0 = [zeros(16,1);0.07;0.07;0.07;0.070];
options = optimoptions(@simulannealbnd,'Display','Iter','ObjectiveLimit',1e-3,'MaxIter',2400, ...
    'PlotFcn',{'saplotf','saplotbestx','saplotstopping','saplottemperature'});
 options.TemperatureFcn = @asatemp;
 options.AnnealingFcn = @asaAnnealingDist;
options.ReannealInterval = 20;
options.InitialTemperature = 10;

% [C,~] = fminsearch(@(C) err_cost(C), C0, options);
[params_star,~] = simulannealbnd(@(params) err_cost(params,y_d,R_id),params_0,[-200*ones(16,1);0.01*ones(4,1)],...
    [200*ones(16,1);0.07*ones(4,1)],options);

%%
    for i = 1:10
    [E(:,:,i),e(:,:,i),W(:,:,i)] = findPose(q_star(i,:),C,D,y_d(i,:),R_id);
    end
% %         q0 = [0 0 0 0];
%         options = optimoptions('lsqnonlin','Display','Iter','FunctionTolerance',1e-8,...
% 'StepTolerance',1e-4,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);
%         options.Algorithm = 'levenberg-marquardt';
% % options = optimset('Display','Iter','ObjectiveLimit',1e-6)
%     for i = 1:10
%         q_star(i,:) = lsqnonlin(@(q) findPose(q,C,D,y_d(i,:),R_d(:,:,i)),q0,[0 0 0 0],[4 4 4 4],options);
% %         q_star = fminsearch(@(q) findPose(q,C,D),q0,options)
%         [E(:,:,i),e(:,:,i),W(:,:,i)] = findPose(q_star(i,:),C,D,y_d(i,:),R_d(:,:,i))    ;
%     end
    %% Just runs pose opt
      options = optimoptions('lsqnonlin','Display','Iter','FunctionTolerance',1e-4,...
          'StepTolerance',1e-4,'MaxFunctionEvaluations',1e4,'MaxIterations',1e4);
     options.Algorithm = 'levenberg-marquardt';

    C = [params_star(1:4)'; 
       params_star(5:8)';
       params_star(9:12)';
       params_star(13:16)'];
%     D = [params(10:12)'; 
%        params(13:15)';
%        params(16:18)'];
D = [0.0035; 0.0035; 0.0035; 0.0035];
s_end = zeros(4,1);
    s_end = round(params_star(17:20)'*1000)/1000;
    L = 0.07;
%%
   % s_end = round(s_end*1000)/1000;
%    L = params(25);
%     L = params_star(16);
%     s_end(4) = L;
     tic
    parfor i = 1:10
        problem = createOptimProblem('lsqnonlin','objective',@(q) findPose(q,C,D,L,s_end,y_d(i,:),R_id(:,:,i)),...
            'x0',[0 0 0 0 0],'lb',[0 0 0 0 -0.03],'ub',[6 6 6 6 0.03],'options',options);
        ms = MultiStart;
        [q_star(i,:),~] = run(ms,problem,30);
    end
    toc
%%
    for i = 1:10
    [E(:,i),e(:,i),W(:,:,i)] = findPose(q_star(i,:),C,D,L,s_end,disp,y_d(i,:),R_d(:,:,i));
    end

    %% ASA Min
    C0 = zeros(12,1);
    [fstar,C_star,~,~,~] = asamin('minimize','err_cost',C0,-100*ones(12,1),100*ones(12,1),-1*ones(12,1));
%% get optimized tube positions
for ii = 1:10
    robot = fkin_terminating_mex(q_star(ii,:),L,C,D,s_end);
     y_good(ii,:) = robot.p(end,:);
     R_good(:,:,ii) = robot.R(:,:,end);
     p_good{ii} = robot.p;
end
%% Get backbone position for unoptimized tubes
for jj = 1:10
    robot = fkin_opt_four_mex(q_start(jj,:),L,C,D);
     y_bad(jj,:) = robot.p(end,:);
     R_bad(:,:,jj) = robot.R(:,:,end);
     p_bad{jj} = robot.p;
end
%%
figure
 for k = 1:10
     draw_coordinate_system2(0.003,R_id(:,:,k),y_start(k,:),'rrr')    
     draw_coordinate_system2(0.003,R_bad(:,:,k),y_start(k,:),'kkk')  
     draw_coordinate_system2(0.003,R_good(:,:,k),y_good(k,:),'bbb')
%      draw_coordinate_system2(0.003,R_x(:,:,k),y_x(k,:),'bbb')
     hold on; grid on; daspect([1 1 1]);
     plot3(p_bad{k}(:,1), p_bad{k}(:,2), p_bad{k}(:,3),'k','LineWidth',3)
%      plot3(pd{k}(:,1), pd{k}(:,2), pd{k}(:,3),'r','LineWidth',3)
     plot3(p_good{k}(:,1), p_good{k}(:,2), p_good{k}(:,3),'b','LineWidth',3)
     drawnow
 end
 
 %%
