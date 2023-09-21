% Brain Optimization

C0=[0 2*pi 0;
pi/4 0 pi/2;
-pi/2 -pi/2 0];



options = optimset('Display','Iter');
[C,~] = fminsearch(@(C) err_cost(C), C0, options);
[e,R3d,X_log] = err_cost(C0);
function [e,R3d,X_log] = err_cost(C)
    load('pBrain_robot.mat','pBrain_robot');
    pBrain_robot = pBrain_robot/1000;
    x1 = 0:0.25:4;
    x2 = 0:0.25:4;
    x3 = 0:0.25:4;
    x4 = 0:0.25:4;
    D = [0.0035;
        0.0035;
        0.0035
];
    [R3d,X_log] = get_ws_opt(x1,x2,x3,C,D);

    l_R = length(R3d);
    l_p = length(pBrain_robot);
    dist = zeros(l_R,l_p);

    parfor i = 1:length(R3d)
        for j = 1:10
            dist(i,j) = norm(pBrain_robot(j,:)-R3d(i,:));
        end
    end
    errs = min(dist);

    e = sum(errs.^2);
end

