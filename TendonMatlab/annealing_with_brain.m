% Simulated Annealing with Brain Points

C = [0.0458893016045499	0.348060460028715	0.784411908136463
-0.345511536077701	-2.30910538769068	-1.16078835879915
-0.864092249752369	-0.322175462222833	-1.75168512657546
-2.67068864101874	-1.15699894252713	-0.513906592795182];

D = [0.0361074642672182	-0.470181722718556	0.206673606697001
0.0399745835825518	-0.932201215401590	1.17364457445820
0.0275147931444065	-0.440826912497184	0.179584382812298
0.0373811673399776	-0.0967567892437073	0.190516681957474];

tic
[e,I,p3d,R_log,X_log] = err_cost_brain(C,D);
toc

function [e,I,p3d,R_log,X_log] = err_cost_brain(C,D)
    load('pBrain_use','pBrain_use');
    x1 = 0:0.25:4.5;
    x2 = 0:0.25:4.5;
    x3 = 0:0.25:4.5;
    x4 = 0:0.25:4.5;
    [p3d,R_log,X_log] = get_ws_opt_full(x1,x2,x3,x4,C,D);

    l_R = length(p3d);
    l_p = length(pBrain_use);
    dist = zeros(l_R,l_p);

    parfor i = 1:length(p3d)
        for j = 1:l_p
            dist(i,j) = norm(pBrain_use(j,:)-p3d(i,:));
        end
    end
    [errs,I ]= min(dist);

    e = sum(errs.^2);
end