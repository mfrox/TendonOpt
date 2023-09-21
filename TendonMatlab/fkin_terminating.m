function robot = fkin_terminating(tau,L,C,D,s_end) %#codegen
  
% Rod parameters 
ro=.686/2/1000;
ri=0.533/2/1000;
I=1/4*pi*(ro^4-ri^4);
Ar=pi*(ro^2-ri^2);
J=2*I;
E=50*10^9;
nu=.33;
Gmod=E/(2*(1+nu));
EI = E*I;
JG = J*Gmod;
disp = tau(5);
T_disp = [eye(3,3) [0 0 disp]'; 0 0 0 1];

tau = tau(1:4);

N_t=size(C,1); %number of tendons
N_a=size(C,2); %number of angle coeffs per tendon
N_m=size(D,2); %number of magnitude coeffs per tendon

% Sort tendon termination, C vals and tension into ascending order
[S_sorted, I] = sort(s_end);
C_old = C;
tau_old = tau;

for k = 1:N_t
    C(k,:) = C_old(I(k),:);
    tau(k) = tau_old(I(k));
end

% Determine number of tendon termination points
N_term = 1;
diff = zeros(3,4);
S_new = S_sorted;

for m = 2:length(S_sorted)
    diff(m-1,m:end) = abs(S_sorted(m-1)-S_sorted(m:end));
    if diff(m-1,m) < 0.002
        N_term = N_term + 1;
%         S_new(m) = max(S_sorted(n:length(diff{m})))
    end
end

if N_term == 4
    for p = 1:N_t
        S_new(p) = max(S_sorted(1:N_t));
    end
else
    for p = 2:length(S_sorted)
        if diff(1,p)<0.002 
            S_new(1:p) = max(S_sorted(1:p));
        end
    end
    for p = 3:length(S_sorted)
        if diff(2,p)<0.002
            S_new(2:p) = max(S_sorted(2:p));
        end
    end
    for p = 4:length(S_sorted)
        if diff(3,p)<0.002
            S_new(3:p) = max(S_sorted(3:p));
        end
    end
end

S_sorted = S_new;
p_length = 0:0.001:max(S_sorted);
L1_length = 0:0.001:S_sorted(1);
L2_length = S_sorted(1):0.001:S_sorted(2);
L3_length = S_sorted(2):0.001:S_sorted(3);
L4_length = S_sorted(3):0.001:S_sorted(4);
p_n = zeros(length(p_length),3);
R_n = zeros(3,3,length(p_length));
L1 = zeros(length(L1_length),1);
L2 = zeros(length(L2_length)-1,1);
L3 = zeros(length(L3_length)-1,1);
L4 = zeros(length(L4_length)-1,1);
s_out = zeros(length(p_length),1);

%Inverses of the stiffness matrices
K_bt=[EI,0,0; 0,EI,0; 0,0,JG]; 
K_bt_inv=[1/EI,0,0;0,1/EI,0;0,0,1/JG];
K_se=[(Gmod*Ar),0,0;0,(Gmod*Ar),0;0,0,(E*Ar)]; 
K_se_inv=[1/(Gmod*Ar),0,0;0,1/(Gmod*Ar),0;0,0,1/(E*Ar)];  

    function [r,r_dot,r_ddot]=get_r_info(s)
        r=zeros(3,N_t);
        r_dot=zeros(3,N_t);
        r_ddot=zeros(3,N_t);
        S_a=zeros(N_a,1);
        S_ad=zeros(N_a,1);
        S_add=zeros(N_a,1);
        S_m=zeros(N_m,1);
        S_md=zeros(N_m,1);
        S_mdd=zeros(N_m,1);
        for j=1:N_a
            S_a(j,1)=s^(j-1);
            S_ad(j,1)=(j-1)*s^(max([0,(j-2)]));
            S_add(j,1)=(j-1)*(j-2)*s^(max([0,(j-3)]));
        end
        for j=1:N_m
            S_m(j,1)=s^(j-1);
            S_md(j,1)=(j-1)*s^(max([0,(j-2)]));
            S_mdd(j,1)=(j-1)*(j-2)*s^(max([0,(j-3)]));
        end
        for j=1:N_t
        r(:,j)=D(j,:)*S_m*[sin(C(j,:)*S_a);cos(C(j,:)*S_a);0];
        r_dot(:,j)=D(j,:)*S_md*[sin(C(j,:)*S_a);cos(C(j,:)*S_a);0]+...
                  +D(j,:)*S_m*[cos(C(j,:)*S_a)*C(j,:)*S_ad;-sin(C(j,:)*S_a)*C(j,:)*S_ad;0];
        r_ddot(:,j)=D(j,:)*S_mdd*[sin(C(j,:)*S_a);cos(C(j,:)*S_a);0]+...
                  +2*D(j,:)*S_md*[cos(C(j,:)*S_a)*C(j,:)*S_ad;-sin(C(j,:)*S_a)*C(j,:)*S_ad;0]+...
                  +D(j,:)*S_m*[cos(C(j,:)*S_a)*C(j,:)*S_add;-sin(C(j,:)*S_a)*C(j,:)*S_add;0]+...
                  +D(j,:)*S_m*[-sin(C(j,:)*S_a)*(C(j,:)*S_ad)^2;-cos(C(j,:)*S_a)*(C(j,:)*S_ad)^2;0];
        end
    end


    function u_hat=hat(u)
        u_hat=[0,-u(3),u(2);u(3),0,-u(1);-u(2),u(1),0];
    end

    function y_dot = deriv(s, y)
        
        %Turn ODE input iN_to N_amed variables
        %p = [y(1);y(2);y(3)];
        R = reshape(y(4:12),3,3);
        v = [y(13);y(14);y(15)];
        u = [y(16);y(17);y(18)];
        int_n = y(24);

        vhat=hat(v);
        uhat=hat(u);
        [r,r_dot,r_ddot]=get_r_info(s);
        
        A=zeros(3,3);B=zeros(3,3);G=zeros(3,3);H=zeros(3,3);
        a=zeros(3,1);b=zeros(3,1);
        si_dot=zeros(N_t,1);
        for j=1:N_t-int_n
            
            rhat=hat(r(:,j+int_n));
            pi_dot_b=uhat*r(:,j+int_n)+r_dot(:,j+int_n)+v;
            si_dot(j+int_n)=norm(pi_dot_b);
        
            Ai=-tau(j+int_n)*hat(pi_dot_b)^2/norm(pi_dot_b)^3;
            Bi=rhat*Ai;
            Gi=-Ai*rhat;
            Hi=-Bi*rhat;
        
            ai=Ai*(uhat*pi_dot_b+uhat*r_dot(:,j+int_n)+r_ddot(:,j+int_n));
            bi=rhat*ai;
            
            A=A+Ai;
            B=B+Bi;
            G=G+Gi;
            H=H+Hi;
            a=a+ai;
            b=b+bi;
        end
         
        c=-uhat*K_bt*(u)-vhat*K_se*(v-[0;0;1])-b;
        d=-uhat*K_se*(v-[0;0;1])-a;
        
        %Compute state variable derivatives 
        p_dot = R*v;
        R_dot = R*uhat; 
        xi_dot= [K_se+A, G; B, K_bt+H]\[d;c];
        v_dot=xi_dot(1:3);
        u_dot=xi_dot(4:6);
        int_dot = 0;
        y_dot = [p_dot;reshape(R_dot,9,1);v_dot;u_dot;norm(v);si_dot;int_dot];
    end
 
%Solve the Problem
% segment 1

if N_term == 4 % Solves problem for all 4 tendons terminating at same point
    [r,r_dot,~]=get_r_info(0); 
    uinit=[0;0;0];
    vinit=[0;0;1];
    for i=1:30
        uinit_hat=hat(uinit);  
        ninit=zeros(3,1);
        minit=zeros(3,1);
        for k=1:N_t
            pi_dot_init=(uinit_hat*r(:,k)+r_dot(:,k)+vinit);
            ninit=ninit-tau(k)*pi_dot_init/norm(pi_dot_init);
            minit=minit-tau(k)*hat(r(:,k))*pi_dot_init/norm(pi_dot_init);
        end
        vnew=K_se_inv*ninit+[0;0;1];
        unew=K_bt_inv*minit;
        %pi_dot_0=(hat(K_bt_inv*(-tau*hat(r)*pi_dot_0/norm(pi_dot_0)))*r+r_dot+K_se_inv*(-tau*pi_dot_0/norm(pi_dot_0))+[0;0;1]);
        if norm(vnew-vinit)/norm(vinit)<1e-5 && norm(unew-uinit)/norm(uinit)<1e-6
            break
        else
            vinit=vnew;
            uinit=unew;
        end
    end
    correct_guess=[vnew;unew];
    
    p_init=[0;0;0];  
    R_init=reshape(eye(3),9,1);
    v_init=correct_guess(1:3); 
    u_init=correct_guess(4:6);
    int_n = 0;
    y_init=[p_init;R_init;v_init;u_init;0;zeros(N_t,1);int_n];        
    [s,y]=ode45(@deriv,0:.001:L,y_init);
    p1 = zeros(size(y,1),3);
    R1 = zeros(3,3,size(y,1));
    v1 = zeros(size(y,1),3);
    u1 =  zeros(size(y,1),3);
    s1 = zeros(size(y,1),1);
    L_1 =  zeros(size(y,1),4);

    
     for i = 1:size(y,1)
        p1(i,1:3) = y(i,1:3);
        R1(:,:,i) = reshape(y(i,4:12),3,3);
        v1(i,1:3) = y(i,13:15);
        u1(i,1:3) = y(i,16:18);
        s1(i,:) = y(i,19);
        L_1(i,1:4) = y(i,20:23);
     end
    p_n = zeros(size(p1,1),3);
    R_n = zeros(3,3,size(p1,1));
    L1 =  zeros(size(y,1),1);
    L2 =  zeros(size(y,1),1);
    L3 =  zeros(size(y,1),1);
    L4 =  zeros(size(y,1),1);

    p_n = p1;

    R_n = R1;
    
    L1 = L_1(:,1);
    L2 = L_1(:,2);
    L3 = L_1(:,3);
    L4 = L_1(:,4);
    s_out = zeros(size(p_n,1),1);
    s_out = s1;

elseif N_term == 3 % Solves problem for two termination points

    % segment 1
    [r,r_dot,~]=get_r_info(0); 
    uinit=[0;0;0];
    vinit=[0;0;1];
    for i=1:30
        uinit_hat=hat(uinit);  
        ninit=zeros(3,1);
        minit=zeros(3,1);
        for k=1:N_t
            pi_dot_init=(uinit_hat*r(:,k)+r_dot(:,k)+vinit);
            ninit=ninit-tau(k)*pi_dot_init/norm(pi_dot_init);
            minit=minit-tau(k)*hat(r(:,k))*pi_dot_init/norm(pi_dot_init);
        end
        vnew=K_se_inv*ninit+[0;0;1];
        unew=K_bt_inv*minit;
        if norm(vnew-vinit)/norm(vinit)<1e-5 && norm(unew-uinit)/norm(uinit)<1e-6
            break
        else
            vinit=vnew;
            uinit=unew;
        end
    end
    correct_guess=[vnew;unew];
    
    p_init=[0;0;0];  
    R_init=reshape(eye(3),9,1);
    v_init=correct_guess(1:3); 
    u_init=correct_guess(4:6);
    int_n = 0;
    y_init=[p_init;R_init;v_init;u_init;0;zeros(N_t,1);int_n];        
    [s_1,y_1]=ode45(@deriv,0:.001:min(S_sorted),y_init);
    p_1 = zeros(size(y_1,1),3);
    R_1 = zeros(3,3,size(y_1,1));
    v_1 = zeros(size(y_1,1),3);
    u_1 =  zeros(size(y_1,1),3);
    s_1 = zeros(size(y_1,1),1);
    L_1 =  zeros(size(y_1,1),4);
    L1 = zeros(size(y_1,1),1);
    
     for i = 1:size(y_1,1)
        p_1(i,1:3) = y_1(i,1:3);
        R_1(:,:,i) = reshape(y_1(i,4:12),3,3);
        v_1(i,1:3) = y_1(i,13:15);
        u_1(i,1:3) = y_1(i,16:18);
        s_1(i,:) = y_1(i,19);
        L_1(i,1:4) = y_1(i,20:23);
    end
     
    % segment 2
    [r,r_dot,~]=get_r_info(min(S_sorted)); 
    uinit2=[0;0;0];
    vinit2=[0;0;1];
    for i=1:30
        uinit2_hat=hat(uinit2);  
        ninit2=zeros(3,1);
        minit2=zeros(3,1);
        if S_sorted(1) == S_sorted(2) && S_sorted(3) == S_sorted(4)
            for k=3:N_t
                pi_dot_init2=(uinit2_hat*r(:,k)+r_dot(:,k)+vinit2);
                ninit2=ninit2-tau(k)*pi_dot_init2/norm(pi_dot_init2);
                minit2=minit2-tau(k)*hat(r(:,k))*pi_dot_init2/norm(pi_dot_init2);
                int_n = 2;
            end
        elseif S_sorted(1) ~= S_sorted(2) 
            for k=2:N_t
                pi_dot_init2=(uinit2_hat*r(:,k)+r_dot(:,k)+vinit2);
                ninit2=ninit2-tau(k)*pi_dot_init2/norm(pi_dot_init2);
                minit2=minit2-tau(k)*hat(r(:,k))*pi_dot_init2/norm(pi_dot_init2);
                int_n = 1;
            end
        elseif S_sorted(3) ~= S_sorted(4)
            for k=4:N_t
                pi_dot_init2=(uinit2_hat*r(:,k)+r_dot(:,k)+vinit2);
                ninit2=ninit2-tau(k)*pi_dot_init2/norm(pi_dot_init2);
                minit2=minit2-tau(k)*hat(r(:,k))*pi_dot_init2/norm(pi_dot_init2);
                int_n = 3;
            end
        end
        v2new=K_se_inv*ninit2+[0;0;1];
        u2new=K_bt_inv*minit2;
        if norm(v2new-vinit2)/norm(vinit2)<1e-5 && norm(u2new-uinit2)/norm(uinit2)<1e-6
            break
        else
            vinit2=v2new;
            uinit2=u2new;
        end
    end
    correct_guess2=[v2new;u2new];
        
    p_init2=[p_1(end,1);p_1(end,2);p_1(end,3)];  
    R_init2=reshape(R_1(:,:,end),9,1);
    v_init2=correct_guess2(1:3); 
    u_init2=correct_guess2(4:6);
    
    y_init2=[p_init2;R_init2;v_init2;u_init2;S_sorted(1);L_1(end,1:4)';int_n];        
    [s_2,y_2]=ode45(@deriv,min(S_sorted):.001:max(S_sorted),y_init2);
    
    p_2 = zeros(size(y_2,1),3);
    R_2 = zeros(3,3,size(y_2,1));
    v_2 = zeros(size(y_2,1),3);
    u_2 =  zeros(size(y_2,1),3);
    s_2 = zeros(size(y_2,1),1);
    L_2 =  zeros(size(y_2,1),4);

    
    for i = 1:size(y_2,1)
        p_2(i,1:3) = y_2(i,1:3);
        R_2(:,:,i) = reshape(y_2(i,4:12),3,3);
        v_2(i,1:3) = y_2(i,13:15);
        u_2(i,1:3) = y_2(i,16:18);
        s_2(i,:) = y_2(i,19);
        L_2(i,1:4) = y_2(i,20:23);
    end

    p_n = zeros(size(p_1,1)+size(p_2,1)-1,3);
    R_n = zeros(3,3,size(p_1,1)+size(p_2,1)-1);
    L1 = zeros(size(y_1,1));
    L2 = zeros(size(y_2,1),1);
    L3 = zeros(size(y_2,1),1);
    L4 = zeros(size(y_2,1),1);

    p_n = [p_1;p_2(2:end,:)];
    l1 = length(p_1);
    l2 = length(p_1)+length(p_2)-1;
    
    R_n = zeros(3,3,l2);
    R_n(:,:,1:l1) = R_1;
    R_n(:,:,l1+1:l2) = R_2(:,:,2:end);
    
    L1 = L_1(:,1);
    L2 = [L_1(:,2);L_2(2:end,2)];
    L3 = [L_1(:,3);L_2(2:end,3)];
    L4 = [L_1(:,4);L_2(2:end,4)];
    s_out = zeros(size(p_n,1),1);
    s_out = [s_1;s_2(2:end)];
    
elseif N_term == 2 % Solves problem for 2 tendons terminating at same point
    
    % segment 1
    [r,r_dot,~]=get_r_info(S_sorted(1)); 
    uinit=[0;0;0];
    vinit=[0;0;1];
    for i=1:30
        uinit_hat=hat(uinit);  
        ninit=zeros(3,1);
        minit=zeros(3,1);
        for k=1:N_t
            pi_dot_init=(uinit_hat*r(:,k)+r_dot(:,k)+vinit);
            ninit=ninit-tau(k)*pi_dot_init/norm(pi_dot_init);
            minit=minit-tau(k)*hat(r(:,k))*pi_dot_init/norm(pi_dot_init);
        end
        vnew=K_se_inv*ninit+[0;0;1];
        unew=K_bt_inv*minit;
        if norm(vnew-vinit)/norm(vinit)<1e-5 && norm(unew-uinit)/norm(uinit)<1e-6
            break
        else
            vinit=vnew;
            uinit=unew;
        end
    end
    correct_guess=[vnew;unew];
        
    p_init=[0;0;0];  
    R_init=reshape(eye(3),9,1);
    v_init=correct_guess(1:3); 
    u_init=correct_guess(4:6);
    int_n = 0;
    y_init=[p_init;R_init;v_init;u_init;0;zeros(N_t,1);int_n];        
    [s_1,y_1]=ode45(@deriv,0:.001:S_sorted(1),y_init);
    p_1 = zeros(size(y_1,1),3);
    R_1 = zeros(3,3,size(y_1,1));
    v_1 = zeros(size(y_1,1),3);
    u_1 =  zeros(size(y_1,1),3);
    s_1 = zeros(size(y_1,1),1);
    L_1 =  zeros(size(y_1,1),4);
    L1 = zeros(size(y_1,1),1);
    
     for i = 1:size(y_1,1)
        p_1(i,1:3) = y_1(i,1:3);
        R_1(:,:,i) = reshape(y_1(i,4:12),3,3);
        v_1(i,1:3) = y_1(i,13:15);
        u_1(i,1:3) = y_1(i,16:18);
        s_1(i,:) = y_1(i,19);
        L_1(i,1:4) = y_1(i,20:23);
        L1(i,:) = L_1(i,1);
    end
     
    % segment 2
    [r,r_dot,~]=get_r_info(S_sorted(2)); 
    uinit2=[0;0;0];
    vinit2=[0;0;1];
    for i=1:30
        uinit2_hat=hat(uinit2);  
        ninit2=zeros(3,1);
        minit2=zeros(3,1);
        if S_sorted(1) == S_sorted(2)
            for k=3:N_t
                pi_dot_init2=(uinit2_hat*r(:,k)+r_dot(:,k)+vinit2);
                ninit2=ninit2-tau(k)*pi_dot_init2/norm(pi_dot_init2);
                minit2=minit2-tau(k)*hat(r(:,k))*pi_dot_init2/norm(pi_dot_init2);
                int_n = 2;
            end
        else
             for k=2:N_t
                pi_dot_init2=(uinit2_hat*r(:,k)+r_dot(:,k)+vinit2);
                ninit2=ninit2-tau(k)*pi_dot_init2/norm(pi_dot_init2);
                minit2=minit2-tau(k)*hat(r(:,k))*pi_dot_init2/norm(pi_dot_init2);
                int_n = 1;
             end
        end
        v2new=K_se_inv*ninit2+[0;0;1];
        u2new=K_bt_inv*minit2;
        if norm(v2new-vinit2)/norm(vinit2)<1e-5 && norm(u2new-uinit2)/norm(uinit2)<1e-6
            break
        else
            vinit2=v2new;
            uinit2=u2new;
        end
    end
    correct_guess2=[v2new;u2new];
        
    p_init2=[p_1(end,1);p_1(end,2);p_1(end,3)];  
    R_init2=reshape(R_1(:,:,end),9,1);
    v_init2=correct_guess2(1:3); 
    u_init2=correct_guess2(4:6);
    y_init2=[p_init2;R_init2;v_init2;u_init2;S_sorted(1);L_1(end,1:4)';int_n];  
    uA = unique(S_sorted);
    [s_2,y_2]=ode45(@deriv,min(S_sorted):.001:uA(2),y_init2);
    
    p_2 = zeros(size(y_2,1),3);
    R_2 = zeros(3,3,size(y_2,1));
    v_2 = zeros(size(y_2,1),3);
    u_2 =  zeros(size(y_2,1),3);
    s_2 = zeros(size(y_2,1),1);
    L_2 =  zeros(size(y_2,1),4);
    
    for i = 1:size(y_2,1)
        p_2(i,1:3) = y_2(i,1:3);
        R_2(:,:,i) = reshape(y_2(i,4:12),3,3);
        v_2(i,1:3) = y_2(i,13:15);
        u_2(i,1:3) = y_2(i,16:18);
        s_2(i,:) = y_2(i,19);
        L_2(i,1:4) = y_2(i,20:23);
    end
    
    % segment 3
    [r,r_dot,~]=get_r_info(S_sorted(2)); 
    uinit3=[0;0;0];
    vinit3=[0;0;1];
    for i=1:30
        uinit3_hat=hat(uinit3);  
        ninit3=zeros(3,1);
        minit3=zeros(3,1);
        if S_sorted(1) == S_sorted(2) || S_sorted(2) == S_sorted(3)
            for k=4:N_t
                pi_dot_init3=(uinit3_hat*r(:,k)+r_dot(:,k)+vinit3);
                ninit3=ninit3-tau(k)*pi_dot_init3/norm(pi_dot_init3);
                minit3=minit3-tau(k)*hat(r(:,k))*pi_dot_init3/norm(pi_dot_init3);
                int_n = 3; 
            end
        else
            for k=3:N_t
                pi_dot_init3=(uinit3_hat*r(:,k)+r_dot(:,k)+vinit3);
                ninit3=ninit3-tau(k)*pi_dot_init3/norm(pi_dot_init3);
                minit3=minit3-tau(k)*hat(r(:,k))*pi_dot_init3/norm(pi_dot_init3);
                int_n = 2;
            end
        end

        v3new=K_se_inv*ninit3+[0;0;1];
        u3new=K_bt_inv*minit3;
        if norm(v3new-vinit3)/norm(vinit3)<1e-5 && norm(u3new-uinit3)/norm(uinit3)<1e-6
            break
        else
            vinit3=v3new;
            uinit3=u3new;
        end
    end
    correct_guess3=[v3new;u3new];
        
    p_init3=[p_2(end,1);p_2(end,2);p_2(end,3)];  
    R_init3=reshape(R_2(:,:,end),9,1);
    v_init3=correct_guess3(1:3); 
    u_init3=correct_guess3(4:6);
    y_init3=[p_init3;R_init3;v_init3;u_init3;S_sorted(2);L_2(end,1:4)';int_n];        
    [s_3,y_3]=ode45(@deriv,uA(2):.001:max(S_sorted),y_init3);
    
    p_3 = zeros(size(y_3,1),3);
    R_3 = zeros(3,3,size(y_3,1));
    v_3 = zeros(size(y_3,1),3);
    u_3 =  zeros(size(y_3,1),3);
    s_3 = zeros(size(y_3,1),1);
    L_3 =  zeros(size(y_3,1),4);
    
    for i = 1:size(y_3,1)
        p_3(i,1:3) = y_3(i,1:3);
        R_3(:,:,i) = reshape(y_3(i,4:12),3,3);
        v_3(i,1:3) = y_3(i,13:15);
        u_3(i,1:3) = y_3(i,16:18);
        s_3(i,:) = y_3(i,19);
        L_3(i,1:4) = y_3(i,20:23);
    end
    p_n = zeros(size(p_1,1)+size(p_2,1)-1+size(p_3,1)-1,3);
    R_n = zeros(3,3,size(p_1,1)+size(p_2,1)-1+size(p_3,1)-1);
    L1 = zeros(size(y_1,1));
    L2 = zeros(size(y_2,1),1);
    L3 = zeros(size(y_3,1),1);
    L4 = zeros(size(y_3,1),1);

    p_n = [p_1;p_2(2:end,:);p_3(2:end,:)];
    l1 = length(p_1);
    l2 = length(p_1)+length(p_2)-1;
    l3 = length(p_1)+length(p_2)+length(p_3)-2;
    
    R_n = zeros(3,3,l3);
    R_n(:,:,1:l1) = R_1;
    R_n(:,:,l1+1:l2) = R_2(:,:,2:end);
    R_n(:,:,l2+1:l3) = R_3(:,:,2:end);
    
    L1 = L_1(:,1);
    if S_sorted(1) == S_sorted(2)
        L2 = [L_1(:,2)];
        L3 = [L_1(:,3);L_2(2:end,3)];
    else
        L2 = [L_1(:,2);L_2(2:end,2)];
        L3 = [L_1(:,3);L_2(2:end,3);L_3(2:end,3)];
    end
    
    L4 = [L_1(:,4);L_2(2:end,4);L_3(2:end,4)];
    s_out = zeros(size(p_n,1),1);
    s_out = [s_1;s_2(2:end);s_3(2:end)];
    
elseif N_term == 1 % 4 different termination points
    % segment 1
    [r,r_dot,~]=get_r_info(0); 
    uinit=[0;0;0];
    vinit=[0;0;1];
    for i=1:30
        uinit_hat=hat(uinit);  
        ninit=zeros(3,1);
        minit=zeros(3,1);
        for k=1:N_t
            pi_dot_init=(uinit_hat*r(:,k)+r_dot(:,k)+vinit);
            ninit=ninit-tau(k)*pi_dot_init/norm(pi_dot_init);
            minit=minit-tau(k)*hat(r(:,k))*pi_dot_init/norm(pi_dot_init);
        end
        vnew=K_se_inv*ninit+[0;0;1];
        unew=K_bt_inv*minit;
        if norm(vnew-vinit)/norm(vinit)<1e-5 && norm(unew-uinit)/norm(uinit)<1e-6
            break
        else
            vinit=vnew;
            uinit=unew;
        end
    end
    correct_guess=[vnew;unew];
        
    p_init=[0;0;0];  
    R_init=reshape(eye(3),9,1);
    v_init=correct_guess(1:3); 
    u_init=correct_guess(4:6);
    int_n = 0;
    y_init=[p_init;R_init;v_init;u_init;0;zeros(N_t,1);int_n];        
    [s_1,y_1]=ode45(@deriv,0:.001:S_sorted(1),y_init);
    p_1 = zeros(size(y_1,1),3);
    R_1 = zeros(3,3,size(y_1,1));
    v_1 = zeros(size(y_1,1),3);
    u_1 =  zeros(size(y_1,1),3);
    s_1 = zeros(size(y_1,1),1);
    L_1 =  zeros(size(y_1,1),4);

     for i = 1:size(y_1,1)
        p_1(i,1:3) = y_1(i,1:3);
        R_1(:,:,i) = reshape(y_1(i,4:12),3,3);
        v_1(i,1:3) = y_1(i,13:15);
        u_1(i,1:3) = y_1(i,16:18);
        s_1(i,:) = y_1(i,19);
        L_1(i,1:4) = y_1(i,20:23);
    end
     
    % segment 2
    [r,r_dot,~]=get_r_info(S_sorted(1)); 
    uinit2=[0;0;0];
    vinit2=[0;0;1];
    for i=1:30
        uinit2_hat=hat(uinit2);  
        ninit2=zeros(3,1);
        minit2=zeros(3,1);
        for k=2:N_t
            pi_dot_init2=(uinit2_hat*r(:,k)+r_dot(:,k)+vinit2);
            ninit2=ninit2-tau(k)*pi_dot_init2/norm(pi_dot_init2);
            minit2=minit2-tau(k)*hat(r(:,k))*pi_dot_init2/norm(pi_dot_init2);
        end
        v2new=K_se_inv*ninit2+[0;0;1];
        u2new=K_bt_inv*minit2;
        if norm(v2new-vinit2)/norm(vinit2)<1e-5 && norm(u2new-uinit2)/norm(uinit2)<1e-6
            break
        else
            vinit2=v2new;
            uinit2=u2new;
        end
    end
    correct_guess2=[v2new;u2new];
        
    p_init2=[p_1(end,1);p_1(end,2);p_1(end,3)];  
    R_init2=reshape(R_1(:,:,end),9,1);
    v_init2=correct_guess2(1:3); 
    u_init2=correct_guess2(4:6);
    int_n = 1; 
    y_init2=[p_init2;R_init2;v_init2;u_init2;S_sorted(1);L_1(end,1:4)';int_n];        
    [s_2,y_2]=ode45(@deriv,S_sorted(1):.001:S_sorted(2),y_init2);
    
    p_2 = zeros(size(y_2,1),3);
    R_2 = zeros(3,3,size(y_2,1));
    v_2 = zeros(size(y_2,1),3);
    u_2 =  zeros(size(y_2,1),3);
    s_2 = zeros(size(y_2,1),1);
    L_2 =  zeros(size(y_2,1),4);
    
    for i = 1:size(y_2,1)
        p_2(i,1:3) = y_2(i,1:3);
        R_2(:,:,i) = reshape(y_2(i,4:12),3,3);
        v_2(i,1:3) = y_2(i,13:15);
        u_2(i,1:3) = y_2(i,16:18);
        s_2(i,:) = y_2(i,19);
        L_2(i,1:4) = y_2(i,20:23);
    end
    
    % segment 3
    [r,r_dot,~]=get_r_info(S_sorted(2)); 
    uinit3=[0;0;0];
    vinit3=[0;0;1];
    for i=1:30
        uinit3_hat=hat(uinit3);  
        ninit3=zeros(3,1);
        minit3=zeros(3,1);
        for k=3:N_t
            pi_dot_init3=(uinit3_hat*r(:,k)+r_dot(:,k)+vinit3);
            ninit3=ninit3-tau(k)*pi_dot_init3/norm(pi_dot_init3);
            minit3=minit3-tau(k)*hat(r(:,k))*pi_dot_init3/norm(pi_dot_init3);
        end
        v3new=K_se_inv*ninit3+[0;0;1];
        u3new=K_bt_inv*minit3;
        if norm(v3new-vinit3)/norm(vinit3)<1e-5 && norm(u3new-uinit3)/norm(uinit3)<1e-6
            break
        else
            vinit3=v3new;
            uinit3=u3new;
        end
    end
    correct_guess3=[v3new;u3new];
        
    p_init3=[p_2(end,1);p_2(end,2);p_2(end,3)];  
    R_init3=reshape(R_2(:,:,end),9,1);
    v_init3=correct_guess3(1:3); 
    u_init3=correct_guess3(4:6);
    int_n = 2; 
    y_init3=[p_init3;R_init3;v_init3;u_init3;S_sorted(2);L_2(end,1:4)'; int_n];        
    [s_3,y_3]=ode45(@deriv,S_sorted(2):.001:S_sorted(3),y_init3);
    
    p_3 = zeros(size(y_3,1),3);
    R_3 = zeros(3,3,size(y_3,1));
    v_3 = zeros(size(y_3,1),3);
    u_3 =  zeros(size(y_3,1),3);
    s_3 = zeros(size(y_3,1),1);
    L_3 =  zeros(size(y_3,1),4);
    
    for i = 1:size(y_3,1)
        p_3(i,1:3) = y_3(i,1:3);
        R_3(:,:,i) = reshape(y_3(i,4:12),3,3);
        v_3(i,1:3) = y_3(i,13:15);
        u_3(i,1:3) = y_3(i,16:18);
        s_3(i,:) = y_3(i,19);
        L_3(i,1:4) = y_3(i,20:23);
    end
    
    % segment 4
    [r,r_dot,~]=get_r_info(S_sorted(3)); 
    uinit4=[0;0;0];
    vinit4=[0;0;1];
    for i=1:30
        uinit4_hat=hat(uinit4);  
        ninit4=zeros(3,1);
        minit4=zeros(3,1);
        for k=4:N_t
            pi_dot_init4=(uinit4_hat*r(:,k)+r_dot(:,k)+vinit4);
            ninit4=ninit4-tau(k)*pi_dot_init4/norm(pi_dot_init4);
            minit4=minit4-tau(k)*hat(r(:,k))*pi_dot_init4/norm(pi_dot_init4);
        end
        v4new=K_se_inv*ninit4+[0;0;1];
        u4new=K_bt_inv*minit4;
        if norm(v4new-vinit4)/norm(vinit4)<1e-5 && norm(u4new-uinit4)/norm(uinit4)<1e-6
            break
        else
            vinit4=v4new;
            uinit4=u4new;
        end
    end
    correct_guess4=[v4new;u4new];
        
    p_init4=[p_3(end,1);p_3(end,2);p_3(end,3)];  
    R_init4=reshape(R_3(:,:,end),9,1);
    v_init4=correct_guess4(1:3); 
    u_init4=correct_guess4(4:6);
    int_n = 3; 
    y_init4=[p_init4;R_init4;v_init4;u_init4;S_sorted(3);L_3(end,1:4)';int_n];        
    [s_4,y_4]=ode45(@deriv,S_sorted(3):.001:S_sorted(4),y_init4);
    
    p_4 = zeros(size(y_4,1),3);
    R_4 = zeros(3,3,size(y_4,1));
    v_4 = zeros(size(y_4,1),3);
    u_4 =  zeros(size(y_4,1),3);
    s_4 = zeros(size(y_4,1),1);
    L_4 =  zeros(size(y_4,1),4);
    
    for i = 1:size(y_4,1)
        p_4(i,1:3) = y_4(i,1:3);
        R_4(:,:,i) = reshape(y_4(i,4:12),3,3);
        v_4(i,1:3) = y_4(i,13:15);
        u_4(i,1:3) = y_4(i,16:18);
        s_4(i,:) = y_4(i,19);
        L_4(i,1) = y_4(i,20:23);
    end
    p_n = zeros(size(p_1,1)+size(p_2,1)-1+size(p_3,1)-1+size(p_4,1)-1,3);
    R_n = zeros(3,3,size(p_1,1)+size(p_2,1)-1+size(p_3,1)-1+size(p_4,1)-1);
    L1 = zeros(size(y_1,1));
    L2 = zeros(size(y_2,1),1);
    L3 = zeros(size(y_3,1),1);
    L4 = zeros(size(y_4,1),1);

    p_n = [p_1;p_2(2:end,:);p_3(2:end,:);p_4(2:end,:)];
    l1 = length(p_1);
    l2 = length(p_1)+length(p_2)-1;
    l3 = length(p_1)+length(p_2)+length(p_3)-2;
    l4 = length(p_1)+length(p_2)+length(p_3)+length(p_4)-3;
    
%     R_n = zeros(3,3,l4);
    R_n(:,:,1:l1) = R_1;
    R_n(:,:,l1+1:l2) = R_2(:,:,2:end);
    R_n(:,:,l2+1:l3) = R_3(:,:,2:end);
    R_n(:,:,l3+1:l4) = R_4(:,:,2:end);
    
    L1 = L_1(:,1);
    L2 = [L_1(:,2);L_2(2:end,2)];
    L3 = [L_1(:,3);L_2(2:end,3);L_3(2:end,3)];
    L4 = [L_1(:,4);L_2(2:end,4);L_3(2:end,4);L_4(2:end,4)];
    s_out = zeros(size(p_n,1),1);
    s_out = [s_1;s_2(2:end);s_3(2:end);s_4(2:end)];

end   

T_pn = zeros(length(p_n),4);
T_pn = cart2hom(p_n);

for i = 1:length(p_n)
    p_n(i,:) = hom2cart((T_disp*T_pn(i,:)')');
end


% L_i = [L_1;L_2;L_3;L_4];
%resulting tendon lengths
% L_i=y(end,end-N_t+1:end);
% 
% p_n = zeros(size(y,1),3);
% R_n = zeros(3,3,size(y,1));
r_n = zeros(3,N_t,length(s_out));
% for i = 1:size(y,1)
% p_n(i,1:3) = y(i,1:3);
% R_n(:,:,i) = reshape(y(i,4:12),3,3);
% 
% end
r_n = zeros(3,N_t,length(s_out));
for i=1:length(s_out)
    [r_n(:,:,i),~,~]=get_r_info(s_out(i));
end
  
robot.r = r_n;
robot.p = p_n;
robot.R = R_n;
% robot.s = s;
robot.L1 = L1;
robot.L2 = L2;
robot.L3 = L3;
robot.L4 = L4;

end