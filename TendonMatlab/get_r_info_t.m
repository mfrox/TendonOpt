function [r,r_dot,r_ddot]=get_r_info_t(s,N,C,D)
        N_t = N.N_t;
        N_a = N.N_a;
        N_m = N.N_m;
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