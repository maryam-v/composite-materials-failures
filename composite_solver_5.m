%@ (M.A.H.S.A)
clc;
clear;
close all;
disp('z_mid-plane = 0 &  Assume the downward direction is positive.  ');
disp('Layers are numbered from top to bottom .  ');
disp('     ');
disp('     ');
disp('________________________                    z===> -');
disp('           k=1          ');
disp('________________________');
disp('           k=2          ');
disp('________________________');
disp('           k=3          ');
disp('________________________');
disp('            ...         ');
disp('________________________________ (mid plane) z=0');
disp('            ...         ');
disp('________________________');
disp('           k=n-2         ');
disp('________________________');
disp('           k=n-1          ');
disp('________________________');
disp('           k=n          ');
disp('________________________                      z===> + ');
disp('     ');
disp('     ');

%start program

A=0;
B=0;
D=0;
Tkol=0;
disp('     ');
disp('     ');
disp('----------------------start program----------------------');

n=input('inter number of ply   : ');
number_delete=input('inter number of deleted ply    : ');
number_ply_q_zero=zeros(1,number_delete);
if number_delete~=0
    for jj=1:number_delete
        number_ply_q_zero(1,jj)=input('shomare laye hazf shode    :  ');
    end
    
    
end
C=cell(19,n); % satre aval= Qbar .... satre dovom= T davaran .... satre sevom o chaharom = tanesh_12_up va tanesh_12_down   panjom tanesh_xy (shishom=tanesh_12= T*tanesh_xy) 7===>Q global   ; 8===>[]; 9===>[]   10va 11 tanesh_xy_up va tanesh_xy_down  12===>strain_xy    13==>strain_12  14=====>strain_xy_top  15====>strain_xy_down  16=====>strain_12_top  17====>strain_12_down
disp('--------------------------------------------');
t=zeros(n,2);  % soton aval zekhamat har laye soton dovo zaviye alyaf har laye


for i=1:n
    disp('--------------------------------------------------------------------------');
    disp('--------------------------------------------------------------------------');
    disp('For Layer:   ');
    disp(i);
    t(i,1)=input('zekhamate laye =   ');
    t(i,2)=input('The Fiber angle (Degree)=    ');
    
    disp('----------------------material----------------------');
    if i==1
        mm=input('               1===> mechanical properties        2===> local stiffness matrix         3===>Default propertie  ');
        
        
        
        if mm==1
            disp('--------------------------------------------');
            E1=input('inter E1=   ');
            E2=input('inter E2=   ');
            G12=input('inter G12=   ');
            v12=input('inter v12=   ');
            
            disp('--------------------------------------------');
            v21=E2/E1*v12;
            S=[1/E1 -v12/E1   0;
                -v21/E2 1/E2   0;
                0     0   1/G12];
            Q=inv(S);
        elseif mm==2
            disp('--------------------------------------------');
            Q11=input('inter Q11=  ');
            Q12=input('inter Q12=  ');
            Q13=input('inter Q13=  ');
            Q22=input('inter Q22=  ');
            Q23=input('inter Q23=  ');
            Q33=input('inter Q33=  ');
            disp('--------------------------------------------');
            disp('--------------------------------------------');
            Q=[Q11 Q12 Q13;Q12 Q22 Q23;Q13 Q23 Q33];
        elseif mm==3
            
            mat=input('properties   :    1 ====>  Glass/epoxy (vf=0.45)     2 ====>  boron/epoxy (vf=0.5)     3 ====>  graphite/epoxy (vf=0.7)         : ');
            disp('--------------------------------------------');
            mat_unit=input('  1 ====>SI           2 ====> Psi   :  '   );
            disp('--------------------------------------------');
            % properties >>      Glass/epoxy
            %--------------------SI Glass/epoxy----vf=0.45-------------------
            if mat==1
                if mat_unit==1
                    E1=38.6e9;
                    E2=8.27e9;
                    G12=4.14e9;
                    v12=0.26;
                    Xt=1062e6;
                    Xc=610e6;
                    Yt=31e6;
                    Yc=118e6;
                    SS=72e6;
                elseif mat_unit==2
                    %---------------------Psi Glass/epoxy------vf=0.45----------------
                    E1=5.6e6;
                    E2=1.2e6;
                    G12=0.6e6;
                    v12=0.26;
                    Xt=154.03e3;
                    Xc=88470;
                    Yt=4496;
                    Yc=17120;
                    SS=10440;
                end
            elseif mat==2
                %----------------------SI boron/epoxy-------vf=0.5--------------
                if mat_unit==1
                    E1=204e9;
                    E2=18.5e9;
                    G12=5.59e9;
                    v12=0.23;
                    Xt=1260e6;
                    Xc=2500e6;
                    Yt=61e6;
                    Yc=202e6;
                    SS=67e6;
                    %----------------------Psi boron/epoxy-------vf=0.5--------------
                elseif mat_unit==2
                    E1=29.59e6;
                    E2=2.683e6;
                    G12=0.811e6;
                    v12=0.23;
                    Xt=182.75e3;
                    Xc=362.6e3;
                    Yt=8.847e3;
                    Yc=29.3e3;
                    SS=9.718e3;
                end
            elseif mat==3
                %----------------------SI graphite/epoxy--------vf=0.7-------------
                if mat_unit==1
                    E1=181e9;
                    E2=10.3e9;
                    G12=7.17e9;
                    v12=0.28;
                    Xt=1500e6;
                    Xc=1500e6;
                    Yt=40e6;
                    Yc=246e6;
                    SS=68e6;
                    %----------------------Psi graphite/epoxy----------vf=0.7-----------
                elseif mat_unit==2
                    E1=26.25e6;
                    E2=1.49e6;
                    G12=1.04e6;
                    v12=0.28;
                    Xt=217.56e3;
                    Xc=217.56e3;
                    Yt=5.802e3;
                    Yc=35.68e3;
                    SS=9.863e3;
                end
            end
            
            
            
            v21=E2/E1*v12;
            S=[1/E1 -v12/E1   0;
                -v21/E2 1/E2   0;
                0     0   1/G12];
            Q=inv(S);
            
        end
    else
        material_de=input('material is similar to the previous layer?              1.===>YES     2.===>NO   :  ');
        if material_de==1
            Q=C{7,i-1};
            C{7,i}=Q;
        else
            mm=input('                1===> mechanical properties        2===> local stiffness matrix(Q)         3===>Default propertie  ');
            
            
            
            if mm==1
                disp('--------------------------------------------');
                E1=input('inter E1=   ');
                E2=input('inter E2=   ');
                v12=input('inter v12=   ');
                G12=input('inter G12=   ');
                disp('--------------------------------------------');
                v21=E2/E1*v12;
                S=[1/E1 -v12/E1   0;
                    -v21/E2 1/E2   0;
                    0     0   1/G12];
                Q=inv(S);
            elseif mm==2
                disp('--------------------------------------------');
                Q11=input('inter Q11=  ');
                Q12=input('inter Q12=  ');
                Q13=input('inter Q13=  ');
                Q22=input('inter Q22=  ');
                Q23=input('inter Q23=  ');
                Q33=input('inter Q33=  ');
                disp('--------------------------------------------');
                disp('--------------------------------------------');
                Q=[Q11 Q12 Q13;Q12 Q22 Q23;Q13 Q23 Q33];
            elseif mm==3
                
                mat=input('properties   :    1 ====>  Glass/epoxy (vf=0.45)     2 ====>  boron/epoxy (vf=0.5)     3 ====>  graphite/epoxy (vf=0.7)         : ');
                disp('--------------------------------------------');
                mat_unit=input('  1 ====>SI           2 ====> Psi   :  '   );
                disp('--------------------------------------------');
                % properties >>      Glass/epoxy
                %--------------------SI Glass/epoxy----vf=0.45-------------------
                if mat==1
                    if mat_unit==1
                        E1=38.6e9;
                        E2=8.27e9;
                        G12=4.14e9;
                        v12=0.26;
                        Xt=1062e6;
                        Xc=610e6;
                        Yt=31e6;
                        Yc=118e6;
                        SS=72e6;
                    elseif mat_unit==2
                        %---------------------Psi Glass/epoxy------vf=0.45----------------
                        E1=5.6e6;
                        E2=1.2e6;
                        G12=0.6e6;
                        v12=0.26;
                        Xt=154.03e3;
                        Xc=88470;
                        Yt=4496;
                        Yc=17120;
                        SS=10440;
                    end
                elseif mat==2
                    %----------------------SI boron/epoxy-------vf=0.5--------------
                    if mat_unit==1
                        E1=204e9;
                        E2=18.5e9;
                        G12=5.59e9;
                        v12=0.23;
                        Xt=1260e6;
                        Xc=2500e6;
                        Yt=61e6;
                        Yc=202e6;
                        SS=67e6;
                        %----------------------Psi boron/epoxy-------vf=0.5--------------
                    elseif mat_unit==2
                        E1=29.59e6;
                        E2=2.683e6;
                        G12=0.811e6;
                        v12=0.23;
                        Xt=182.75e3;
                        Xc=362.6e3;
                        Yt=8.847e3;
                        Yc=29.3e3;
                        SS=9.718e3;
                    end
                elseif mat==3
                    %----------------------SI graphite/epoxy--------vf=0.7-------------
                    if mat_unit==1
                        E1=181e9;
                        E2=10.3e9;
                        G12=7.17e9;
                        v12=0.28;
                        Xt=1500e6;
                        Xc=1500e6;
                        Yt=40e6;
                        Yc=246e6;
                        SS=68e6;
                        %----------------------Psi graphite/epoxy----------vf=0.7-----------
                    elseif mat_unit==2
                        E1=26.25e6;
                        E2=1.49e6;
                        G12=1.04e6;
                        v12=0.28;
                        Xt=217.56e3;
                        Xc=217.56e3;
                        Yt=5.802e3;
                        Yc=35.68e3;
                        SS=9.863e3;
                    end
                end
                
                
                
                v21=E2/E1*v12;
                S=[1/E1 -v12/E1   0;
                    -v21/E2 1/E2   0;
                    0     0   1/G12];
                Q=inv(S);
                
            end
            
            
        end
    end
    
    C{7,i}=Q;
    
    
    Tkol=Tkol+t(i,1);
    disp('--------------------------------------------');
end
disp('--------------------------------------------');
disp('--------------------------------------------');
mideplane=Tkol/2;
z=zeros(n+1,1);
z(1,1)=mideplane;
z(n+1,1)=-mideplane;
tt=0;
for i=2:n
    tt=tt+t(i-1,1);
    z(i,1)=mideplane-tt;
end
hold on

for i=1:n+1
    
    x=[-1 +1];
    y=[z(i,1) z(i,1)];
    figure(1)
    line(x,y)
    axis([-2 2  -Tkol +Tkol]);
end
x=[-1 +1];
y=[0 0];
figure(1)
plot(x,y,'--r')


for i=1:n
    Q=C{7,i};
    tetha=t(i,2);
    c=cosd(tetha);
    s=sind(tetha);
    R=[1 0 0;0 1 0;0 0 2];
    T=[c^2 s^2 2*s*c;s^2 c^2 -2*s*c;-s*c s*c c^2-s^2];
    Q_barr=inv(T)*Q*R*T*inv(R);     %{sigma_xy}=[Q_barr]*{epsilon_xy}
    C{1,i}=Q_barr;
    C{2,i}=T;
    z1=-z(i,1);
    z2=-z(i+1,1);
    
    %     disp('----------------------hazf_layeha----------------------');
    if number_delete~=0
        for jjj=1:number_delete
            
            C{1,number_ply_q_zero(1,jjj)}=zeros(3,3);
            
        end
    end
    %     disp('----------------------payane hazf_layeha----------------------');
    Q_barr=C{1,i};
    A=Q_barr*(z2-z1)+A;
    B=Q_barr*(z2^2-z1^2)/2+B;
    D=Q_barr*(z2^3-z1^3)/3+D;
end
ABD=[A B;B D];
NT=zeros(3,1);
MT=zeros(3,1);
NH=zeros(3,1);
MH=zeros(3,1);
alfa=[0;0;0];
beta=[0;0;0];
alfa1=0;
alfa2=0;
beta1=0;
beta2=0;
deltaT=0;
deltaM=0;
for i=1:n
C{19,i}=[0;0;0];
C{18,i}=[0;0;0];
end
disp('----------------------------------------Temprature stress------------------------------------------');
TT=input('tanesh hararti darid 1 ra vared konid= ');
disp('--------------------------------------------');
if TT==1
    TTCH=input('Alfa local dari 1 global dari 2 = ');
    if TTCH==1
        alfa1=input('Alfa1= ');
        alfa2=input('Alfa2= ');
        deltaT=input('DeltaT= ');
        disp('--------------------------------------------');
    else
        alfax=input('AlfaX= ');
        alfay=input('AlfaY= ');
        alfaxy=input('AlfaXY= ');
        deltaT=input('DeltaT= ');
        disp('--------------------------------------------');
    end
    
    for i=1:n
        if TTCH==1
            T=C{2,i};
            alfa11=inv(T)*[alfa1;alfa2;0];
            alfa=[ alfa11(1,1);alfa11(2,1);alfa11(3,1)*2];
            C{18,i}=alfa;
            
        else
            alfa=[alfax;alfay;alfaxy];
            C{18,i}=alfa;
        end
        
        Q=C{1,i};
        z1=-z(i,1);
        z2=-z(i+1,1);
        NT1=Q*alfa*(z2-z1)*deltaT;
        NT=NT+NT1;
        MT1=0.5*Q*alfa*((z2)^2-(z1)^2)*deltaT;
        MT=MT+MT1;
    end
end
disp('---------------------------------------Hygroscopic stress----------------------------------------');
HH=input('tanesh rotoobati darid 1 ra vared konid= ');
disp('--------------------------------------------');
if HH==1
    TTCH=input('Beta local dari 1 global dari 2 = ');
    if TTCH==1
        beta1=input('beta1= ');
        beta2=input('beta2= ');
        deltaM=input('DeltaM= ');
        disp('--------------------------------------------');
    else
        betax=input('betaX= ');
        betay=input('betaY= ');
        betaxy=input('betaXY= ');
        deltaM=input('DeltaM= ');
        disp('--------------------------------------------');
    end
    for i=1:n
        if TTCH==1
            T=C{2,i};
            beta11=inv(T)*[beta1;beta2;0];
            beta=[ beta11(1,1);beta11(2,1);beta11(3,1)*2];
            C{19,i}=beta;
            
        else
            beta=[betax;betay;betaxy];
            C{19,i}=beta;
        end
        Q=C{1,i};
        z1=-z(i,1);
        z2=-z(i+1,1);
        NH1=Q*beta*(z2-z1)*deltaM;
        NH=NH+NH1;
        MH1=0.5*Q*beta*((z2)^2-(z1)^2)*deltaM;
        MH=MH+MH1;
    end
end
disp('----------------------------------------------------------------------------------');
mmm=input('agar Niroo mikhahid adade 1 ya agar kornesh mikhahid addade 2 ra vared konid= ');
if mmm==1
    disp('--------------------------------------------');
    epsilonX=input('inter epsilonX= ');
    epsilonY=input('inter epsilonY= ');
    epsilonXY=input('inter epsilonXY= ');
    KX=input('inter Kx= ');
    KY=input('inter KY= ');
    KXY=input('inter KXY= ');
    disp('--------------------------------------------');
    Epsilon=[epsilonX;epsilonY;epsilonXY;KX;KY;KXY];
    N=ABD*Epsilon;
elseif mmm==2
    disp('--------------------------------------------');
    Nx=input('inter Nx= ');
    Ny=input('inter Ny= ');
    Nxy=input('inter Nxy= ');
    Mx=input('inter Mx= ');
    My=input('inter My= ');
    Mxy=input('inter Mxy= ');
    disp('--------------------------------------------');
    NM=[Nx;Ny;Nxy;Mx;My;Mxy];
    NM=NM+[NT;MT]+[NH;MH];
    Epsilon=inv(ABD)*NM;
end
syms d

for i=1:n
    T=C{2,i};
    alfa=C{18,i};
    beta=C{19,i};
    alfa11=inv(T)*[alfa1;alfa2;0];
    Qbar=C{1,i};
    Stress=Qbar*([Epsilon(1,1);Epsilon(2,1);Epsilon(3,1)]+d*[Epsilon(4,1);Epsilon(5,1);Epsilon(6,1)]-(deltaT*alfa)-(deltaM*beta));
    C{5,i}=vpa(Stress);
    Stress=T*Stress;
    C{6,i}=vpa(Stress);
    strain_xy=([Epsilon(1,1);Epsilon(2,1);Epsilon(3,1)]+d*[Epsilon(4,1);Epsilon(5,1);Epsilon(6,1)]);
    strain_12=T*strain_xy;
    C{12,i}=vpa(strain_xy);
    C{13,i}=vpa(strain_12);
end



% disp('-------------------rasme tanesh ha-------------------------');
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(2);
subplot(2,2,1)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(2);
subplot(2,2,2)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(2);
subplot(2,2,3)
plot(y,x,'--g')
a22=zeros(3,n);
for j=1:n
    z_meti=-z;
    k_a=z_meti(j,1);
    k_e=z_meti(j+1,1);
    increment=(k_e-k_a)/4;
    k_b=(z_meti(j,1)+increment);
    k_c=z_meti(j,1)+2*increment;
    k_d=z_meti(j,1)+3*increment;
    
    Stress_a=subs(C{5,j},d,k_a);
    Stress_b=subs(C{5,j},d,k_b);
    Stress_c=subs(C{5,j},d,k_c);
    Stress_d=subs(C{5,j},d,k_d);
    Stress_e=subs(C{5,j},d,k_e);
    
    C{10,j}=vpa(Stress_a);
    C{11,j}=vpa(Stress_e);
    
    sigma_x_a=simplify(vpa(Stress_a(1,1)));
    sigma_y_a=simplify(Stress_a(2,1));
    sigma_xy_a=simplify(Stress_a(3,1));
    
    sigma_x_b=simplify(vpa(Stress_b(1,1)));
    sigma_y_b=simplify(Stress_b(2,1));
    sigma_xy_b=simplify(Stress_b(3,1));
    
    sigma_x_c=simplify(vpa(Stress_c(1,1)));
    sigma_y_c=simplify(Stress_c(2,1));
    sigma_xy_c=simplify(Stress_c(3,1));
    
    sigma_x_d=simplify(vpa(Stress_d(1,1)));
    sigma_y_d=simplify(Stress_d(2,1));
    sigma_xy_d=simplify(Stress_d(3,1));
    
    sigma_x_e=simplify(vpa(Stress_e(1,1)));
    sigma_y_e=simplify(Stress_e(2,1));
    sigma_xy_e=simplify(Stress_e(3,1));
    a22(1,j)=(max([abs(sigma_x_a),abs(sigma_x_b),abs(sigma_x_c),abs(sigma_x_d),abs(sigma_x_e)])*1.5)+1;
    a22(2,j)=(max([abs(sigma_y_a),abs(sigma_y_b),abs(sigma_y_c),abs(sigma_y_d),abs(sigma_y_e)])*1.5)+1;
    a22(3,j)=(max([abs(sigma_xy_a),abs(sigma_xy_b),abs(sigma_xy_c),abs(sigma_xy_d),abs(sigma_xy_e)])*1.5)+1;
    
    
    x01=[sigma_x_a sigma_x_b];
    y01=[-k_a -k_b];
    x02=[sigma_x_b sigma_x_c];
    y02=[-k_b -k_c];
    x03=[sigma_x_c sigma_x_d];
    y03=[-k_c -k_d];
    x04=[sigma_x_d sigma_x_e];
    y04=[-k_d -k_e];
    figure(2);
    subplot(2,2,1)
    title('\sigma_x')
    hold on
    line(x01,y01);
    line(x02,y02);
    line(x03,y03);
    line(x04,y04);
    
    % disp('----------------------------------------------------------------');
    x012=[sigma_y_a sigma_y_b];
    y012=[-k_a -k_b];
    x022=[sigma_y_b sigma_y_c];
    y022=[-k_b -k_c];
    x032=[sigma_y_c sigma_y_d];
    y032=[-k_c -k_d];
    x042=[sigma_y_d sigma_y_e];
    y042=[-k_d -k_e];
    figure(2);
    subplot(2,2,2)
    title('\sigma_y')
    hold on
    line(x012,y012);
    line(x022,y022);
    line(x032,y032);
    line(x042,y042);
    
    % disp('----------------------------------------------------------------');
    % disp('----------------------------------------------------------------');
    x011=[sigma_xy_a sigma_xy_b];
    y011=[-k_a -k_b];
    x021=[sigma_xy_b sigma_xy_c];
    y021=[-k_b -k_c];
    x031=[sigma_xy_c sigma_xy_d];
    y031=[-k_c -k_d];
    x041=[sigma_xy_d sigma_xy_e];
    y041=[-k_d -k_e];
    figure(2);
    subplot(2,2,3)
    title('\sigma_x_y')
    hold on
    line(x011,y011);
    line(x021,y021);
    line(x031,y031);
    line(x041,y041);
    x=[-1 +1];
    y=[0 0];
    plot(x,y,'--r')
    % disp('----------------------------------------------------------------');
    
    
end
for i=1:n+1
    
    x=[-1 +1];
    y=[z(i,1) z(i,1)];
    subplot(2,2,4)
    title('composite layer')
    hold on
    line(x,y)
    axis([-2 2  -Tkol +Tkol]);
end
x=[-1 +1];
y=[0 0];
subplot(2,2,4)
title('composite layer')
hold on
plot(x,y,'--r')
for i=1:n
    a1(i,1)=a22(1,i);
    b1(i,1)=a22(2,i);
    c1(i,1)=a22(3,i);
end
max_stress_xy_x=max(a1);
max_stress_xy_y=max(b1);
max_stress_xy_xy=max(c1);
for i=1:n+1
    x=[-max_stress_xy_x +max_stress_xy_x];
    y=[z(i,1) z(i,1)];
    figure(2);
    subplot(2,2,1)
    line(x,y)
    
end
x=[-max_stress_xy_x +max_stress_xy_x];
y=[0 0];
figure(2);
hold on
subplot(2,2,1)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_xy_y +max_stress_xy_y];
    y=[z(i,1) z(i,1)];
    figure(2);
    subplot(2,2,2)
    hold on
    line(x,y)
    
end
x=[-max_stress_xy_y +max_stress_xy_y];
y=[0 0];
figure(2);
subplot(2,2,2)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_xy_xy +max_stress_xy_xy];
    y=[z(i,1) z(i,1)];
    figure(2);
    subplot(2,2,3)
    hold on
    line(x,y)
    
end
x=[-max_stress_xy_xy +max_stress_xy_xy];
y=[0 0];
figure(2);
subplot(2,2,3)
plot(x,y,'--r')
% [sigma_taghribi,ertefa_taghribi] = ginput(5);

% disp([vpa(sigma_taghribi) vpa(ertefa_taghribi)])


% disp('----------------------------------------------------------------');
% disp('----------------------------------------------------------------');

x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(3);
subplot(2,2,1)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(3);
subplot(2,2,2)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(3);
subplot(2,2,3)
plot(y,x,'--g')
a22=zeros(3,n);
for j=1:n
    z_meti=-z;
    k_a=z_meti(j,1);
    k_e=z_meti(j+1,1);
    increment=(k_e-k_a)/4;
    k_b=(z_meti(j,1)+increment);
    k_c=z_meti(j,1)+2*increment;
    k_d=z_meti(j,1)+3*increment;
    T=C{2,j};
    stress_global=C{5,j};
    stress_local=T*stress_global;
    
    Stress_a=subs(stress_local,d,k_a);
    Stress_b=subs(stress_local,d,k_b);
    Stress_c=subs(stress_local,d,k_c);
    Stress_d=subs(stress_local,d,k_d);
    Stress_e=subs(stress_local,d,k_e);
    
    C{3,j}=vpa(Stress_a);
    C{4,j}=vpa(Stress_e);
    
    sigma_x_a=simplify(vpa(Stress_a(1,1)));
    sigma_y_a=simplify(Stress_a(2,1));
    sigma_xy_a=simplify(Stress_a(3,1));
    
    sigma_x_b=simplify(vpa(Stress_b(1,1)));
    sigma_y_b=simplify(Stress_b(2,1));
    sigma_xy_b=simplify(Stress_b(3,1));
    
    sigma_x_c=simplify(vpa(Stress_c(1,1)));
    sigma_y_c=simplify(Stress_c(2,1));
    sigma_xy_c=simplify(Stress_c(3,1));
    
    sigma_x_d=simplify(vpa(Stress_d(1,1)));
    sigma_y_d=simplify(Stress_d(2,1));
    sigma_xy_d=simplify(Stress_d(3,1));
    
    sigma_x_e=simplify(vpa(Stress_e(1,1)));
    sigma_y_e=simplify(Stress_e(2,1));
    sigma_xy_e=simplify(Stress_e(3,1));
    a22(1,j)=(max([abs(sigma_x_a),abs(sigma_x_b),abs(sigma_x_c),abs(sigma_x_d),abs(sigma_x_e)])*1.5)+1;
    a22(2,j)=(max([abs(sigma_y_a),abs(sigma_y_b),abs(sigma_y_c),abs(sigma_y_d),abs(sigma_y_e)])*1.5)+1;
    a22(3,j)=(max([abs(sigma_xy_a),abs(sigma_xy_b),abs(sigma_xy_c),abs(sigma_xy_d),abs(sigma_xy_e)])*1.5)+1;
    
    
    x01=[sigma_x_a sigma_x_b];
    y01=[-k_a -k_b];
    x02=[sigma_x_b sigma_x_c];
    y02=[-k_b -k_c];
    x03=[sigma_x_c sigma_x_d];
    y03=[-k_c -k_d];
    x04=[sigma_x_d sigma_x_e];
    y04=[-k_d -k_e];
    figure(3);
    subplot(2,2,1)
    title('\sigma_1')
    hold on
    line(x01,y01);
    line(x02,y02);
    line(x03,y03);
    line(x04,y04);
    
    % disp('----------------------------------------------------------------');
    x012=[sigma_y_a sigma_y_b];
    y012=[-k_a -k_b];
    x022=[sigma_y_b sigma_y_c];
    y022=[-k_b -k_c];
    x032=[sigma_y_c sigma_y_d];
    y032=[-k_c -k_d];
    x042=[sigma_y_d sigma_y_e];
    y042=[-k_d -k_e];
    figure(3);
    subplot(2,2,2)
    title('\sigma_2')
    hold on
    line(x012,y012);
    line(x022,y022);
    line(x032,y032);
    line(x042,y042);
    
    % disp('----------------------------------------------------------------');
    
    % disp('----------------------------------------------------------------');
    x011=[sigma_xy_a sigma_xy_b];
    y011=[-k_a -k_b];
    x021=[sigma_xy_b sigma_xy_c];
    y021=[-k_b -k_c];
    x031=[sigma_xy_c sigma_xy_d];
    y031=[-k_c -k_d];
    x041=[sigma_xy_d sigma_xy_e];
    y041=[-k_d -k_e];
    figure(3);
    subplot(2,2,3)
    title('\sigma_1_2')
    hold on
    line(x011,y011);
    line(x021,y021);
    line(x031,y031);
    line(x041,y041);
    x=[-1 +1];
    y=[0 0];
    plot(x,y,'--r')
    % disp('----------------------------------------------------------------');
    
    
end
for i=1:n+1
    
    x=[-1 +1];
    y=[z(i,1) z(i,1)];
    subplot(2,2,4)
    title('composite layer')
    hold on
    line(x,y)
    axis([-2 2  -Tkol +Tkol]);
end
x=[-1 +1];
y=[0 0];
subplot(2,2,4)
title('composite layer')
hold on
plot(x,y,'--r')
for i=1:n
    a1(i,1)=a22(1,i);
    b1(i,1)=a22(2,i);
    c1(i,1)=a22(3,i);
end

max_stress_12_1=max(a1);
max_stress_12_2=max(b1);
max_stress_12_12=max(c1);
for i=1:n+1
    x=[-max_stress_12_1 +max_stress_12_1];
    y=[z(i,1) z(i,1)];
    figure(3);
    subplot(2,2,1)
    line(x,y)
    
end
x=[-max_stress_12_1 +max_stress_12_1];
y=[0 0];
figure(3);
hold on
subplot(2,2,1)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_12_2 +max_stress_12_2];
    y=[z(i,1) z(i,1)];
    figure(3);
    subplot(2,2,2)
    hold on
    line(x,y)
    
end
x=[-max_stress_12_2 +max_stress_12_2];
y=[0 0];
figure(3);
subplot(2,2,2)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_12_12 +max_stress_12_12];
    y=[z(i,1) z(i,1)];
    figure(3);
    subplot(2,2,3)
    hold on
    line(x,y)
    
end
x=[-max_stress_12_12 +max_stress_12_12];
y=[0 0];
figure(3);
subplot(2,2,3)
plot(x,y,'--r')
% [sigma_taghribi,ertefa_taghribi] = ginput(5);

% disp([vpa(sigma_taghribi) vpa(ertefa_taghribi)])




% disp('-------------------rasme strain ha-------------------------');
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(4);
subplot(2,2,1)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(4);
subplot(2,2,2)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(4);
subplot(2,2,3)
plot(y,x,'--g')
a22=zeros(3,n);
for j=1:n
    z_meti=-z;
    k_a=z_meti(j,1);
    k_e=z_meti(j+1,1);
    increment=(k_e-k_a)/4;
    k_b=(z_meti(j,1)+increment);
    k_c=z_meti(j,1)+2*increment;
    k_d=z_meti(j,1)+3*increment;
    
    Stress_a=subs(C{12,j},d,k_a);
    Stress_b=subs(C{12,j},d,k_b);
    Stress_c=subs(C{12,j},d,k_c);
    Stress_d=subs(C{12,j},d,k_d);
    Stress_e=subs(C{12,j},d,k_e);
    
    C{14,j}=vpa(Stress_a);
    C{15,j}=vpa(Stress_e);
    
    sigma_x_a=simplify(vpa(Stress_a(1,1)));
    sigma_y_a=simplify(Stress_a(2,1));
    sigma_xy_a=simplify(Stress_a(3,1));
    
    sigma_x_b=simplify(vpa(Stress_b(1,1)));
    sigma_y_b=simplify(Stress_b(2,1));
    sigma_xy_b=simplify(Stress_b(3,1));
    
    sigma_x_c=simplify(vpa(Stress_c(1,1)));
    sigma_y_c=simplify(Stress_c(2,1));
    sigma_xy_c=simplify(Stress_c(3,1));
    
    sigma_x_d=simplify(vpa(Stress_d(1,1)));
    sigma_y_d=simplify(Stress_d(2,1));
    sigma_xy_d=simplify(Stress_d(3,1));
    
    sigma_x_e=simplify(vpa(Stress_e(1,1)));
    sigma_y_e=simplify(Stress_e(2,1));
    sigma_xy_e=simplify(Stress_e(3,1));
    a22(1,j)=(max([abs(sigma_x_a),abs(sigma_x_b),abs(sigma_x_c),abs(sigma_x_d),abs(sigma_x_e)])*1.5);
    a22(2,j)=(max([abs(sigma_y_a),abs(sigma_y_b),abs(sigma_y_c),abs(sigma_y_d),abs(sigma_y_e)])*1.5);
    a22(3,j)=(max([abs(sigma_xy_a),abs(sigma_xy_b),abs(sigma_xy_c),abs(sigma_xy_d),abs(sigma_xy_e)])*1.5);
    
    
    x01=[sigma_x_a sigma_x_b];
    y01=[-k_a -k_b];
    x02=[sigma_x_b sigma_x_c];
    y02=[-k_b -k_c];
    x03=[sigma_x_c sigma_x_d];
    y03=[-k_c -k_d];
    x04=[sigma_x_d sigma_x_e];
    y04=[-k_d -k_e];
    figure(4);
    subplot(2,2,1)
    title('\epsilon_x')
    hold on
    line(x01,y01);
    line(x02,y02);
    line(x03,y03);
    line(x04,y04);
    
    % disp('----------------------------------------------------------------');
    x012=[sigma_y_a sigma_y_b];
    y012=[-k_a -k_b];
    x022=[sigma_y_b sigma_y_c];
    y022=[-k_b -k_c];
    x032=[sigma_y_c sigma_y_d];
    y032=[-k_c -k_d];
    x042=[sigma_y_d sigma_y_e];
    y042=[-k_d -k_e];
    figure(4);
    subplot(2,2,2)
    title('\epsilon_y')
    hold on
    line(x012,y012);
    line(x022,y022);
    line(x032,y032);
    line(x042,y042);
    
    % disp('----------------------------------------------------------------');
    
    % disp('----------------------------------------------------------------');
    x011=[sigma_xy_a sigma_xy_b];
    y011=[-k_a -k_b];
    x021=[sigma_xy_b sigma_xy_c];
    y021=[-k_b -k_c];
    x031=[sigma_xy_c sigma_xy_d];
    y031=[-k_c -k_d];
    x041=[sigma_xy_d sigma_xy_e];
    y041=[-k_d -k_e];
    figure(4);
    subplot(2,2,3)
    title('\epsilon_x_y')
    hold on
    line(x011,y011);
    line(x021,y021);
    line(x031,y031);
    line(x041,y041);
    %     x=[-1 +1];
    %     y=[0 0];
    %     plot(x,y,'--r')
    % disp('----------------------------------------------------------------');
    
    
end
for i=1:n+1
    
    x=[-1 +1];
    y=[z(i,1) z(i,1)];
    subplot(2,2,4)
    title('composite layer')
    hold on
    line(x,y)
    axis([-2 2  -Tkol +Tkol]);
end
x=[-1 +1];
y=[0 0];
subplot(2,2,4)
title('composite layer')
hold on
plot(x,y,'--r')
for i=1:n
    a1(i,1)=a22(1,i);
    b1(i,1)=a22(2,i);
    c1(i,1)=a22(3,i);
end
max_stress_xy_x=max(a1)*1.5;
max_stress_xy_y=max(b1)*1.5;
max_stress_xy_xy=max(c1)*1.5;

if max_stress_12_12==0
    max_stress_xy_xy=max([max_stress_12_2 max_stress_12_1]);
end

for i=1:n+1
    x=[-max_stress_xy_x +max_stress_xy_x];
    y=[z(i,1) z(i,1)];
    figure(4);
    subplot(2,2,1)
    line(x,y)
    
end
x=[-max_stress_xy_x +max_stress_xy_x];
y=[0 0];
figure(4);
hold on
subplot(2,2,1)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_xy_y +max_stress_xy_y];
    y=[z(i,1) z(i,1)];
    figure(4);
    subplot(2,2,2)
    hold on
    line(x,y)
    
end
x=[-max_stress_xy_y +max_stress_xy_y];
y=[0 0];
figure(4);
subplot(2,2,2)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_xy_xy +max_stress_xy_xy];
    y=[z(i,1) z(i,1)];
    figure(4);
    subplot(2,2,3)
    hold on
    line(x,y)
    
end
x=[-max_stress_xy_xy +max_stress_xy_xy];
y=[0 0];
figure(4);
subplot(2,2,3)
plot(x,y,'--r')
% [sigma_taghribi,ertefa_taghribi] = ginput(5);

% disp([vpa(sigma_taghribi) vpa(ertefa_taghribi)])


% disp('----------------------------------------------------------------');
% disp('----------------------------------------------------------------');

x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(5);
subplot(2,2,1)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(5);
subplot(2,2,2)
plot(y,x,'--g')
x=[-Tkol/2 +Tkol/2];
y=[0 0];
figure(5);
subplot(2,2,3)
plot(y,x,'--g')
a22=zeros(3,n);
for j=1:n
    z_meti=-z;
    k_a=z_meti(j,1);
    k_e=z_meti(j+1,1);
    increment=(k_e-k_a)/4;
    k_b=(z_meti(j,1)+increment);
    k_c=z_meti(j,1)+2*increment;
    k_d=z_meti(j,1)+3*increment;
    T=C{2,j};
    stress_local=C{13,j};
    
    Stress_a=subs(stress_local,d,k_a);
    Stress_b=subs(stress_local,d,k_b);
    Stress_c=subs(stress_local,d,k_c);
    Stress_d=subs(stress_local,d,k_d);
    Stress_e=subs(stress_local,d,k_e);
    
    C{16,j}=vpa(Stress_a);
    C{17,j}=vpa(Stress_e);
    
    sigma_x_a=simplify(vpa(Stress_a(1,1)));
    sigma_y_a=simplify(Stress_a(2,1));
    sigma_xy_a=simplify(Stress_a(3,1));
    
    sigma_x_b=simplify(vpa(Stress_b(1,1)));
    sigma_y_b=simplify(Stress_b(2,1));
    sigma_xy_b=simplify(Stress_b(3,1));
    
    sigma_x_c=simplify(vpa(Stress_c(1,1)));
    sigma_y_c=simplify(Stress_c(2,1));
    sigma_xy_c=simplify(Stress_c(3,1));
    
    sigma_x_d=simplify(vpa(Stress_d(1,1)));
    sigma_y_d=simplify(Stress_d(2,1));
    sigma_xy_d=simplify(Stress_d(3,1));
    
    sigma_x_e=simplify(vpa(Stress_e(1,1)));
    sigma_y_e=simplify(Stress_e(2,1));
    sigma_xy_e=simplify(Stress_e(3,1));
    a22(1,j)=(max([abs(sigma_x_a),abs(sigma_x_b),abs(sigma_x_c),abs(sigma_x_d),abs(sigma_x_e)])*1.5);
    a22(2,j)=(max([abs(sigma_y_a),abs(sigma_y_b),abs(sigma_y_c),abs(sigma_y_d),abs(sigma_y_e)])*1.5);
    a22(3,j)=(max([abs(sigma_xy_a),abs(sigma_xy_b),abs(sigma_xy_c),abs(sigma_xy_d),abs(sigma_xy_e)])*1.5);
    
    
    x01=[sigma_x_a sigma_x_b];
    y01=[-k_a -k_b];
    x02=[sigma_x_b sigma_x_c];
    y02=[-k_b -k_c];
    x03=[sigma_x_c sigma_x_d];
    y03=[-k_c -k_d];
    x04=[sigma_x_d sigma_x_e];
    y04=[-k_d -k_e];
    figure(5);
    subplot(2,2,1)
    title('\epsilon_1')
    hold on
    line(x01,y01);
    line(x02,y02);
    line(x03,y03);
    line(x04,y04);
    
    % disp('----------------------------------------------------------------');
    x012=[sigma_y_a sigma_y_b];
    y012=[-k_a -k_b];
    x022=[sigma_y_b sigma_y_c];
    y022=[-k_b -k_c];
    x032=[sigma_y_c sigma_y_d];
    y032=[-k_c -k_d];
    x042=[sigma_y_d sigma_y_e];
    y042=[-k_d -k_e];
    figure(5);
    subplot(2,2,2)
    title('\epsilon_2')
    hold on
    line(x012,y012);
    line(x022,y022);
    line(x032,y032);
    line(x042,y042);
    
    % disp('----------------------------------------------------------------');
    
    % disp('----------------------------------------------------------------');
    x011=[sigma_xy_a sigma_xy_b];
    y011=[-k_a -k_b];
    x021=[sigma_xy_b sigma_xy_c];
    y021=[-k_b -k_c];
    x031=[sigma_xy_c sigma_xy_d];
    y031=[-k_c -k_d];
    x041=[sigma_xy_d sigma_xy_e];
    y041=[-k_d -k_e];
    figure(5);
    subplot(2,2,3)
    title('\epsilon_1_2')
    hold on
    line(x011,y011);
    line(x021,y021);
    line(x031,y031);
    line(x041,y041);
    %     x=[-1 +1];
    %     y=[0 0];
    %     plot(x,y,'--r')
    % disp('----------------------------------------------------------------');
    
    
end
for i=1:n+1
    
    x=[-1 +1];
    y=[z(i,1) z(i,1)];
    subplot(2,2,4)
    title('composite layer')
    hold on
    line(x,y)
    axis([-2 2  -Tkol +Tkol]);
end
x=[-1 +1];
y=[0 0];
subplot(2,2,4)
title('composite layer')
hold on
plot(x,y,'--r')
for i=1:n
    a1(i,1)=a22(1,i);
    b1(i,1)=a22(2,i);
    c1(i,1)=a22(3,i);
end

max_stress_12_1=max(a1)*1.5;
max_stress_12_2=max(b1)*1.5;
max_stress_12_12=max(c1)*1.5;
if max_stress_12_12==0
    max_stress_12_12=max([max_stress_12_2 max_stress_12_1]);
end
for i=1:n+1
    x=[-max_stress_12_1 +max_stress_12_1];
    y=[z(i,1) z(i,1)];
    figure(5);
    subplot(2,2,1)
    line(x,y)
    
end
x=[-max_stress_12_1 +max_stress_12_1];
y=[0 0];
figure(5);
hold on
subplot(2,2,1)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_12_2 +max_stress_12_2];
    y=[z(i,1) z(i,1)];
    figure(5);
    subplot(2,2,2)
    hold on
    line(x,y)
    
end
x=[-max_stress_12_2 +max_stress_12_2];
y=[0 0];
figure(5);
subplot(2,2,2)
plot(x,y,'--r')
for i=1:n+1
    
    x=[-max_stress_12_12 +max_stress_12_12];
    y=[z(i,1) z(i,1)];
    figure(5);
    subplot(2,2,3)
    hold on
    line(x,y)
    
end
x=[-max_stress_12_12 +max_stress_12_12];
y=[0 0];
figure(5);
subplot(2,2,3)
plot(x,y,'--r')
figure(2);
title('\sigma_x_y_ _m_e_c_h_a_n_i_c_a_l')
figure(3);
title('\sigma_1_2_ _m_e_c_h_a_n_i_c_a_l')
figure(4);
title('\epsilon_x_y_ _t_o_t_a_l')
figure(5);
title('\epsilon_1_2_ _t_o_t_a_l')

disp('-----------------------------------------criteria-----------------------------------------');

% disp('------------------------------------criteria-------------------------------');
% disp('----------------------------------modifiled hill------------------------------------');
disp('--------------------------------------------');
criteria_method=input('select criteria_method :         modify hill  =====> 1     Tsai_Wu =====> 2   hashin =====> 3    :  ');
if criteria_method==1
    if mm==3
        xc=input('1====>  new_mat       2====> default mat');
        if  xc==1
            Xt=input('please inter Xt        : ');
            
            Xc=input('please inter Xc   (mosbatesho )     : ');
            
            Yt=input('please inter Yt        : ');
            
            Yc=input('please inter Yc    (mosbatesho )    : ');
            
            SS=input('please inter S        : ');
        end
    else
        Xt=input('please inter Xt        : ');
        
        Xc=input('please inter Xc   (mosbatesho )     : ');
        
        Yt=input('please inter Yt        : ');
        
        Yc=input('please inter Yc    (mosbatesho )    : ');
        
        SS=input('please inter S        : ');
    end
    disp('--------------------------------------------');
    for i=1:n
        Stress_up=C{3,i};
        Stress_down= C{4,i};
        Stress=(Stress_down+Stress_up)/2;  %miyangin tanesh bala va paein laye
        
        if Stress(1,1)>0
            X1=Xt;
        else
            X1=Xc;
        end
        if Stress(2,1)>0
            X2=Xt;
            Y=Yt;
        else
            X2=Xc;
            Y=Yc;
        end
        Hill(i,1)=(Stress(1,1)/X1)^2 - (Stress(1,1)*Stress(2,1))/X2^2 + (Stress(2,1)/Y)^2 + (Stress(3,1)/SS)^2 ;
    end
    
    [max_hill,layer_num]=max(Hill);
    layers_num=find(Hill==max_hill);
    SR=1/sqrt(max_hill);
    disp('--------------------------------------------');
    disp(['          max_modify_hill','          SR_modify_Hill     ','           Layer_num']);
    disp([vpa(max_hill),vpa(SR),[layers_num]]);
    disp('--------------------------------------------');
    
    % disp('----------------------------------------------------------------------');
    % disp('----------------------------------TSAI_WU------------------------------------');
elseif criteria_method==2
    
    if mm==3
        xc=input('1====>  new_mat       2====> default mat');
        if  xc==1
            Xt=input('please inter Xt        : ');
            
            Xc=input('please inter Xc   (mosbatesho )     : ');
            
            Yt=input('please inter Yt        : ');
            
            Yc=input('please inter Yc    (mosbatesho )    : ');
            
            SS=input('please inter S        : ');
        end
    else
        Xt=input('please inter Xt        : ');
        
        Xc=input('please inter Xc   (mosbatesho )     : ');
        
        Yt=input('please inter Yt        : ');
        
        Yc=input('please inter Yc    (mosbatesho )    : ');
        
        SS=input('please inter S        : ');
    end
    disp('--------------------------------------------------------------------------------------');
    F12=input('F12    =   ');
    disp('--------------------------------------------------------------------------------------');
    %     zzz=input('agar zarib bayad hesab koni 1 dar gheir in soorat 2   : ');
    disp('--------------------------------------------------------------------------------------');
    %     if zzz==1
    for i=1:n
        Stress_up=C{3,i};
        syms s_mojaz
        Stress_down= C{4,i};
        Stress=s_mojaz*(Stress_down+Stress_up)/2;  %miyangin tanesh bala va paein laye
        F1=(1/Xt)-(1/Xc);
        F2=(1/Yt)-(1/Yc);
        F11=1/(Xt*Xc);
        F22=1/(Yt*Yc);
        F66=1/(SS)^2;
        F6=0;
        %    TSAI_WU(i,1)=([F1 F2 F6]*[Stress(1,1);Stress(2,1);Stress(3,1)])+([Stress(1,1) Stress(2,1) Stress(3,1)]*[F11 F12 0;F12 F22 0;0 0 F66]*[Stress(1,1); Stress(2,1);Stress(3,1)]);
        TSAI_WU(i,1)=(F1*Stress(1,1))+(F2*Stress(2,1))+(F6*Stress(3,1))+(F11*(Stress(1,1))^2)+(F22*(Stress(2,1))^2)+(F66*(Stress(3,1))^2)+(2*F12*(Stress(1,1)*(Stress(2,1))));
        if        TSAI_WU(i,1)==0
            baze_mojaz(i,:)=0;
        else
            baze_mojaz(i,:)=solve(TSAI_WU(i,1)==1,s_mojaz);
        end
        
    end
  disp('SR) = ');
    disp( baze_mojaz(:,2));
    %     elseif zzz==2
    %         for i=1:n
    %             Stress_up=C{3,i};
    %             Stress_down= C{4,i};
    %             Stress=(Stress_down+Stress_up)/2;  %miyangin tanesh bala va paein laye
    %             F1=(1/Xt)-(1/Xc);
    %             F2=(1/Yt)-(1/Yc);
    %             F11=1/(Xt*Xc);
    %             F22=1/(Yt*Yc);
    %             F66=1/(SS)^2;
    %             F6=0;
    %             %TSAI_WU(i,1)=([F1 F2 F6]*[Stress(1,1);Stress(2,1);Stress(3,1)])+([Stress(1,1) Stress(2,1) Stress(3,1)]*[F11 F12 0;F12 F22 0;0 0 F66]*[Stress(1,1); Stress(2,1);Stress(3,1)]);
    %             TSAI_WU(i,1)=(F1*Stress(1,1))+(F2*Stress(2,1))+(F6*Stress(3,1))+(F11*(Stress(1,1))^2)+(F22*(Stress(2,1))^2)+(F66*(Stress(3,1))^2)+(2*F12*(Stress(1,1)*(Stress(2,1))));
    %         end
    TSAI_WU=subs(TSAI_WU,s_mojaz,1);
    [max_TSAI_WU,layer_num]=max(TSAI_WU);
    layers_num=find(TSAI_WU==max_TSAI_WU);
    SR=1/(max_TSAI_WU);
    disp('--------------------------------------------');
    disp(['   MAX_value_TSAI_WU     ','     Layer_num']);
    disp([vpa(max_TSAI_WU),[layers_num]]);
    disp('--------------------agar chand satr dashtim Avalin Adad neshan dade shode max tanesh baraye tamam layehahaye sotun dovom ast. ------------------------');
    %     end
elseif criteria_method==3
    if mm==3
        xc=input('1====>  new_mat       2====> default mat');
        if  xc==1
            Xt=input('please inter Xt        : ');
            
            Xc=input('please inter Xc   (mosbatesho )     : ');
            
            Yt=input('please inter Yt        : ');
            
            Yc=input('please inter Yc    (mosbatesho )    : ');
            
            SS=input('please inter S        : ');
        end
    else
        Xt=input('please inter Xt        : ');
        
        Xc=input('please inter Xc   (mosbatesho )     : ');
        
        Yt=input('please inter Yt        : ');
        
        Yc=input('please inter Yc    (mosbatesho )    : ');
        
        SS=input('please inter S        : ');
    end
    for i=1:n
        Stress_up=C{3,i};
        Stress_down= C{4,i};
        Stress=(Stress_down+Stress_up)/2;  %miyangin tanesh bala va paein laye
        if  Stress(1,1)>=0
            hashin(1,i)=(Stress(1,1)/Xt)^2+(Stress(3,1)/SS)^2;
        else
            hashin(1,i)=(Stress(1,1)/Xc)^2;
        end
        if Stress(2,1)>=0
            hashin(2,i)=(Stress(2,1)/Yt)^2+(Stress(3,1)/SS)^2;
        else
            hashin(2,i)=((Yc/(2*SS))^2-1)*(Stress(2,1)/Yc)+(Stress(2,1)/(2*SS))^2+(Stress(3,1)/SS)^2;
        end
    end
    ma=max(hashin);
    maa=max(ma);
    nn=find(hashin==maa);
    nnn=round(nn/2);
    disp('--------------------------------------------');
    disp('laye pokide shode   :');
    disp(nnn);
    disp('--------------------------------------------');
    Stress_up=C{3,nnn};
    disp('--------------------------------------------');
    Stress_down= C{4,nnn};
    Stress=(Stress_down+Stress_up)/2;  %miyangin tanesh bala va paein laye
    if rem(nn,2)==1
        
        if Stress(1,1)>=0
            disp('--------------------------------------------');
            disp('mode takhrib  :  tensile fiber')
            mode_takhrib=1;
            disp('--------------------------------------------');
            disp('--------------------------------------------');
        else
            disp('--------------------------------------------');
            disp('mode takhrib  :  compressive fiber')
            mode_takhrib=2;
            disp('--------------------------------------------');
            disp('--------------------------------------------');
        end
    else
        if Stress(2,1)>=0
            disp('--------------------------------------------');
            disp('mode takhrib  :  tensile matrix')
            mode_takhrib=3;
            disp('--------------------------------------------');
            disp('--------------------------------------------');
        else
            disp('--------------------------------------------');
            disp('mode takhrib  :  compressive matrix')
            mode_takhrib=4;
            disp('--------------------------------------------');
            disp('--------------------------------------------');
        end
    end
    %     zzzz=input('agar zarib bayad hesab koni 1 dar gheir in soorat 2   : ');
    disp('--------------------------------------------------------------------------------------');
    %     if zzzz==1
    if mode_takhrib==4
        syms SQ
        SR_Hashin=solve(((Yc/(2*SS))^2-1)*(Stress(2,1)/Yc)*SQ+(Stress(2,1)/(2*SS))^2*SQ^2+(Stress(3,1)/SS)^2*SQ^2==1,SQ);
    else
        SR_Hashin=1/sqrt(maa);
    end
    %         disp('--------------------------------------------');
    %         disp(['   SR      ','     Layer_num']);
    %         disp([vpa(SR_Hashin),nnn]);
    disp('--------------------------------------------');
    %     else
    %         SR_Hashin=1/maa;
    disp('--------------------------------------------');
    disp(['             SR               ','             max_value_hashin    ','      fail_Layer_num']);
    disp([vpa(SR_Hashin),vpa(maa),nnn]);
    disp('--------------------agar chand satr dashtim Avalin Adad neshan dade shode max tanesh baraye tamam layehahaye sotun dovom ast. ------------------------');
    disp('-----------------------End---------------------');
    
    %     end
end



max_stress_xy_x=max(a1)*1.5;
max_stress_xy_y=max(b1)*1.5;
max_stress_xy_xy=max(c1)*1.5;



if max_stress_xy_x==0
    max_stress_xy_x=1;
end

if max_stress_xy_y==0
    max_stress_xy_y=1;
end

if max_stress_xy_xy==0
    max_stress_xy_xy=1;
end

