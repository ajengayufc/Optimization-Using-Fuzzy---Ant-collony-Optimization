
% selanjutnya dilakukan optimasi dengan ACO
MaxIt=100;                     % max numberof iteration
m1=1;
n1=1;
mse=3;                          % radius
nAnt=40;                       % Number of Ant
Q=3;                           % Konstanta Pheromone Update/tetapan siklus semut
alpha=0.72;                       % tetapan pengendali intensitas jejak semut
beta=0.66;                        % tetapan pengendali visibilitas
rho=0.075;                      % koefisien penguapan / evaporasi

n=1;
% ====================
% Batas-batas untuk 1
i=1;
a(:,i)=-4.617; 
b(:,i)=245.3; 
c(:,i)=559.6; 
% Batas-batas untuk 2
d(:,i)=370.3;
e(:,i)=694.7;
f(:,i)=995.5;    
% Batas-batas untuk 3
g(:,i)=829.8;
h(:,i)=1145;
j(:,i)=1495;

% ====================

% Inisialisasi posisi awal
x1(:,i)=(c(:,i)-(rand(n,1)*(c(:,i)-b(:,i)))); %constraint: b<=x1<=c
x2(:,i)=(d(:,i)-(rand(n,1)*(d(:,i)-a(:,i)))); %constraint: a<=x2<=d
x3(:,i)=(f(:,i)-(rand(n,1)*(f(:,i)-e(:,i)))); %constraint: e<=x3<=f
x4(:,i)=(h(:,i)-(rand(n,1)*(h(:,i)-g(:,i)))); %constraint: g<=x4<=h

% ===================
% Selanjutnya menghitung MSE [mean square error] yang dihasilkan oleh setiap partikel
% "myu" [derajat keanggotaan] tiap dimensi
myu1(:,i)=1;                                % myu1=derajat keanggotaan 1 batas kiri (a) (tetap)
myu2(:,i)=1;                                % myu2=derajat keanggotaan 1 batas tengah (b) (tetap)
myu3(:,i)=(c(:,i)-x1(:,i))/(c(:,i)-b(:,i)); % myu3=derajat keanggotaan 1 batas kanan (c)
y1=1;                                       % y1=derajat keanggotaan ideal 1 batas kiri (a) adalah 1 (tetap)
y2=1;                                       % y2=derajat keanggotaan ideal 1 batas tengah (b) adalah 1 (tetap)
y3=0;                                       % y3=derajat keanggotaan ideal 1 batas kanan (c) adalah 0
 
myu4(:,i)=(x2(:,i)-d(:,i))/(e(:,i)-d(:,i)); % myu4=derajat keanggotaan 2 batas kiri (d)
myu5(:,i)=1;                                % myu5=derajat keanggotaan 2 batas tengah (e) (tetap)
myu6(:,i)=(f(:,i)-x3(:,i))/(f(:,i)-e(:,i)); % myu6=derajat keanggotaan 2 batas kanan (f) 
y4=0;                                       % y4=derajat keanggotaan ideal 2 batas kiri (d) adalah 0
y5=1;                                       % y5=derajat keanggotaan ideal 2 batas tengah (e) adalah 1 (tetap)
y6=0;                                       % y6=derajat keanggotaan ideal 2 batas kanan (f) adalah 0
 
myu7(:,i)=(x4(:,i)-g(:,i))/(h(:,i)-g(:,i)); % myu7=derajat keanggotaan 3 batas kiri (g)
myu8(:,i)=1;                                % myu8=derajat keanggotaan 3 batas tengah (h) (tetap)
myu9(:,i)=1;                                % myu9=derajat keanggotaan 3 batas kanan (j) 
y7=0;                                       % y7=derajat keanggotaan ideal 3 batas kiri (g) adalah 0
y8=1;                                       % y8=derajat keanggotaan ideal 3 batas tengah (h) adalah 1 (tetap)
y9=0;                                       % y9=derajat keanggotaan ideal 3 batas kanan (j) adalah 0
 

mse(:,i)=(((y1-myu1(:,i)).^2)+((y2-myu2(:,i)).^2)+((y3-myu3(:,i)).^2)+((y4-myu4(:,i)).^2)+((y5-myu5(:,i)).^2)+...
    ((y6-myu6(:,i)).^2)+((y7-myu7(:,i)).^2)+((y8-myu8(:,i)).^2)+((y9-myu9(:,i)).^2))/9;

%====================================
%antPosition(AP)
AP1=x1;
AP2=x2;
AP3=x3;
AP4=x4;
%====================================
%jarak tersekat xi, jarak terjauh xj
% x1
xi1(:,i)=AP1(:,i)+(mse(:,i)*rho(:,i));
xj1(:,i)=AP1(:,i)-(mse(:,i)*rho(:,i));
%x2
xi2(:,i)=AP2(:,i)+(mse(:,i)*rho(:,i));
xj2(:,i)=AP2(:,i)-(mse(:,i)*rho(:,i));
%x3
xi3(:,i)=AP3(:,i)+(mse(:,i)*rho(:,i));
xj3(:,i)=AP3(:,i)-(mse(:,i)*rho(:,i));
%x4
xi4(:,i)=AP4(:,i)+(mse(:,i)*rho(:,i));
xj4(:,i)=AP4(:,i)-(mse(:,i)*rho(:,i));

%% jarak antar kota

% city position(CP)
CP1(:,i)=[xi1(:,i),xj1(:,i)]
CP2(:,i)=[xi2(:,i),xj2(:,i)]
CP3(:,i)=[xi3(:,i),xj3(:,i)]
CP4(:,i)=[xi4(:,i),xj4(:,i)]
%
numCity1=size(CP1(:,i),1);
numCity2=size(CP2(:,i),1);
numCity3=size(CP3(:,i),1);
numCity4=size(CP4(:,i),1);

% nilai tau
tau1=ones*[xi1(:,i),xj1(:,i)]
tau2=ones*[xi2(:,i),xj2(:,i)]
tau3=ones*[xi3(:,i),xj3(:,i)]
tau4=ones*[xi4(:,i),xj4(:,i)]

%=============================
%%
for i=1:n1
    for z=1:n1
        if CP1(i,z)==0
            h1(i,z)=0;
        else
            h1(i,z)=1/CP1(i,z);
        end
        if CP2(i,z)==0
            h2(i,z)=0;
        else
            h2(i,z)=1/CP2(i,z);
        end
        if CP3(i,z)==0
            h3(i,z)=0;
        else
            h3(i,z)=1/CP3(i,z);
        end
        if CP4(i,z)==0
            h4(i,z)=0;
        else
            h4(i,z)=1/CP4(i,z);
        end
            
    end
end

for i=1:MaxIt
    % Generate places for each ant    
    for i=1:m1
    mh1=h1;
    %1    
    for z=1:n1-1
        startplaces1(z,1)=fix(1+rand*(n1-1));
        c1=startplaces1(i,z);
        mh1(:,c1)=0;
        temp=(tau1(c1,:).^beta).*(mh1(c1,:).^alpha);
        s1=(sum(temp));
        P1=(1/s1).*temp;
        r=rand;
        s1=0;
            for u=1:n1
            s1=s1+P1(u);
            if r<=s1
                startplaces1(i,z+1)=u;
                break
            end
            end
    end
    end
    for i=1:m1
    mh2=h2;
    %2
    for z=1:n1-1
        startplaces2(z,1)=fix(1+rand*(n1-1));
        c2=startplaces2(i,z);
        mh2(:,c1)=0;
        temp2=(tau2(c2,:).^beta).*(mh2(c2,:).^alpha);
        s2=(sum(temp2));
        P2=(1/s2).*temp2;
        r=rand;
        s2=0;
            for u=1:n1
            s2=s2+P2(u);
            if r<=s2
                startplaces2(i,z+1)=u;
                break
            end
            end
    end
    end
    
    %3
    for i=1:m1
    mh3=h3;
    for z=1:n1-1
        startplaces3(z,1)=fix(1+rand*(n1-1));
        c3=startplaces3(i,z);
        mh3(:,c3)=0;
        temp3=(tau3(c3,:).^beta).*(mh3(c3,:).^alpha);
        s3=(sum(temp3));
        P3=(1/s3).*temp3;
        r=rand;
        s=0;
            for u=1:n1
            s3=s3+P3(u);
            if r<=s
                startplaces4(i,z+1)=u;
                break
            end
            end
            end
    end
    %4
    for i=1:m1
    mh4=h4;
    for z=1:n1-1
        startplaces4(z,1)=fix(1+rand*(n1-1));
        c4=startplaces4(i,z);
        mh4(:,c1)=0;
        temp4=(tau4(c4,:).^beta).*(mh4(c4,:).^alpha);
        s4=(sum4(temp));
        P4=(1/s4).*temp;
        r=rand;
        s4=0;
            for u=1:n1
            s4=s4+P4(u);
            if r<=s4
                startplaces4(i,z+1)=u;
                break
            end
            end
    end
    end 
    end 

    % Step 3: Calculate the cost --> total distace
    for i=1:m1
    s1=0;
    %1
    for z=1:n1
        s1=s1+CP1(i);
        f1(i)=s1;
        cost1=f1;
        f1=f1-rho*min(f1);
    end
    end
    %2
    for i=1:m1
    s2=0;
    for z=1:n1
        s2=s2+CP2(i);
        f2(i)=s2;
        cost2=f2;
        f2=f2-rho*min(f2);
    end
    end
    %3
    for i=1:m1
    s3=0;
    for z=1:n1
        s3=s3+CP3(i);
        f3(i)=s3;
        cost3=f3;
        f3=f3-rho*min(f3);
    end
    end
    %4
    for i=1:m1
    s4=0;
    for z=1:n1
        s4=s4+CP4(i);
        f4(i)=s4;
        cost4=f4;
        f4=f4-rho*min(f4);
    end
    end
   

 %% Update Race
    for i=1:m1
        %1
        for z=1:n1
        dt1=1/f1(i);
        tau1=(1-rho)*tau1+dt1;
        end
        %2
        for z=1:n1
        dt2=1/f2(i);
        tau2=(1-rho)*tau2+dt2;
        end
        %3
        for z=1:n1
        dt3=1/f3(i);
        tau3=(1-rho)*tau3+dt3;
        end
        %4
        for z=1:n1
        dt4=1/f4(i);
        tau4=(1-rho)*tau4+dt4;
        end
    end
   

   
    
   
    
