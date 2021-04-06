
%  optimasi dengan ACO
MaxIt=100;                     % max numberof iteration
m1=1;
n1=1;
nAnt=40;                       % Number of Ant
alpha=0.72;                    % tetapan pengendali intensitas jejak semut
beta=0.66;                     % tetapan pengendali visibilitas
rho=0.75;                      % koefisien penguapan / evaporasi
n=1;
% ====================
% Batas-batas untuk 1
i=1;
a(:,i)=-205.8; 
b(:,i)=19.3; 
c(:,i)=244.4; 
% Batas-batas untuk 2
d(:,i)=19.3;
e(:,i)=244.4;
f(:,i)=469.5;    
% Batas-batas untuk 3
g(:,i)=244.4;
h(:,i)=469.5;
j(:,i)=694.7;
% Batas-batas untuk 4
k(:,i)=469.5;
l(:,i)=694.6;
m(:,i)=919.8;
% Batas-batas untuk 5
o(:,i)=694.6;
p(:,i)=919.8;
q(:,i)=1145;

r(:,i)=919.9;
s(:,i)=1145;
t(:,i)=13470;

v(:,i)=1145;
w(:,i)=1370;
x(:,i)=1595;
% ====================

% Inisialisasi posisi awal
x1(:,i)=(c(:,i)-(rand(n,1)*(c(:,i)-b(:,i)))); %constraint: b<=x1<=c
x2(:,i)=(d(:,i)-(rand(n,1)*(d(:,i)-a(:,i)))); %constraint: a<=x2<=d
x3(:,i)=(f(:,i)-(rand(n,1)*(f(:,i)-e(:,i)))); %constraint: e<=x3<=f
x4(:,i)=(h(:,i)-(rand(n,1)*(h(:,i)-g(:,i)))); %constraint: g<=x4<=h
x5(:,i)=(j(:,i)-(rand(n,1)*(j(:,i)-h(:,i)))); %constraint: h<=x5<=j
x6(:,i)=(l(:,i)-(rand(n,1)*(l(:,i)-k(:,i)))); %constraint: k<=x6<=l
x7(:,i)=(m(:,i)-(rand(n,1)*(m(:,i)-l(:,i)))); %constraint: l<=x7<=m
x8(:,i)=(p(:,i)-(rand(n,1)*(p(:,i)-o(:,i)))); %constraint: o<=x8<=p
x9(:,i)=(q(:,i)-(rand(n,1)*(q(:,i)-p(:,i)))); %constraint: p<=x7<=q
x10(:,i)=(s(:,i)-(rand(n,1)*(s(:,i)-r(:,i)))); %constraint: r<=x8<=s
x11(:,i)=(w(:,i)-(rand(n,1)*(w(:,i)-v(:,i)))); %constraint: v<=x7<=w
x12(:,i)=(t(:,i)-(rand(n,1)*(t(:,i)-s(:,i)))); %constraint: s<=x8<=t
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
myu9(:,i)=(j(:,i)-x5(:,i))/(j(:,i)-h(:,i)); % myu9=derajat keanggotaan 3 batas kanan (j) 
y7=0;                                       % y7=derajat keanggotaan ideal 3 batas kiri (g) adalah 0
y8=1;                                       % y8=derajat keanggotaan ideal 3 batas tengah (h) adalah 1 (tetap)
y9=0;                                       % y9=derajat keanggotaan ideal 3 batas kanan (j) adalah 0
 
myu10(:,i)=(x6(:,i)-k(:,i))/(l(:,i)-k(:,i)); % myu10=derajat keanggotaan 4 batas kiri (k)
myu11(:,i)=1;                                % myu11=derajat keanggotaan 4 batas tengah (l) (tetap)
myu12(:,i)=(m(:,i)-x7(:,i))/(m(:,i)-l(:,i)); % myu12=derajat keanggotaan 4 batas kanan (m) 
y10=0;                                       % y10=derajat keanggotaan ideal 4 batas kiri (k) adalah 0
y11=1;                                       % y11=derajat keanggotaan ideal 4 batas tengah (l) adalah 1 (tetap) 
y12=0;                                       % y12=derajat keanggotaan ideal P4 batas kanan (m) adalah 0
 
myu13(:,i)=(x8(:,i)-o(:,i))/(p(:,i)-o(:,i)); % myu13=derajat keanggotaan 5 batas kiri (o)
myu14(:,i)=1;                                % myu14=derajat keanggotaan 5 batas tengah (p) (tetap)
myu15(:,i)=(q(:,i)-x9(:,i))/(q(:,i)-p(:,i)); % myu15=derajat keanggotaan 5 batas kanan (q) (tetap) 
y13=0;                                       % y13=derajat keanggotaan ideal 5 batas kiri (o) adalah 0
y14=1;                                       % y14=derajat keanggotaan ideal 5 batas tengah (p) adalah 1 (tetap)
y15=0;                                       % y15=derajat keanggotaan ideal 5 batas kanan (q) adalah 1 (tetap)

myu16(:,i)=(x10(:,i)-r(:,i))/(s(:,i)-r(:,i)); % myu13=derajat keanggotaan 5 batas kiri (o)
myu17(:,i)=1;                                % myu14=derajat keanggotaan 5 batas tengah (p) (tetap)
myu18(:,i)=(t(:,i)-x11(:,i))/(q(:,i)-s(:,i));;                                % myu15=derajat keanggotaan 5 batas kanan (q) (tetap) 
y16=0;                                       % y13=derajat keanggotaan ideal 5 batas kiri (o) adalah 0
y17=1;                                       % y14=derajat keanggotaan ideal 5 batas tengah (p) adalah 1 (tetap)
y18=0;                                       % y15=derajat keanggotaan ideal 5 batas kanan (q) adalah 1 (tetap)

myu19(:,i)=(x12(:,i)-v(:,i))/(w(:,i)-v(:,i)); % myu13=derajat keanggotaan 5 batas kiri (o)
myu20(:,i)=1;                                % myu14=derajat keanggotaan 5 batas tengah (p) (tetap)
myu21(:,i)=1;                                % myu15=derajat keanggotaan 5 batas kanan (q) (tetap) 
y19=0;                                       % y13=derajat keanggotaan ideal 5 batas kiri (o) adalah 0
y20=1;                                       % y14=derajat keanggotaan ideal 5 batas tengah (p) adalah 1 (tetap)
y21=1;                                       % y15=derajat keanggotaan ideal 5 batas kanan (q) adalah 1 (tetap)

mse(:,i)=(((y1-myu1(:,i)).^2)+((y2-myu2(:,i)).^2)+((y3-myu3(:,i)).^2)+((y4-myu4(:,i)).^2)+((y5-myu5(:,i)).^2)+...
    ((y6-myu6(:,i)).^2)+((y7-myu7(:,i)).^2)+((y8-myu8(:,i)).^2)+((y9-myu9(:,i)).^2)+((y10-myu10(:,i)).^2)+((y11-myu11(:,i)).^2)+...
    ((y12-myu12(:,i)).^2)+((y13-myu13(:,i)).^2)+((y14-myu14(:,i)).^2)+((y15-myu15(:,i)).^2)+((y16-myu16(:,i)).^2)+((y17-myu17(:,i)).^2)+((y18-myu18(:,i)).^2)+((y19-myu19(:,i)).^2)+((y20-myu20(:,i)).^2)+((y21-myu21(:,i)).^2))/21;

%====================================
%antPosition(AP)
AP1=x1;
AP2=x2;
AP3=x3;
AP4=x4;
AP5=x5;
AP6=x6;
AP7=x7;
AP8=x8;
AP9=x9;
AP10=x10;
AP11=x11;
AP12=x12;
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
%x5
xi5(:,i)=AP5(:,i)+(mse(:,i)*rho(:,i));
xj5(:,i)=AP5(:,i)-(mse(:,i)*rho(:,i));
%x6
xi6(:,i)=AP6(:,i)+(mse(:,i)*rho(:,i));
xj6(:,i)=AP6(:,i)-(mse(:,i)*rho(:,i));
%x7
xi7(:,i)=AP7(:,i)+(mse(:,i)*rho(:,i));
xj7(:,i)=AP7(:,i)-(mse(:,i)*rho(:,i));
%x8
xi8(:,i)=AP8(:,i)+(mse(:,i)*rho(:,i));
xj8(:,i)=AP8(:,i)-(mse(:,i)*rho(:,i));
%x9
xi9(:,i)=AP9(:,i)+(mse(:,i)*rho(:,i));
xj9(:,i)=AP9(:,i)-(mse(:,i)*rho(:,i));
%x10
xi10(:,i)=AP10(:,i)+(mse(:,i)*rho(:,i));
xj10(:,i)=AP10(:,i)-(mse(:,i)*rho(:,i));
%x11
xi11(:,i)=AP11(:,i)+(mse(:,i)*rho(:,i));
xj11(:,i)=AP11(:,i)-(mse(:,i)*rho(:,i));
%x12
xi12(:,i)=AP12(:,i)+(mse(:,i)*rho(:,i));
xj12(:,i)=AP12(:,i)-(mse(:,i)*rho(:,i));

% city position(CP)
CP1(:,i)=[xi1(:,i),xj1(:,i)]
CP2(:,i)=[xi2(:,i),xj2(:,i)]
CP3(:,i)=[xi3(:,i),xj3(:,i)]
CP4(:,i)=[xi4(:,i),xj4(:,i)]
CP5(:,i)=[xi5(:,i),xj5(:,i)]
CP6(:,i)=[xi6(:,i),xj6(:,i)]
CP7(:,i)=[xi7(:,i),xj7(:,i)]
CP8(:,i)=[xi8(:,i),xj8(:,i)]
CP9(:,i)=[xi9(:,i),xj9(:,i)]
CP10(:,i)=[xi10(:,i),xj10(:,i)]
CP11(:,i)=[xi11(:,i),xj11(:,i)]
CP12(:,i)=[xi12(:,i),xj12(:,i)]

% nilai tau atau pheromone
tau1=ones*[xi1(:,i),xj1(:,i)]
tau2=ones*[xi2(:,i),xj2(:,i)]
tau3=ones*[xi3(:,i),xj3(:,i)]
tau4=ones*[xi4(:,i),xj4(:,i)]
tau5=ones*[xi5(:,i),xj5(:,i)]
tau6=ones*[xi6(:,i),xj6(:,i)]
tau7=ones*[xi7(:,i),xj7(:,i)]
tau8=ones*[xi8(:,i),xj8(:,i)]
tau9=ones*[xi9(:,i),xj9(:,i)]
tau10=ones*[xi10(:,i),xj10(:,i)]
tau11=ones*[xi11(:,i),xj11(:,i)]
tau12=ones*[xi12(:,i),xj12(:,i)]

%=============================
%%  jalur semut
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
        if CP5(i,z)==0
            h5(i,z)=0;
        else
            h5(i,z)=1/CP5(i,z);
        end
        if CP6(i,z)==0
            h6(i,z)=0;
        else
            h6(i,z)=1/CP6(i,z);
        end
        if CP7(i,z)==0
            h7(i,z)=0;
        else
            h7(i,z)=1/CP7(i,z);
        end
        if CP8(i,z)==0
            h8(i,z)=0;
        else
            h8(i,z)=1/CP8(i,z);
        end  
        if CP9(i,z)==0
            h9(i,z)=0;
        else
            h9(i,z)=1/CP9(i,z);
        end
        if CP10(i,z)==0
            h10(i,z)=0;
        else
            h10(i,z)=1/CP10(i,z);
        end
        if CP11(i,z)==0
            h11(i,z)=0;
        else
            h11(i,z)=1/CP11(i,z);
        end
        if CP12(i,z)==0
            h12(i,z)=0;
        else
            h12(i,z)=1/CP12(i,z);
        end
    end
end

for i=1:MaxIt
    % Posisi tiap semut     
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
    %5
    for i=1:m1
    mh5=h5;
    for z=1:n1-1
        startplaces5(z,1)=fix(1+rand*(n1-1));
        c5=startplaces5(i,z);
        mh5(:,c5)=0;
        temp5=(tau5(c5,:).^beta).*(mh5(c5,:).^alpha);
        s5=(sum(temp5));
        P5=(1/s5).*temp5;
        r=rand;
        s5=0;
            for u=1:n1
            s5=s5+P1(u);
            if r<=s5
                startplaces5(i,z+1)=u;
                break
            end
            end
    end 
    end
    %6
    for i=1:m1
    mh6=h6;
    for z=1:n1-1
        startplaces6(z,1)=fix(1+rand*(n1-1));
        c6=startplaces6(i,z);
        mh6(:,c6)=0;
        temp6=(tau6(c6,:).^beta).*(mh6(c6,:).^alpha);
        s6=(sum(temp));
        P6=(1/s6).*temp;
        r=rand;
        s6=0;
            for u=1:n1
            s6=s6+P1(u);
            if r<=s6
                startplaces6(i,z+1)=u;
                break
            end
            end
    end
    end
    %7
    for i=1:m1
    mh7=h7;
    for z=1:n1-1
        startplaces7(z,1)=fix(1+rand*(n1-1));
        c7=startplaces7(i,z);
        mh7(:,c7)=0;
        temp7=(tau7(c7,:).^beta).*(mh7(c7,:).^alpha);
        s7=(sum(temp7));
        P7=(1/s7).*temp7;
        r=rand;
        s7=0;
            for u=1:n1
            s7=s7+P7(u);
            if r<=s7
                startplaces7(i,z+1)=u;
                break
            end
            end
    end
    end
    %8
    for i=1:m1
    mh8=h8;
    for z=1:n1-1
        startplaces8(z,1)=fix(1+rand*(n1-1));
        c8=startplaces8(i,z);
        mh8(:,c8)=0;
        temp8=(tau8(c8,:).^beta).*(mh8(c8,:).^alpha);
        s8=(sum(temp8));
        P8=(1/s8).*temp8;
        r=rand;
        s=0;
            for u=1:n1
            s8=s8+P8(u);
            if r<=s8
                startplaces8(i,z+1)=u;
                break
            end
            end
    end 
    end
    %9
    for i=1:m1
    mh9=h9;
    for z=1:n1-1
        startplaces9(z,1)=fix(1+rand*(n1-1));
        c7=startplaces9(i,z);
        mh9(:,c9)=0;
        temp9=(tau9(c9,:).^beta).*(mh9(c9,:).^alpha);
        s9=(sum(temp9));
        P9=(1/s9).*temp9;
        r=rand;
        s9=0;
            for u=1:n1
            s9=s9+P9(u);
            if r<=s9
                startplaces7(i,z+1)=u;
                break
            end
            end
    end
    end
    %10
    for i=1:m1
    mh10=h10;
    for z=1:n1-1
        startplaces10(z,1)=fix(1+rand*(n1-1));
        c10=startplaces10(i,z);
        mh10(:,c10)=0;
        temp10=(tau10(c10,:).^beta).*(mh10(c10,:).^alpha);
        s10=(sum(temp10));
        P10=(1/s10).*temp10;
        r=rand;
        s10=0;
            for u=1:n1
            s10=s10+P10(u);
            if r<=s10
                startplaces10(i,z+1)=u;
                break
            end
            end
    end
    end
    %11
    for i=1:m1
    mh11=h11;
    for z=1:n1-1
        startplaces11(z,1)=fix(1+rand*(n1-1));
        c11=startplaces11(i,z);
        mh11(:,c11)=0;
        temp11=(tau11(c11,:).^beta).*(mh11(c11,:).^alpha);
        s11=(sum(temp11));
        P11=(1/s11).*temp11;
        r=rand;
        s11=0;
            for u=1:n1
            s11=s11+P11(u);
            if r<=s11
                startplaces11(i,z+1)=u;
                break
            end
            end
    end
    end
    %12
    for i=1:m1
    mh12=h12;
    for z=1:n1-1
        startplaces12(z,1)=fix(1+rand*(n1-1));
        c12=startplaces12(i,z);
        mh12(:,c12)=0;
        temp12=(tau12(c12,:).^beta).*(mh12(c12,:).^alpha);
        s12=(sum(temp12));
        P12=(1/s12).*temp12;
        r=rand;
        s12=0;
            for u=1:n1
            s12=s12+P12(u);
            if r<=s12
                startplaces12(i,z+1)=u;
                break
            end
            end
    end
    end
end   
    % Menghitung  distance
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
    %5
    for i=1:m1
    s5=0;
    for z=1:n1
        s5=s5+CP5(i);
        f5(i)=s5;
        cost5=f5;
        f5=f5-rho*min(f5);
    end
    end
    %6
    for i=1:m1
    s6=0;
    for z=1:n1
        s6=s6+CP6(i);
        f6(i)=s6;
        cost6=f6;
        f6=f6-rho*min(f6);
    end
    end
     %7
    for i=1:m1
    s7=0;
    for z=1:n1
        s7=s7+CP7(i);
        f7(i)=s7;
        cost7=f7;
        f7=f7-rho*min(f7);
    end
    end
     %8
    for i=1:m1
    s8=0;
    for z=1:n1
        s8=s8+CP8(i);
        f8(i)=s8;
        cost8=f8;
        f8=f8-rho*min(f8);
    end
    end
    %9
    for i=1:m1
    s9=0;
    for z=1:n1
        s9=s9+CP9(i);
        f9(i)=s9;
        cost9=f9;
        f9=f9-rho*min(f9);
    end
    end
    %10
    for i=1:m1
    s10=0;
    for z=1:n1
        s10=s10+CP10(i);
        f10(i)=s10;
        cost10=f10;
        f10=f10-rho*min(f10);
    end
    end
    %11
    for i=1:m1
    s11=0;
    for z=1:n1
        s11=s11+CP11(i);
        f11(i)=s11;
        cost11=f11;
        f11=f11-rho*min(f11);
    end
    end
    %12
    for i=1:m1
    s12=0;
    for z=1:n1
        s12=s12+CP12(i);
        f12(i)=s12;
        cost12=f12;
        f12=f12-rho*min(f12);
    end
    end
    
    

 %% Update Jalur semut dan pheromone
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
        %5
        for z=1:n1
        dt5=1/f5(i);
        tau5=(1-rho)*tau5+dt5;
        end
        %6
        for z=1:n1
        dt6=1/f6(i);
        tau6=(1-rho)*tau6+dt6;
        end
        %7
        for z=1:n1
        dt7=1/f7(i);
        tau7=(1-rho)*tau7+dt7;
        end
        %8
        for z=1:n1
        dt8=1/f8(i);
        tau8=(1-rho)*tau8+dt8;
        end
        %9
        for z=1:n1
        dt9=1/f9(i);
        tau9=(1-rho)*tau9+dt9;
        end
        %10
        for z=1:n1
        dt10=1/f10(i);
        tau10=(1-rho)*tau10+dt10;
        end
        %11
        for z=1:n1
        dt11=1/f11(i);
        tau11=(1-rho)*tau11+dt11;
        end
        %12
        for z=1:n1
        dt12=1/f12(i);
        tau12=(1-rho)*tau12+dt12;
        end
    end
   

    
   
    
