function [tv]=UP_BiP_complete
clc

tf=7200;             %Final time (s)
%--------Initial conditions-------%
C0=0.0502467;           %Cytosolic Ca2+ (uM)            
Cre0=506.395;          %ER Ca2+ (uM)
UP0=98; BiP0=66.5684; BUP0 =100-BiP0;PUP0=0.309452;
Ip0=0.03682;Pp0=0.0389;ATF60=0.0706;AB0=0.0367011;Ac0=0.2498;
[t,x]=ode23s(@UPRI2,[0:1:tf],[C0 Cre0 UP0 BiP0 BUP0 PUP0 Ip0 Pp0 ATF60 AB0 Ac0]);

tv=[t,x(:,4)];
ca=x(:,1);
cre=x(:,2);
up=x(:,3);
bip=x(:,4);
bup=x(:,5);
pup=x(:,6);
ip=x(:,7);
pb=x(:,8);
atf6=x(:,9);
ab=x(:,10);
ac=x(:,11);

figure(1)
hold on
plot(t/60,cre,'b','LineWidth',2)
xlabel('Time (min)','FontSize',15);ylabel('[Ca^{2+}]_{ER} ({\mu}M)','FontSize',15)
set(gca,'FontSize',15)

figure(2)
hold on
plot(t/60,up,'LineWidth',1.5)
xlabel('Time(min)','FontSize',15);ylabel('[UP] ({\mu}M)','FontSize',15)
set(gca,'FontSize',15)
box on

figure(3)
hold on
plot(t/60,ip/0.249,'b','LineWidth',2)
xlabel('Time (min)','FontSize',15);ylabel('[IRE1-p] ({\mu}M)','FontSize',15)
set(gca,'FontSize',15)

figure(4)
hold on
plot(t/60,(pb/(0.249)),'LineWidth',1.5)
xlabel('Time(min)','FontSize',15);ylabel('[PERK-p] ({\mu}M)','FontSize',15)
set(gca,'FontSize',15)
box on

figure(5)
hold on
plot(t/60,ac/0.2498,'LineWidth',1.5)
xlabel('Time(min)','FontSize',15);ylabel('[Ac] ({\mu}M)','FontSize',15)
set(gca,'FontSize',15)
box on

function dydt = UPRI2(t,y)

dydt=zeros(11,1);

gamma=10;

Ks=400;

%====tBu parameters======%
KI=2;
%====Other parameters====%
Ccyt=y(1);Cer=y(2);UP=y(3);BiP=y(4);BUP=y(5);
Ip=y(7);Pp=y(8);ATF6=y(9);AB=y(10);Ac=y(11);
na=1.09;k1=15;k2=0.7;ker=243;vmb=0.0804;vprod=1649;kd=0.087;tp=1772;
kout=0.0005;
Vs=0.02;
Ve=1.82;
Vp=5;
Ke=0.125;
Kp=1.5;

tBu=30;
tb=((KI^2)/(KI^2+tBu^2));

It=0.249;
KIp=0.7631;
KIB=51;
nib=5.38;
tIB=1579;
Pt=0.249;
KPp=1.2213;
KPB=14.35;
npb=1.3;
tPB=185.6;
kABplus=0.3025;kABminus=4.54;kAc=5.08E-04;kdAc=1.71E-05;
%-----------------------------------------------------------%
%--------------------------Expresiones----------------------%
%-----------------------------------------------------------%
%ATPase pumps
Jpmca=Vp*((Ccyt^2)/((Kp^2) + (Ccyt^2)));         %PMCA                                   
Jserca=Ve*((Ccyt^2)/((Ke^2) + (Ccyt^2)))*tb;     %SERCA
Jsoce=Vs*(Ks^4)/((Ks^4) + (Cer^(4)));           %SOCE
Jleak=kout*(Cer-Ccyt);                          %ER leak
%---------------------------------------------------------------------%
%-------------------------Mathematical Model--------------------------%
%---------------------------------------------------------------------%
dydt(1) = Jleak + Jsoce - Jpmca - Jserca;
dydt(2) = gamma*(Jserca-Jleak);
dydt(3) = vprod*y(6) - (vmb*(BiP)*(UP))- kd*UP + k2*BUP; %UP
dydt(4) = (k1 + k2)*(BUP) - (vmb*(BiP)*(UP));       %BiP
dydt(5) = (vmb*(BiP)*(UP)) - (k1 + k2)*(BUP); %BUP
dydt(6)=(((ker^na)/((Cer^na) + (ker^na)))-y(6))/tp; %PUP delay UP
dydt(7) = (It*KIp*((KIB^nib)/((BiP^nib)+(KIB^nib)))-Ip)/tIB;  %I-p
dydt(8) = (Pt*KPp*((KPB^npb)/((BiP^npb)+(KPB^npb)))-Pp)/tPB;  %P-p
dydt(9) = kABminus*AB-kABplus*ATF6*BiP-ATF6*kAc;            %ATF6
dydt(10) = kABplus*ATF6*BiP-kABminus*AB;                           %ATF6-BiP
dydt(11) = kAc*ATF6-kdAc*Ac;                                      %ATF6-c
end
end