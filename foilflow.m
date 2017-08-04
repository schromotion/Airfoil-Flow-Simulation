%--------------------------------------------------------------------------
%                        AEM 313: Aerodynamics 1
%                               Project #2
%   Lift and Drag Calculation for a Conventional Wing-Tail Configuration
%
%                Kellen Schroeter & Lineu Ghaiasso Parra
%--------------------------------------------------------------------------

clc
clear all
close all

%                ***********Open Text Files*************
%Read the 'airfoil' text file
geometry = load('airfoil.txt');
n = size(geometry);
pts = n(1);
x = geometry(:,1);
zup = geometry(:,2);
zlo = geometry(:,3);

%Read the 'wing' text file
wing = load('wing.txt');
Cr_wing = wing(1,1);
b_wing_ratio = wing(2,1);
tr_wing = wing(3,1);
sweep_wing = wing(4,1);
thick_wing_ratio= wing(5,1);
cg = wing(6,1);

%Read the 'tail' text file
tail = load('tail.txt');
Cr_tail_ratio = tail(1,1);
b_tail_ratio = tail(2,1);
tr_tail = tail(3,1);
sweep_tail = tail(4,1);
thick_tail_ratio = tail(5,1);
offset_tail = tail(6,1);
incidence_tail = tail(7,1);

%Read the 'conditions' text file
conditions = load('conditions.txt');
U = conditions(1,1);
rho = conditions(2,1);
P = conditions(3,1);

%Read the 'tau.txt' text file
tau = load('tau.txt');

%Read the 'e.txt' text file
e = load('e.txt');




%                ***********Mach Number************* 

%Perfect gas equation
R  = 287.05;
temp = P / (rho * R);

%Speed of sound
gamma = 1.4;
a = sqrt(gamma * R * temp);

%Mach number
M = U/a;




%                ***********Aspect Ratios************* 

%Tip Chord
Ct_wing = Cr_wing * tr_wing;
Cr_tail = Cr_wing * Cr_tail_ratio;
Ct_tail = Cr_tail * tr_tail;

%Dimensional Wingspan
b_wing = Cr_wing * b_wing_ratio;
b_tail = Cr_tail * b_tail_ratio;

%Area
S_wing = b_wing * (Cr_wing + Ct_wing) / 2;
S_tail = b_tail * (Cr_tail + Ct_tail) / 2;

%Aspect Ratio
AR_wing = b_wing * b_wing / S_wing;
AR_tail = b_tail * b_tail / S_tail;




%                    ***********MACS*************

%MAC for trapezoidal wing
MAC_wing = (2/3) * Cr_wing * (((tr_wing^2) + tr_wing + 1) / (tr_wing + 1));
MAC_tail = (2/3) * Cr_tail * (((tr_tail^2) + tr_tail + 1) / (tr_tail + 1));




%              ***********Reynolds Numbers*************

%Sutherland's Equation
Ci = 1.458*10^-6;
Cii = 110.4;
mu = Ci * (temp^1.5) / (temp + Cii);

%Reynolds Number
Re_wing = rho * U * MAC_wing / mu;
Re_tail = rho * U * MAC_tail / mu;




%           ***********Zero Lift Angle of Attack************

%Calculate the Mean Camber Line
MCL = ones(pts,1);
for i=1:1:pts;
    MCL(i,1) = (zup(i) + zlo(i))/2;
end

%Integrate for the Zero Lift AOA
alpha0L = 0;
A0 = 0;
A1 = 0;
A2 = 0;
for i=2:1:pts;
    dzdx = (MCL(i,1)-MCL(i-1,1))/(x(i)-x(i-1));
    costheta2 = 1 - 2*x(i);
    costheta1 = 1 - 2*x(i-1);
    dtheta = acos(costheta2) - acos(costheta1);
    sweep_tip_tail = (-1/pi)*dzdx*(1/2)*dtheta*((costheta2-1)+(costheta1-1));
    alpha0L = sweep_tip_tail + alpha0L;
    
    aux2 = -(1/pi)*dzdx*dtheta;
    A0 = aux2 + A0;
    
    aux3 = (2/pi)*dzdx*(1/2)*dtheta*(costheta1+costheta2);
    A1 = aux3 + A1;
    
    aux4 = (2/pi)*dzdx*(1/2)*dtheta*((-1+2*costheta1^2)+(-1+2*costheta2^2));
    A2 = aux4 +A2;
end




%               ***********Lift Curve Slope************

%Glauert Correction Factor
aspect_ratio = [ 16 12 08 04 02 ];
taper_ratio = [ 1.0 0.8 0.6 0.4 0.2 0.0 ];
tau1 = interp2 ( aspect_ratio , taper_ratio , tau , AR_wing , tr_wing );
tau2 = interp2 ( aspect_ratio , taper_ratio , tau , AR_tail , tr_tail );

%Lift Curve Slope
a_o = 2 * pi;
lc_slope_wing = (a_o) / (1 + (a_o * (1 + tau1) / pi / AR_wing));
lc_slope_tail = (a_o) / (1 + (a_o * (1 + tau2) / pi / AR_tail));




%                 ***********Cl vs. Cd************

%Lift Coefficient
alpha = -6:.5:16;
Cl_wing = lc_slope_wing * (pi/180) * (alpha - (180/pi)*alpha0L);
Cl_tail = lc_slope_tail * (pi/180) * (alpha + incidence_tail);
CL = Cl_wing + Cl_tail * (S_tail / S_wing);

%Vortex Drag Coefficient
e1 = interp2 ( aspect_ratio , taper_ratio , e , AR_wing , tr_wing );
e2 = interp2 ( aspect_ratio , taper_ratio , e , AR_tail , tr_tail);
Cdv_wing = Cl_wing.^2 / pi / AR_wing / e1;
Cdv_tail = Cl_tail.^2 / pi / AR_tail / e2;

%Form Drag Coefficient Wing
K1 = (5.46 * (thick_wing_ratio)^2) + ((1.55 - sin(sweep_wing * pi /180)) * (thick_wing_ratio)) + 1;
if Re_wing <= 500000
    Cf1 = 1.328 / (Re_wing^(.5));
end
if Re_wing > 500000
    Cf1 = (0.455 / (log10(Re_wing)^(2.58))) - (1700 / Re_wing);
end
S_wet_wing = 2 * (1 + .2 * thick_wing_ratio) * S_wing;
Cdo_wing = K1 * Cf1 * S_wet_wing / S_wing;

%Form Drag Coefficient Tail
K2 = (5.46 * (thick_tail_ratio)^2) + ((1.55 - sin(sweep_tail * pi /180)) * (thick_tail_ratio)) + 1;
if Re_tail <= 500000
    Cf2 = 1.328 / (Re_tail^(.5));
end
if Re_tail > 500000
    Cf2 = (0.455 / (log10(Re_tail)^(2.58))) - (1700 / Re_tail);
end
S_wet_tail = 2 * (1 + .2 * thick_tail_ratio) * S_tail;
Cdo_tail = K2 * Cf2 * S_wet_tail / S_wing;

%Total Drag Coefficient
Cdo = Cdo_wing + Cdo_tail;
Cdv = Cdv_wing + Cdv_tail*(S_tail/S_wing);
Cd = Cdv + Cdo;



%               ***********Cm,cg vs. AOA************

%Cmc4_wing = (pi/4)*(A2-A1);
%Cmc4_tail = Cl_tail/4;
%Cmcg_tail_cg = -pi/2*(pi/180)*alpha - x_tail * Cl_tail * (S_tail/S_wing);

x_wing = (1/4)-cg;
x_tail = offset_tail + x_wing - (1/4) + (1/4)*Cr_tail;
Cmcg_wing_cg = -Cl_wing*x_wing + Cl_wing/4 - (pi/2)*(A0+A1-(A2/2));
Cmcg_tail_cg = -Cl_tail/4 - x_tail * Cl_tail * (S_tail/S_wing);
Cmcg = Cmcg_wing_cg + Cmcg_tail_cg;




%                ***********L/D vs. AOA************

%Lift-to-drag ratio
lift_drag_ratio = CL./Cd;




%                ********* Draw wing and tail *******
%Tail
Tip_length_tail = Cr_tail * tr_tail;
semispan_tail = ( b_tail_ratio * Cr_tail ) / 2;
sweep_tip_tail = - tan(sweep_tail * pi/180) * semispan_tail;
l = offset_tail*Cr_wing*Cr_tail_ratio;
D1 = [-semispan_tail
    0
    semispan_tail
    semispan_tail
    0
    -semispan_tail
    -semispan_tail];
D2 = [sweep_tip_tail+(1/4)*Tip_length_tail-l
    (1/4)*Cr_tail-l
    sweep_tip_tail+(1/4)*Tip_length_tail-l
    sweep_tip_tail-(3/4)*Tip_length_tail-l
    -(3/4)*Cr_tail-l
    sweep_tip_tail-(3/4)*Tip_length_tail-l
    sweep_tip_tail+(1/4)*Tip_length_tail-l];

%Wing
Tip_length_wing = Cr_wing * tr_wing;
semispan_wing = ( b_wing_ratio * Cr_wing ) / 2;
sweep_tip_wing = - tan(sweep_wing*pi/180) * semispan_wing;
D3 = [-semispan_wing
    0
    semispan_wing
    semispan_wing
    0
    -semispan_wing
    -semispan_wing];
D4 = [sweep_tip_wing+(1/4)*Tip_length_wing 
    (1/4)*Cr_wing
    sweep_tip_wing+(1/4)*Tip_length_wing
    sweep_tip_wing-(3/4)*Tip_length_wing
    -(3/4)*Cr_wing
    sweep_tip_wing-(3/4)*Tip_length_wing
    sweep_tip_wing+(1/4)*Tip_length_wing];




%                ********** Warning Messages*********

%Text File Errors
b=false;
if sweep_wing > 30
    disp('**Warning: The wing sweep is beyond the sweep limit (Sweep > 30)')
    b=true;
end
if sweep_tail > 30
    disp('**Warning: The tail sweep is beyond the sweep limit (Sweep > 30)')
    b=true;
end
if M > .3
    disp ('**Warning: This mach number is beyond the compressible limit (M>.3)')
    b=true;
end
if AR_wing < 4
    disp ('**Warning: The Aspect Ratio of the wing is less than the AR limit (AR<4)')
    b=true;
end
if AR_tail < 4
    disp ('**Warning: The Aspect Ratio of the tail is less than the AR limit (AR<4)')
    b=true;
end
if Re_wing < 200000
    disp('**Warning: The Reynolds Number for the wing is less than the Re limit (Re<200000)')
    b=true;
end
if Re_tail < 200000
    disp('**Warning: The Reynolds Number for the tail is less than the Re limit (Re<200000)')
    b=true;
end
if b==true
    fprintf('\n\n')
end



%              *********** Displays ***************

fprintf('Mach Number = %5.3f\n', M)
fprintf('Alpha(L=0) = %4.2f degrees ---> TAT (Thin Airfoil Theory)\n', (180/pi)*alpha0L)
fprintf('Lift Curve Slope = %4.4f degrees^-1\n\n', (pi/180)*lc_slope_wing)


%Wing
fprintf('Wing:\n')
fprintf('AR = %4.2f\n', AR_wing)
fprintf('MAC = %4.2f meters\n', MAC_wing)
fprintf('Re = %8.2e\n\n', Re_wing)

%Wing
fprintf('Tail:\n')
fprintf('AR = %4.2f\n', AR_tail)
fprintf('MAC = %4.2f meters\n', MAC_tail)
fprintf('Re = %8.2e\n', Re_tail)




%              *********** Plots ***************

figure('units','normalized','outerposition',[0 0 1 1])
%Plot The Airfoil Shape
subplot(2,3,1);
plot( x , zup , 'b-' , x , zlo , 'b-' )
title('Airfoil Shape')
xlabel('x/c')
ylabel('z/c')
grid on 
box on
axis equal
axis([0 1 -0.2 0.2])

%Plot Cl vs. Cd
subplot(2,3,3)
plot( Cd , CL , 'b-')
title('Cl vs. Cd')
xlabel('Cd')
ylabel('Cl')
grid on 
box on

%Plot Cmcg vs. AoA
subplot(2,3,4)
hold on
plot(alpha, Cmcg )
title('Cmcg vs. AoA')
xlabel('AoA (Degrees)')
ylabel('Cmcg')
grid on 
box on
leg = legend('Cmcg wing and tail');
set(leg, 'Position', [ 0.230    0.36    0.0918    0.0694])

%Plot L/D vs. AoA
subplot(2,3,5)
plot( alpha, lift_drag_ratio , 'b-')
title('L/D vs. AoA')
xlabel('AoA (Degrees)')
ylabel('L/D')
grid on 
box on

%Plot wing and tail
subplot(2,3,6)
hold on
axis equal
grid on
box on
title('Wing, Tail and CG location')
plot(0 , (Cr_wing/4) - cg * Cr_wing , 'x', 'MarkerSize', 12);
axis([-semispan_wing-1 semispan_wing+1 -(3/4)*Cr_tail-l-1 (1/4)*Cr_wing+1]);
line(D1, D2);
line(D3, D4);

%Plot Cl wing and tail
sub2 = subplot(2,3,2);
title('Cl vs. AoA')
xlabel('AoA (Degrees)')
ylabel('Cl')
hold on
grid on
box on
plot(alpha, Cl_wing, 'r')
plot(alpha, Cl_tail, 'g')
plot(alpha, CL)
leg = legend('Cl wing','Cl tail', 'Cl wing and tail');
get2 = get(sub2, 'Position');
set(leg, 'Position', [ get2(1)+(1/20)*get2(3) get2(2)+(3/4)*get2(4) 0.0918 0.0694])