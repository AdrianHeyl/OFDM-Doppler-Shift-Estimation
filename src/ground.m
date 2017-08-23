clc

%% VARIABLES

V = .15;
L = 2 * 1.38 + 2 * pi * .3;
T = L / V;
STEP = .1;

MIC1 = [0.55; -1];
MIC2 = [0; 0];

%% TIME STEPS

t1 = 0 : STEP : 1.38 / V;
t2a = 0 : STEP : pi * .15 / V;
t2 = t2a + t1(length(t1));
t3 = t1 + t2(length(t2));
t4 = t2a + t3(length(t3));
t = [t1, t2(2:length(t2)), t3(2:length(t3)), t4(2:length(t4))];

one1 = ones(1, length(t1));

%% POSITIONS

stp1_1 = [-.15 * cos(pi * t2a / t2a(end)) + .55; .15 * sin(pi * t2a / t2a(end)) + .69];
stp1_2 = [.7 * one1; .69 - V * t1];
stp1_3 = [.15 * cos(pi * t2a / t2a(end)) + .55; -.15 * sin(pi * t2a / t2a(end)) - .69];
stp1_4 = [.4 * one1; -.69 + V * t1];

%% DISTANCE AND VELOCITY

d1 = [sqrt((stp1_1(1, 1:length(stp1_1)) - MIC1(1)) .^ 2 + (stp1_1(2, 1:length(stp1_1)) - MIC1(2)) .^ 2), ...
    sqrt((stp1_2(1, 2:length(stp1_2)) - MIC1(1)) .^ 2 + (stp1_2(2, 2:length(stp1_2)) - MIC1(2)) .^ 2), ...
    sqrt((stp1_3(1, 2:length(stp1_3)) - MIC1(1)) .^ 2 + (stp1_3(2, 2:length(stp1_3)) - MIC1(2)) .^ 2), ...
    sqrt((stp1_4(1, 2:length(stp1_4)) - MIC1(1)) .^ 2 + (stp1_4(2, 2:length(stp1_4)) - MIC1(2)) .^ 2)];

dd1 = diff(d1) ./ STEP;

d2 = [sqrt((stp1_1(1, 1:length(stp1_1)) - MIC2(1)) .^ 2 + (stp1_1(2, 1:length(stp1_1)) - MIC2(2)) .^ 2), ...
    sqrt((stp1_2(1, 2:length(stp1_2)) - MIC2(1)) .^ 2 + (stp1_2(2, 2:length(stp1_2)) - MIC2(2)) .^ 2), ...
    sqrt((stp1_3(1, 2:length(stp1_3)) - MIC2(1)) .^ 2 + (stp1_3(2, 2:length(stp1_3)) - MIC2(2)) .^ 2), ...
    sqrt((stp1_4(1, 2:length(stp1_4)) - MIC2(1)) .^ 2 + (stp1_4(2, 2:length(stp1_4)) - MIC2(2)) .^ 2)];

dd2 = diff(d2) ./ STEP;

%% DISPLAY RESULTS

figure;
subplot(2, 3, [1,4])
plot(stp1_1(1, 1:length(stp1_1)), stp1_1(2, 1:length(stp1_1)), stp1_2(1, 1:length(stp1_2)), stp1_2(2, 1:length(stp1_2)), stp1_3(1, 1:length(stp1_3)), stp1_3(2, 1:length(stp1_3)), stp1_4(1, 1:length(stp1_4)), stp1_4(2, 1:length(stp1_4)))
title('Top-down view')
xlabel('x [m]')
ylabel('y [m]')
axis([0 1 -1 1])
pbaspect([1 2 1])

subplot(2, 3, 2)
plot(t, d1)
title(strcat('MIC1: (', strcat(num2str(MIC1(1)), strcat(';', strcat(num2str(MIC1(2)), ')')))))
xlabel('time [s]')
ylabel('distance [m]')

subplot(2, 3, 5)
plot(t, [0, dd1])
xlabel('time [s]')
ylabel('velocity [m/s]')

subplot(2, 3, 3)
plot(t, d2)
title(strcat('MIC2: (', strcat(num2str(MIC2(1)), strcat(';', strcat(num2str(MIC2(2)), ')')))))
xlabel('time [s]')
ylabel('distance [m]')

subplot(2, 3, 6)
plot(t, [0, dd2])
xlabel('time [s]')
ylabel('velocity [m/s]')