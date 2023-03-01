% This script will plot stress strain data for the standard dogbones

% Giordano Liska
% MAE 171A Solids Lab
% Created 2-24-2023

clear; clc;

RawPlotOn = false;
PlotOn = true;

%% Parameters
T = 0.0032; % Acrylic thickness [m]
W = 0.013; % Acrylic width [m]
L = 0.057; % Acrylic Length [m]
CrossA = T*W; % Cross sectional area [m^2]

%% Import Data

cs1.data = readmatrix('PMMA_CustomSpecimenFracture1_02242023.xlsx');
cs2.data = readmatrix('PMMA_CustomSpecimenFracture2_02242023.xlsx');
cs3.data = readmatrix('PMMA_CustomSpecimenFracture3_02242023.xlsx');

ssf1.data = readmatrix('PMMA_StandardDogboneFracture1_02102023.csv');
ssf2.data = readmatrix('PMMA_StandardDogboneFracture2_02102023.csv');

ssu.data = readmatrix('PMMA_StandardDogboneUnloading_02102023.csv');

if RawPlotOn == true

%% Plot Standard Dogbone Fracture Results

figure(1)
subplot(2,1,1)
plot(ssf1.data(:,2), ssf1.data(:,3), Linewidth = 1.25, Color = [.9 0 .1]);
xlabel('Displacement [mm]');
ylabel('Force [kN]');
axis([min(ssf1.data(:,2)) max(ssf1.data(:,2))+.25 min(ssf1.data(:,3)) max(ssf1.data(:,3))+.15]);

subplot(2,1,2)
plot(ssf2.data(:,2), ssf2.data(:,3), Linewidth = 1.25, Color = [0 .1 .9]);
xlabel('Displacement [mm]');
ylabel('Force [kN]');
axis([min(ssf2.data(:,2)) max(ssf2.data(:,2))+.25 min(ssf2.data(:,3)) max(ssf2.data(:,3))+.15]);

%% Plot Standard Dogbone Unloading Results

figure(2)
plot(ssu.data(:,2), ssu.data(:,3), Linewidth = 1.25, Color = [.9 0 .1]);
xlabel('Displacement [mm]');
ylabel('Force [kN]');
axis([min(ssu.data(:,2)) max(ssu.data(:,2))+.25 min(ssu.data(:,3)) max(ssu.data(:,3))+.15]);

%% Plot Custom Dogbone Fracture Results

figure(3)
subplot(3,1,1)
plot(cs1.data(:,2), cs1.data(:,3), Linewidth = 1.25, Color = [.9 0 .1]);
xlabel('Displacement [mm]');
ylabel('Force [kN]');
axis([min(cs1.data(:,2)) max(cs1.data(:,2))+.05 min(cs1.data(:,3)) max(cs1.data(:,3))+.05]);

subplot(3,1,2)
plot(cs2.data(:,2), cs2.data(:,3), Linewidth = 1.25, Color = [0 .1 .9]);
xlabel('Displacement [mm]');
ylabel('Force [kN]');
axis([min(cs2.data(:,2)) max(cs2.data(:,2))+.05 min(cs2.data(:,3)) max(cs2.data(:,3))+.05]);

subplot(3,1,3)
plot(cs3.data(:,2), cs3.data(:,3), Linewidth = 1.25, Color = [0 .55 .45]);
xlabel('Displacement [mm]');
ylabel('Force [kN]');
axis([min(cs3.data(:,2)) max(cs3.data(:,2))+.05 min(cs3.data(:,3)) max(cs3.data(:,3))+.05]);

end

%% Calculate Stress

cs1.EngStress = cs1.data(:,3)./(1000*CrossA);
cs2.EngStress = cs2.data(:,3)./(1000*CrossA);
cs3.EngStress = cs3.data(:,3)./(1000*CrossA);

ssf1.EngStress = ssf1.data(:,3)./(1000*CrossA);
ssf2.EngStress = ssf2.data(:,3)./(1000*CrossA);

ssu.EngStress = ssu.data(:,3)./(1000*CrossA);

%% Calculate Strain

cs1.EngStrain = cs1.data(:,2)./L;
cs2.EngStrain = cs2.data(:,2)./L;
cs3.EngStrain = cs3.data(:,2)./L;

ssf1.EngStrain = ssf1.data(:,2)./L;
ssf2.EngStrain = ssf2.data(:,2)./L;

ssu.EngStrain = ssu.data(:,2)./L;

%% Plot Stress-Strain Curves

if PlotOn == true

figure(4)
subplot(3,1,1)
plot(cs1.EngStrain, cs1.EngStress, Linewidth = 1.25, Color = [.9 0 .1]);
xlabel('Strain [unitless]');
ylabel('Stress [MPa]');
axis([min(cs1.EngStrain) max(cs1.EngStrain)+1 min(cs1.EngStress) max(cs1.EngStress)+.5]);

subplot(3,1,2)
plot(cs2.EngStrain, cs2.EngStress, Linewidth = 1.25, Color = [0 .1 .9]);
xlabel('Strain [unitless]');
ylabel('Stress [MPa]');
axis([min(cs2.EngStrain) max(cs2.EngStrain)+1 min(cs2.EngStress) max(cs2.EngStress)+.5]);

subplot(3,1,3)
plot(cs3.EngStrain, cs3.EngStress, Linewidth = 1.25, Color = [0 .55 .45]);
xlabel('Strain [unitless]');
ylabel('Stress [MPa]');
axis([min(cs3.EngStrain) max(cs3.EngStrain)+1 min(cs3.EngStress) max(cs3.EngStress)+.5]);

figure(5)


end




