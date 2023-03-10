% This script will plot stress strain data for the specimens
% Giordano Liska
% MAE 171A Solids Lab
% Created 2-24-2023

clear; clc;

RawPlotOn = false;
PlotOn = true;

%% Parameters
T = 0.0032; % Acrylic thickness [m]
W = 0.013; % Acrylic width [m]
L = 0.057*1000; % Acrylic Length [mm]
CrossA = T*W; % Cross sectional area [m^2]
DR = 0.4064; % Material density ratio [unitless]

%% Import Data


cs1.data = readmatrix('PMMA_CustomSpecimenFracture1_02242023.xlsx');
for i = 1:7
    filename = sprintf('cs1-DICe_%i.txt', i);
    dat = readmatrix(filename);
    cs1.dice(:,:,i) = dat;
end

cs2.data = readmatrix('PMMA_CustomSpecimenFracture2_02242023.xlsx');
for i = 1:7
    filename = sprintf('cs2-DICe_%i.txt', i);
    dat = readmatrix(filename);
    cs2.dice(:,:,i) = dat;
end

cs3.data = readmatrix('PMMA_CustomSpecimenFracture3_02242023.xlsx');
for i = 1:7
    filename = sprintf('cs3-DICe_%i.txt', i);
    dat = readmatrix(filename);
    cs3.dice(:,:,i) = dat;
end

ssf1.data = readmatrix('PMMA_StandardDogboneFracture1_02102023.csv');
for i = 1:18
    filename = sprintf('DICe_solution_%i', i+2267);
    dat = readmatrix(filename);
    ssf1.dice(:,:,i) = dat;
end

ssf2.data = readmatrix('PMMA_StandardDogboneFracture2_02102023.csv');
for i = 1:20
    filename = sprintf('DICe_solution_%i', i+2295);
    dat = readmatrix(filename);
    ssf2.dice(:,:,i) = dat;
end

ssu.data = readmatrix('PMMA_StandardDogboneUnloading_02102023.csv');

%% Calculate Stress

cs1.EngStress = cs1.data(:,3)./(DR*1000*CrossA); % Engineering stress for custom specimen 1
cs2.EngStress = cs2.data(:,3)./(DR*1000*CrossA); % Engineering stress for custom specimen 2
cs3.EngStress = cs3.data(:,3)./(DR*1000*CrossA); % Engineering stress for custom specimen 3

cs_TS = [max(cs1.EngStress), max(cs2.EngStress), max(cs3.EngStress)]; % Max stress for custom specimens

cs_avgTS = mean(cs_TS); % Average max stress for custom specimens

cs_TS_SD = std(cs_TS); % Standard deviation of max stress for custom specimens

ssf1.EngStress = ssf1.data(:,3)./(1000*CrossA); % Engineering stress for bulk specimen 1
ssf2.EngStress = ssf2.data(:,3)./(1000*CrossA); % Engineering stress for bulk specimen 2

ssf_TS = [max(ssf1.EngStress), max(ssf2.EngStress)]; % Max stress for bulk specimens

ssf_avgTS = mean(ssf_TS); % Average max stress for bulk specimens

ssf_TS_SD = std(ssf_TS); % Standard deviation of max stress for bulk specimens

s_TSerr = [cs_TS_SD/2 ssf_TS_SD/2]; % Vector of one standard deviation for both specimens

avgTS = [cs_avgTS, ssf_avgTS]; % Vector of average max stresses for bulk and custom specimens

ssu.EngStress = ssu.data(:,3)./(1000*CrossA); % Engineering stress for bulk specimen 3

%% Calculate Strain

cs1.EngStrain = cs1.data(:,2)./L; % Fracture length for custom specimen 1
cs2.EngStrain = cs2.data(:,2)./L; % Fracture length for custom specimen 2
cs3.EngStrain = cs3.data(:,2)./L; % Fracture length for custom specimen 3

cs_FL = [cs1.EngStrain(end), cs2.EngStrain(end), cs3.EngStrain(end)]; % Fracture length for custom specimens

cs_avgFL = mean(cs_FL); % Average fracture length of custom specimens

cs_FL_SD = std(cs_FL); % Standard deviation of fracture length of custom specimens

ssf1.EngStrain = ssf1.data(:,2)./L; % Fracutre length for bulk specimen 1
ssf2.EngStrain = ssf2.data(:,2)./L; % Fracutre length for bulk specimen 2

ssf_FL = [ssf1.EngStrain(end), ssf2.EngStrain(end)]; % Fracutre length for bulk specimens

ssf_avgFL = mean(ssf_FL); % Average fracture length of bulk specimens

ssf_FL_SD = std(ssf_FL); % Standard deviation of fracture length of bulk specimens

s_FLerr = [cs_FL_SD/2 ssf_FL_SD/2]; % Vector of one standard deviation for both specimens

avgFL = [cs_avgFL, ssf_avgFL]; % Vector of average fracture lengths for both specimens

ssu.EngStrain = ssu.data(:,2)./L;

%% Calculate Young's Modulus

cs1.YoungMod = (cs1.EngStress(10)-cs1.EngStress(5))/(cs1.EngStrain(10)-cs1.EngStrain(5)); % Young's modulus for custom specimen 1
cs2.YoungMod = (cs2.EngStress(10)-cs2.EngStress(5))/(cs2.EngStrain(10)-cs2.EngStrain(5)); % Young's modulus for custom specimen 2
cs3.YoungMod = (cs3.EngStress(10)-cs3.EngStress(5))/(cs3.EngStrain(10)-cs3.EngStrain(5)); % Young's modulus for custom specimen 3

ssf1.YoungMod = (ssf1.EngStress(10)-ssf1.EngStress(5))/(ssf1.EngStrain(10)-ssf1.EngStrain(5)); % Young's modulus for bulk specimen 1
ssf2.YoungMod = (ssf1.EngStress(10)-ssf1.EngStress(5))/(ssf1.EngStrain(10)-ssf1.EngStrain(5)); % Young's modulus for bulk specimen 2

YoungMod_Vec = [cs1.YoungMod, cs2.YoungMod, cs3.YoungMod, ssf1.YoungMod, ssf2.YoungMod]; % Vector of young's modulus for both specimens

%% Calculate Poisson's Ratio

for i = 1:7
    cs1_strainyy_vec(i) = mean(cs1.dice(:,12,1));
    cs1_strainxx_vec(i) = mean(cs1.dice(:,11,1));
end

cs1.poisson = mean(cs1_strainyy_vec)/(mean(cs1_strainxx_vec))/DR;

for i = 1:7
    cs2_strainyy_vec(i) = mean(cs2.dice(:,12,1));
    cs2_strainxx_vec(i) = mean(cs2.dice(:,11,1));
end

cs2.poisson = mean(cs2_strainyy_vec)/(mean(cs2_strainxx_vec))/DR;

for i = 1:7
    cs3_strainyy_vec(i) = mean(cs3.dice(:,12,1));
    cs3_strainxx_vec(i) = mean(cs3.dice(:,11,1));
end

cs3.poisson = mean(cs3_strainyy_vec)/(mean(cs3_strainxx_vec))/DR;

for i = 1:18
    ssf1_strainyy_vec(i) = mean(ssf1.dice(:,12,1));
    ssf1_strainxx_vec(i) = mean(ssf1.dice(:,11,1));
end

ssf1.poisson = mean(ssf1_strainyy_vec)/mean(ssf1_strainxx_vec);

for i = 1:20
    ssf2_strainyy_vec(i) = mean(ssf2.dice(:,12,1));
    ssf2_strainxx_vec(i) = mean(ssf2.dice(:,11,1));
end

ssf2.poisson = mean(ssf2_strainyy_vec)/mean(ssf2_strainxx_vec);

poisson_vec = [cs1.poisson, cs2.poisson, cs3.poisson, ssf1.poisson, ssf2.poisson];

%% Plot Stress-Strain Curves

if PlotOn == true

figure(4)
plot(cs1.EngStrain, cs1.EngStress, Linewidth = 1.25, Color = [.9 0 .1]);
hold on
plot(cs2.EngStrain, cs2.EngStress, Linewidth = 1.25, Color = [0 .1 .9]);
plot(cs3.EngStrain, cs3.EngStress, Linewidth = 1.25, Color = [0 .55 .45]);
legend('Custom Specimen 1', 'Custom Specimen 2', 'Custom Specimen 3',...
    location = 'northwest');
xlabel('Strain [unitless]');
ylabel('Stress [MPa]');

figure(5)
plot(ssf1.EngStrain, ssf1.EngStress, Linewidth = 1.25, Color = [.9 0 .1]);
hold on
plot(ssf2.EngStrain, ssf2.EngStress, Linewidth = 1.25, Color = [0 .1 .9]);
legend('Bulk Acrylic 1', 'Bulk Acrylic 2', location = 'northwest');
xlabel('Strain [unitless]');
ylabel('Stress [MPa]');


figure(6)
xname = categorical({'Custom','Bulk'});
xname = reordercats(xname,{'Custom','Bulk'});
bar(xname, avgTS)
hold on
er1 = errorbar(xname, avgTS, s_TSerr, s_TSerr);
er1.Color = [0 0 0];
er1.LineStyle = 'none';
er1.LineWidth = 1;
er1.CapSize = 30;
xlabel('Specimen');
ylabel('Tensile Strength [MPa]');

figure(7)
bar(xname, avgFL);
hold on
er2 = errorbar(xname, avgFL, s_FLerr, s_FLerr);
er2.Color = [0 0 0];
er2.LineStyle = 'none';
er2.LineWidth = 1;
er2.CapSize = 30;
xlabel('Specimen');
ylabel('Strain at Fracture [unitless]');

figure(8)
xnname = categorical({'Custom 1', 'Custom 2', 'Custom 3', 'Bulk 1', 'Bulk 2'});
xnname = reordercats(xnname, {'Custom 1', 'Custom 2', 'Custom 3', 'Bulk 1', 'Bulk 2'});
bar(xnname, YoungMod_Vec);
xlabel('Specimen');
ylabel('Young''s Modulus [MPa]');

figure(9)
bar(xnname, poisson_vec);
xlabel('Specimen');
ylabel('Poisson''s Ratio');

end




