%% MAE 171A: Heat Transfer Lab
% Created 1/27/2023

clear; clc; close all
%% Plot Commands
    PlotOn = false;
    PlotOn2 = true;

%% Parameters
    nu = 2.317e-5; % Kinematic viscocity of air at 25 C [m^2/s]
    k = 3.186e-2; % Thermal conductivity of copper [W/(mK)]
    rho = 0.9413; % Density of air [kg/m^3]
    c_p = 1010.6; % heat capacity of air [J/(Kg*K)]
    Pr = 0.7; % Prandtl Number
    L = 0.009; % Distance from duct to plate [m]
    Sens_Res = 5; % Resistance of sensor
    A = pi*(0.0762/2)^2; % Area of plate [m^2]
    temp = 75:25:150;
    temp2 = 50:25:150;
    spd = [35, 45 55, 75, 85, 100];

%% Anemometer Calibration

% Anemometer 0 in
    AnemCal.Recorded_Flow.D0 = [9, 18, 27, 36, 43]; % Data reported air speed [unitless]
    AnemCal.Recorded_PercentSpd.D0 = [20, 40, 60, 80, 100]; % Recorded fan speed [%]
    AnemCal.Anem_Flow.D0 = [1.1, 2.8, 4.6, 6.3, 7.8]; % Anemometer air velocity [m/s]

% Anemometer 0.4 in
    AnemCal.Recorded_Flow.Dpt4 = [9, 18, 27, 36, 43]; % Data reported air speed [unitless]
    AnemCal.Recorded_PercentSpd.Dpt4 = [20, 40, 60, 80, 100]; % Recorded fan speed [%]
    AnemCal.Anem_Flow.Dpt4 = [1.1, 2.8, 4.9, 6.5, 7.9]; % Anemometer air velocity [m/s]

% Anemometer 0.8 in
    AnemCal.Recorded_Flow.Dpt8 = [9, 18, 27, 36, 43]; % Data reported air speed [unitless]
    AnemCal.Recorded_PercentSpd.Dpt8 = [20, 40, 60, 80, 100]; % Recorded fan speed [%]
    AnemCal.Anem_Flow.Dpt8 = [1.1, 2.7, 4.3, 6.3, 7.6]; % Anemometer air velocity [m^3/s]

% Anemometer 1.2 in
    AnemCal.Recorded_Flow.D1pt2 = [9, 18, 27, 36, 43]; % Data reported air speed [unitless]
    AnemCal.Recorded_PercentSpd.D1pt2 = [20, 40, 60, 80, 100]; % Recorded fan speed [%]
    AnemCal.Anem_Flow.D1pt2 = [1, 2.5, 4.1, 6.1, 7.3]; % Anemometer air velocity [m^3/s]
%% Plate Calibration

% Plate Temperature

    for i = 1:length(temp2)
        filename = sprintf('Heat Transfer %iC', temp2(i));
        matrix = readmatrix(filename);
        PlateCal.Temps(i) = mean(matrix(end-30:end, 2));
        PlateCal.Pad_Voltage(i) = mean(matrix(end-30:end, 3));
        PlateCal.Sens_Voltage(i) = mean(matrix(end-30:end, 4));
        PlateCal.Sens_Curr = PlateCal.Sens_Voltage./Sens_Res; %Sensor Current
    end

% Plate Power
    PlateCal.Power = (PlateCal.Pad_Voltage - PlateCal.Sens_Voltage) .* PlateCal.Sens_Curr; % Plate power
    PlateCal.HeatFlux = PlateCal.Power./A; % Heat Flux

if PlotOn == true

        % Plot Anemometer 0 in
            figure(1)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0,...
                Linewidth = 1.25, Marker = "diamond" , Color = [.85 0 .14]);
            grid on
            % title('Air Speed vs. Percent Fan Speed for Anemometer at 0 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.D0, AnemCal.Anem_Flow.D0,...
                Linewidth = 1.25, Marker = "diamond" , Color = [0 .1 .9]);
            grid on
            % title('Air Speed vs. Recorded Fan Speed for Anemometer at 0 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);
               
        % Plot Anemometer 0.4 in
            figure(2)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4,...
                Linewidth = 1.25, Marker = "diamond" , Color = [.85 0 .14]);
            grid on
            % title('Air Speed vs. Percent Fan Speed for Anemometer at 0.4 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.Dpt4, AnemCal.Anem_Flow.Dpt4,...
                Linewidth = 1.25, Marker = "diamond" , Color = [0 .1 .9]);
            grid on
            % title('Air Speed vs. Recorded Fan Speed for Anemometer at 0.4 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);

        % Plot Anemometer 0.8 in
            figure(3)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8,...
                Linewidth = 1.25, Marker = "diamond" , Color = [.85 0 .14]);
            grid on
            % title('Air Speed vs. Percent Fan Speed for Anemometer at 0.8 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.Dpt8, AnemCal.Anem_Flow.Dpt8,...
                Linewidth = 1.25, Marker = "diamond" , Color = [0 .1 .9]);
            grid on
            % title('Air Speed vs. Recorded Fan Speed for Anemometer at 0.8 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);
                
        % Plot Anemometer 1.2 in
            figure(4)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2,...
                Linewidth = 1.25, Marker = "diamond" , Color = [.85 0 .14]);
            grid on
            % title('Air Speed vs. Percent Fan Speed for Anemometer at 1.2 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.D1pt2, AnemCal.Anem_Flow.D1pt2,...
                Linewidth = 1.25, Marker = "diamond" , Color = [0 .1 .9]);
            grid on
            % title('Air Speed vs. Recorded Fan Speed for Anemometer at 1.2 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);

        % Plot Plate Power Calibration
            figure(5)
            plot(PlateCal.Temps, PlateCal.Power,...
                Linewidth = 1.25, Marker = "diamond" , Color = [.85 0 .14]);
            grid on
            % title('Plate Power vs. Plate Temperature');
            xlabel('Plate Temperature [C]');
            ylabel('Plate Power [W]');
end

% 35% Flow Speed

    Exp.Flow35.spd0 = interp1(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, 35);
    Exp.Flow35.spdpt4 = interp1(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, 35);
    Exp.Flow35.spdpt8 = interp1(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, 35);
    Exp.Flow35.spd1pt2 = interp1(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, 35);
    
    Exp.Flow35.AvgSpd = (Exp.Flow35.spd0 + Exp.Flow35.spdpt4 + Exp.Flow35.spdpt8 + ...
        Exp.Flow35.spd1pt2) / 4; 
    
    Exp.Flow35.Re = Exp.Flow35.AvgSpd * .009 / nu; % Reynold's Number for Exp at 35% fan speed
    Exp.Flow35.NuL_Theoretical = 0.763*sqrt(Exp.Flow35.Re)*Pr^(0.4); % Theoretical Nusselt Number

% 45% Flow Speed
    
    Exp.Flow45.spd0 = interp1(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, 45);
    Exp.Flow45.spdpt4 = interp1(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, 45);
    Exp.Flow45.spdpt8 = interp1(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, 45);
    Exp.Flow45.spd1pt2 = interp1(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, 45);
    
    Exp.Flow45.AvgSpd = (Exp.Flow45.spd0 + Exp.Flow45.spdpt4 + Exp.Flow45.spdpt8 + ...
       Exp.Flow45.spd1pt2) / 4; 
    
    Exp.Flow45.Re = Exp.Flow45.AvgSpd * .009 / nu; % Reynold's Number for Exp at 45% fan speed
    Exp.Flow45.NuL_Theoretical = 0.764*sqrt(Exp.Flow45.Re)*Pr^(0.4); % Theoretical Nusslet Number

% 55% Flow Speed

    Exp.Flow55.spd0 = interp1(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, 55);
    Exp.Flow55.spdpt4 = interp1(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, 55);
    Exp.Flow55.spdpt8 = interp1(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, 55);
    Exp.Flow55.spd1pt2 = interp1(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, 55);
    
    Exp.Flow55.AvgSpd = (Exp.Flow55.spd0 + Exp.Flow55.spdpt4 + Exp.Flow55.spdpt8 + ...
        Exp.Flow55.spd1pt2) / 4; 
    
    Exp.Flow55.Re = Exp.Flow55.AvgSpd * .009 / nu; % Reynold's Number for Exp at 55% fan speed
    Exp.Flow55.NuL_Theoretical = 0.764*sqrt(Exp.Flow55.Re)*Pr^(0.4); % Theoretical Nusslet Number

% 75% Flow Speed

    Exp.Flow75.spd0 = interp1(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, 75);
    Exp.Flow75.spdpt4 = interp1(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, 75);
    Exp.Flow75.spdpt8 = interp1(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, 75);
    Exp.Flow75.spd1pt2 = interp1(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, 75);
    
    Exp.Flow75.AvgSpd = (Exp.Flow75.spd0 + Exp.Flow75.spdpt4 + Exp.Flow75.spdpt8 + ...
        Exp.Flow75.spd1pt2) / 4; 
    
    Exp.Flow75.Re = Exp.Flow75.AvgSpd * .009 / nu; % Reynold's Number for Exp at 75% fan speed
    Exp.Flow75.NuL_Theoretical = 0.764*sqrt(Exp.Flow75.Re)*Pr^(0.4); % Theoretical Nusslet Number

% 85% Flow Speed

    Exp.Flow85.spd0 = interp1(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, 85);
    Exp.Flow85.spdpt4 = interp1(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, 85);
    Exp.Flow85.spdpt8 = interp1(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, 85);
    Exp.Flow85.spd1pt2 = interp1(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, 85);
    
    Exp.Flow85.AvgSpd = (Exp.Flow85.spd0 + Exp.Flow85.spdpt4 + Exp.Flow85.spdpt8 + ...
        Exp.Flow85.spd1pt2) / 4; 
    
    Exp.Flow85.Re = Exp.Flow85.AvgSpd * .009 / nu; % Reynold's Number for Exp at 85% fan speed
    Exp.Flow85.NuL_Theoretical = 0.764*sqrt(Exp.Flow85.Re)*Pr^(0.4); % Theoretical Nusslet Number

% 100% Flow Speed

    Exp.Flow100.spd0 = interp1(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, 100);
    Exp.Flow100.spdpt4 = interp1(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, 100);
    Exp.Flow100.spdpt8 = interp1(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, 100);
    Exp.Flow100.spd1pt2 = interp1(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, 100);
    
    Exp.Flow100.AvgSpd = (Exp.Flow100.spd0 + Exp.Flow100.spdpt4 + Exp.Flow100.spdpt8 + ...
        Exp.Flow100.spd1pt2) / 4; 
    
    Exp.Flow100.Re = Exp.Flow100.AvgSpd * .009 / nu; % Reynold's Number for Exp at 100% fan speed
    Exp.Flow100.NuL_Theoretical = 0.764*sqrt(Exp.Flow100.Re)*Pr^(0.4); % Theoretical Nusslet Number

% Vector of Reynold's Number

Exp.Re_vec = [Exp.Flow35.Re, Exp.Flow45.Re, Exp.Flow55.Re, Exp.Flow75.Re,...
    Exp.Flow85.Re, Exp.Flow100.Re]; % Vector of Reynold's Numbers

%% Import Data

Temp.ActualPlateTemp = zeros(4, 6);
Temp.ActualAirTemp = zeros(4,6);

for i = 1:length(temp)
    for j = 1:length(spd)
       filename = sprintf('temp%iCwind%i', temp(i), spd(j));
       matrix = readmatrix(filename);
       Temp.ActualPlateTemp(i, j) = mean(matrix(end-30:end, 2));
       Temp.ActualAirTemp(i,j) = mean(matrix(end-30:end, 6));
    end
end


%% Calculating Nu and H

Temp.h_coeff = zeros(4,6);
Temp.NuL = zeros(4,6);

for i = 1:length(temp)-1
        Temp.h_coeff(i,:) = interp1(PlateCal.Temps, PlateCal.HeatFlux, temp(i))...
            ./(Temp.ActualPlateTemp(i,:) - Temp.ActualAirTemp(i,:)); % Heat transfer coeff
        Temp.NuL(i,:) = interp1(PlateCal.Temps, PlateCal.HeatFlux, temp(i)).*L./(k.*...
            (Temp.ActualPlateTemp(i,:) - Temp.ActualAirTemp(i,:))); % Nusselt number
end

% Since the plate never reaches 150C, the for loop cannot interpolate the
% heat flux, so a calculation is needed just for the 150C trial

Temp.h_coeff(end,:) = PlateCal.HeatFlux(end)./...
    (Temp.ActualPlateTemp(end,:) - Temp.ActualAirTemp(end,:)); % Heat transfer coeff for exp 4
Temp.NuL(end,:) = PlateCal.HeatFlux(end).*L./(k.*...
        (Temp.ActualPlateTemp(end,:) - Temp.ActualAirTemp(end,:))); % Nusselt number for exp 4

%% Plotting Experimental Data

if PlotOn2 == true

    for i = 1:length(temp)
    
        figure(i)
        subplot(2,1,1)
        plot(sqrt(Exp.Re_vec), Temp.h_coeff(i,:), Linewidth = 1.25, Color = [.9 0 .1], Marker = "diamond");
        xlabel('Square Root of Reynold''s Number');
        ylabel('Heat Transfer Coefficient [h]');
        axis([28 56 min(Temp.h_coeff(i,:))-.25 max(Temp.h_coeff(i,:))+.25]);
    
        subplot(2,1,2)
        plot(sqrt(Exp.Re_vec), Temp.NuL(i,:), Linewidth = 1.25, Color = [0 .1 .9], Marker = "diamond");
        xlabel('Square Root of Reynold''s Number');
        ylabel('Heat Transfer Coefficient [h]');
        axis([28 56 min(Temp.NuL(i,:))-.25 max(Temp.NuL(i,:))+.25]);

    end
end
