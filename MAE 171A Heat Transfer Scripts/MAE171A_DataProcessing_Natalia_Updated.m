%% MAE 171A: Heat Transfer Lab
% Created 1/27/2023

clear; clc; close all
%% Import Data
    PlotOn = false;
    PlotOn2 = true;

% Parameters
    nu = 2.317e-5; % Kinemativ viscocity of air at 25 C [m^2/s]
    k = 3.186e-2; % Thermal conductivity of copper [W/(mK)]
    rho = 0.9413; % Density of air [kg/m^3]
    c_p = 1010.6; % heat capacity of air [J/(Kg*K)]
    Pr = 0.7; % Prandtl Number
    L = 0.009; % Distance from duct to plate [m]
    Sens_Res = 5; % Resistance of sensor
    A = pi*(0.0762/2)^2; % Area of plate [m^2]

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

% Plate Temperature
    PlateCal.Temps = [0, 47.8224, 70.956, 95.256, 124.42, 148.91]; % Plate temperature
    PlateCal.Pad_Voltage = [0, 174.586, 173.188, 170.991, 167.795, 166.396]; % Heating pad voltage
    PlateCal.Sens_Voltage = [0, 0.0733, 0.1613, 0.279, 0.4301, 0.489]; % Sensor voltage
    PlateCal.Sens_Curr = PlateCal.Sens_Voltage./Sens_Res; % Sensor current

%Plate Power
    PlateCal.Power = (PlateCal.Pad_Voltage - PlateCal.Sens_Voltage) .* PlateCal.Sens_Curr; % Plate power
    PlateCal.HeatFlux = PlateCal.Power./A; % Heat Flux

if PlotOn == true

        % Plot Anemometer 0 in
            figure(1)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.D0, AnemCal.Anem_Flow.D0, Color = [.9 0 .1]);
            title('Air Speed vs. Percent Fan Speed for Anemometer at 0 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.D0, AnemCal.Anem_Flow.D0, Color = [0 0 1]);
            title('Air Speed vs. Recorded Fan Speed for Anemometer at 0 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);
               
        % Plot Anemometer 0.4 in
            figure(2)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.Dpt4, AnemCal.Anem_Flow.Dpt4, Color = [.9 0 .1]);
            title('Air Speed vs. Percent Fan Speed for Anemometer at 0.4 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.Dpt4, AnemCal.Anem_Flow.Dpt4, Color = [0 0 1]);
            title('Air Speed vs. Recorded Fan Speed for Anemometer at 0.4 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);

        % Plot Anemometer 0.8 in
            figure(3)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.Dpt8, AnemCal.Anem_Flow.Dpt8, Color = [.9 0 .1]);
            title('Air Speed vs. Percent Fan Speed for Anemometer at 0.8 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.Dpt8, AnemCal.Anem_Flow.Dpt8, Color = [0 0 1]);
            title('Air Speed vs. Recorded Fan Speed for Anemometer at 0.8 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);
                
        % Plot Anemometer 1.2 in
            figure(4)
            subplot(2,1,1)
            plot(AnemCal.Recorded_PercentSpd.D1pt2, AnemCal.Anem_Flow.D1pt2, Color = [.9 0 .1]);
            title('Air Speed vs. Percent Fan Speed for Anemometer at 1.2 in. From Center');
            xlabel('Fan Speed [%]');
            ylabel('Air Speed [m/s]');
        
            subplot(2,1,2)
            plot(AnemCal.Recorded_Flow.D1pt2, AnemCal.Anem_Flow.D1pt2, Color = [0 0 1]);
            title('Air Speed vs. Recorded Fan Speed for Anemometer at 1.2 in. From Center');
            xlabel('Fan Speed [Unitless]');
            ylabel('Air Speed [m/s]');
            xlim([9 43]);

        % Plot Plate Power Calibration
            figure(5)
            plot(PlateCal.Temps, PlateCal.Power);
            title('Plate Power vs. Plate Temperature');
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

% Experiment 1
    
    Temp75.ActualPlateTemp = [72.9, 72.32, 72.51, 71.34, 71.34, 71.15];
    Temp75.ActualAirTemp = [21.97, 22.74, 22.36, 23.52, 23.13, 23.72];

    Temp75.h_coeff = interp1(PlateCal.Temps, PlateCal.HeatFlux, 75)...
        ./(Temp75.ActualPlateTemp - Temp75.ActualAirTemp); % Heat transfer coefficient for experiment 1
    Temp75.NuL = interp1(PlateCal.Temps, PlateCal.HeatFlux, 75).*L./(k.*...
        (Temp75.ActualPlateTemp - Temp75.ActualAirTemp)); % Nusselt number for experiment 1

% Experiment 2

    Temp100.ActualPlateTemp = [97.98, 97.01, 96.62, 95.45, 94.87, 94.67]; % Recorded plate temp at various flow speeds [C]
    Temp100.ActualAirTemp = [23.91, 23.71, 24.49, 22.94, 24.11, 24.49]; % Recorded air temp at various flow speeds [C]

    Temp100.h_coeff = interp1(PlateCal.Temps, PlateCal.HeatFlux, 100)...
        ./(Temp100.ActualPlateTemp - Temp100.ActualAirTemp); % Heat transfer coefficient for experiment 2
    Temp100.NuL = interp1(PlateCal.Temps, PlateCal.HeatFlux, 100).*L./(k.*...
        (Temp100.ActualPlateTemp - Temp100.ActualAirTemp)); % Nusselt number for experiment 2

% Experiment 3

    Temp125.ActualPlateTemp = [121.69, 120.92, 119.94, 118.78, 118.39,116.0]; % Recorded plate temp at various flow speeds [C]
    Temp125.ActualAirTemp = [24.11, 23.91, 23.13, 23.72, 24.88,24.52]; % Recorded air temp at various flow speeds [C]

    Temp125.h_coeff = interp1(PlateCal.Temps, PlateCal.HeatFlux, 125)...
        ./(Temp125.ActualPlateTemp - Temp125.ActualAirTemp); % Heat transfer coefficient for experiment 3
    Temp125.NuL = interp1(PlateCal.Temps, PlateCal.HeatFlux, 125).*L./(k.*...
        (Temp125.ActualPlateTemp - Temp125.ActualAirTemp)); % Nusselt number for experiment 3

% Experiment 4

    Temp150.ActualPlateTemp = [145.22, 144.05, 143.66, 141.52, 141.33, 140.16]; % Recorded plate temp at various flow speeds [C]
    Temp150.ActualAirTemp = [23.33, 22.74, 23.13, 22.55, 22.94, 22.74]; % Recorded air temp at various flow speeds [C]

    Temp150.h_coeff = PlateCal.HeatFlux(6)...
        ./(Temp150.ActualPlateTemp - Temp150.ActualAirTemp); % Heat transfer coefficient for experiment 4
    Temp150.NuL = PlateCal.HeatFlux(6).*L./(k.*...
        (Temp150.ActualPlateTemp - Temp150.ActualAirTemp)); % Nusselt number for experiment 4

if PlotOn2 == true

% Experiment 1 Plots
    figure(6)
    subplot(2,1,1)
    plot(sqrt(Exp.Re_vec), Temp75.h_coeff, Linewidth = 1, Color = [.9 0 .1]);
%     title('Heat Transfer Coefficient [H] vs. Square Root of Reynold''s Number at 75 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('h');
    subplot(2,1,2)
    plot(sqrt(Exp.Re_vec), Temp75.NuL, Linewidth = 1, Color = [.1 0 .9]);
%     title('Experimental Nusselt Number [Nu_L] vs. Square Root of Reynold''s Number at 75 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('Nu_L');

% Experiment 2 Plots
    figure(7)
    subplot(2,1,1)
    plot(sqrt(Exp.Re_vec), Temp100.h_coeff, Linewidth = 1, Color = [.9 0 .1]);
%     title('Heat Transfer Coefficient [H] vs. Square Root of Reynold''s Number at 100 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('h');
    subplot(2,1,2)
    plot(sqrt(Exp.Re_vec), Temp100.NuL, Linewidth = 1, Color = [.1 0 .9]);
%     title('Experimental Nusselt Number [Nu_L] vs. Square Root of Reynold''s Number at 100 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('Nu_L');

% Experiment 3 Plots
    figure(8)
    subplot(3,1,1)
%     err = [5 8 2 9 3 3];
    plot(sqrt(Exp.Re_vec), Temp125.h_coeff, Linewidth = 1, Color = [.9 0 .1], Marker='diamond');
    x=sqrt(Exp.Re_vec);
    y=Temp125.h_coeff;
%     err=.05*y;
%     neg=err./2
%     pos=err./2
    errba=y*.02
    errba=mean(errba)
    err = errba.*ones(size(y));
    neg=err./2
    pos=err./2
     errorbar(x,y,neg,pos,'-ro')
%     title('Heat Transfer Coefficient [H] vs. Square Root of Reynold''s Number at 125 Degrees C');
%     xlabel('Square Root of Reynold''s Number');
    ylabel('h');

    subplot(3,1,2)
%     plot(sqrt(Exp.Re_vec), Temp125.NuL, Linewidth = 1, Color = [.1 0 .9], Marker='diamond');
    x=sqrt(Exp.Re_vec);
    y=Temp125.NuL;
%     err=.05*y;
%     neg=err./2
%     pos=err./2
    errba2=y*.03
    errba2=mean(errba2)
    err = errba2.*ones(size(y));
    neg=err./2
    pos=err./2
% pos=NaN
     errorbar(x,y,neg,pos,'-bo')
% errorbar(x,y,err,'-bo')
    ylabel(' Experimental Nu_L');
% xlabel('Square Root of Reynold''s Number');

    subplot(3,1,3)
    exp_val=[19.7467 23.1481 26.2126 31.5600 33.6351 36.1095]
    num=length((sqrt(Exp.Re_vec)))
    x=linspace(0,55,num)
    plot(sqrt(Exp.Re_vec),exp_val,Linewidth = 1,Color='black', Marker='diamond')
%     title('Experimental Nusselt Number [Nu_L] vs. Square Root of Reynold''s Number at 125 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('Theoretical Nu_L');




% Experiment 4 Plots
    figure(9)
    subplot(2,1,1)
    plot(sqrt(Exp.Re_vec), Temp150.h_coeff, Linewidth = 1, Color = [.9 0 .1]);
%     title('Heat Transfer Coefficient [H] vs. Square Root of Reynold''s Number at 150 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('h');
    subplot(2,1,2)
    plot(sqrt(Exp.Re_vec), Temp150.NuL, Linewidth = 1, Color = [.1 0 .9]);
%     title('Experimental Nusselt Number [Nu_L] vs. Square Root of Reynold''s Number at 150 Degrees C');
    xlabel('Square Root of Reynold''s Number');
    ylabel('Nu_L');

end
%%
PlotOn==false
if PlotOn == true
    figure(10)
    data=readtable('temp125Cwind45.csv','NumHeaderLines',1);
    time=data.Var1;
    temp=data.Var2;
    amb_temp=data.Var6;
    % plot(data.Var1,data.Var2)
    secondstime= datetime(time)
    time_vec=1:length(secondstime);
    % start=secondstime(1)
    % endtime=secondstime(end)
    % counting=secondstime-start
    p=plot(time_vec,temp,'k')
    % hold on
    % plot(time_vec,amb_temp)
    p.LineWidth=1
    xlabel('Time (seconds)')
    ylabel(['Temperature (' char(176) 'C)'])
end

