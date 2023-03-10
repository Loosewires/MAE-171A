% Plotting Open Loop

%Combined Mass

for i = 1:5
    filename = sprintf('MAE171A_OpenLoopStep_CombinedMass%i.txt', i);
    mat = readecp(filename);
    t = mat(:, 2);
    ed = mat(:, 3);
    
    figure(i)
    plot(t, ed);
end

% Mass 1

for i = 1:5
    filename = sprintf('MAE171A_OpenLoopStep_1Mass%i.txt', i);
    mat = readecp(filename);
    t = mat(:, 2);
    ed = mat(:, 3);
    
    figure(i+5)
    plot(t, ed);
end

% Mass 2

for i = 1:5
    filename = sprintf('MAE171A_OpenLoopStep_2Mass%i.txt', i);
    mat = readecp(filename);
    t = mat(:, 2);
    ed = mat(:, 4);
    
    figure(i+10)
    plot(t, ed)
end

