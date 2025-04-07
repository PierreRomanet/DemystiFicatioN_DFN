%% 
close all
clear all
clc



% Name of the file
save_name = 'network_rate_strengthening';

% Path of all the simulation
direc = ['../problems/' save_name '/'];

% Load data
f1 = loadandprocessdata(direc,'V','P','theta'); % Load slip velocity, state, and pressure
input = load_input(direc); % Load input file



% Reorganise data according to each fault
id_beg = 1;
for fault_id = 1:length(input)

    % Set id_end
    id_end = id_beg + length(input(fault_id).a)-1;

    % Create structure
    f(fault_id).V = f1.V(:,id_beg:id_end);
    f(fault_id).theta = f1.theta(:,id_beg:id_end);
    f(fault_id).P = f1.P(:,id_beg:id_end);
    % f(fault_id).Sn = f1.Sn(:,id_beg:id_end);
    f(fault_id).node = input(fault_id).node;
    f(fault_id).time = f1.time;

    % Calculate element
    f(fault_id).element =  (f(fault_id).node(:,2:end)+f(fault_id).node(:,1:end-1))/2;



    id_beg = id_end + 2;

end





%% Plot maximum velocity on each fault
figure('Position',[1 1 1000 800])
for i=1:length(f)
semilogy((f(1).time-10*365.25*86400)/60, max(f(i).V,[],2),'linewidth',2)
hold on

end


xlim([0 86400/60])
ylabel('Maximum slip velocity (m/s)')
xlabel('Time (min)')

set(gca,'FontSize',20)
set(gcf,'Color','w')


%% Plot pressure evolution on fault 10
% Create curvilinear length along the fault
f(10).ds = zeros(1,size(f(10).element,2)); 
f(10).ds(1) = 0  
ds = sqrt((f(10).node(1,2:end)-f(10).node(1,1:end-1)).^2+(f(10).node(2,2:end)-f(10).node(2,1:end-1)).^2);
for i=2:length(f(10).element)
    f(10).ds(i) = f(10).ds(i-1) + ds(i-1)/2 + ds(i)/2;
end

% Find index of beginning of injection
[~,idx] = min(abs(f(10).time-86400*10*365.25));

% skip time index 
skip = 50;

% Figure
figure('Position',[1 1 1000 800])
plot(f(10).ds, f(10).P(idx:skip:end,:)/1e6,'LineWidth',2)

ylabel('Pressure (MPa')
xlabel('Position along the fault (m)')
set(gca,'FontSize',20)
set(gcf,'Color','w')


