%% Load data
clear all
clc
close all
addpath '/Users/pierre/Dropbox/PapersinPreparation/Romanet2025_tripleJunction/matlab'
addpath '/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_v07/wrapper'




% Name of the file
save_name = 'network_rate_strengthening';
version = 'v07';


% Path of all the simulation
direc1 = ['/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_' version '/problems/' save_name '/'];

% Base level
base_level = 1e-20;
max_level = 1e0;
skip = 10;
unit = 1.0; % meter
unit_P = 1e6;

% Load data
f1 = loadandprocessdata(direc1,'V','P');
input = load_input(direc1);

% Reorganise data according to each fault
id_beg = 1;
for fault_id = 1:length(input)

    % Set id_end
    id_end = id_beg + length(input(fault_id).a)-1;

    % Create structure
    f(fault_id).V = f1.V(:,id_beg:id_end);
    % f(fault_id).theta = f1.theta(:,id_beg:id_end);
    f(fault_id).P = f1.P(:,id_beg:id_end);
    % f(fault_id).Sn = f1.Sn(:,id_beg:id_end);
    f(fault_id).node = input(fault_id).node;
    f(1).time = f1.time;

    % Calculate element
    f(fault_id).element =  (f(fault_id).node(:,2:end)+f(fault_id).node(:,1:end-1))/2;
    f(fault_id).nb_element = size(f(fault_id).element,2);
    
    % Calculate fault.ds
    f(fault_id).ds = sqrt((f(fault_id).node(1,1:end-1)-f(fault_id).node(1,2:end)).^2 + ...
                     (f(fault_id).node(2,1:end-1)-f(fault_id).node(2,2:end)).^2);

    % 
    if input(fault_id).nb_source>=1
        f(fault_id).index_injection = input(fault_id).index_injection;
    end

    id_beg = id_end + 2;

end
nb_fault = length(f);

%% Find all nodes 

nb_crossing = 0;
for fault_id1 =1:nb_fault
    for fault_id2=fault_id1+1:nb_fault
        %
        % for each pair of node
        for node_id1=1:f(fault_id1).nb_element+1
            for node_id2=1:f(fault_id2).nb_element+1
                % Calculate the distance between the two nodes
                dist = sqrt((f(fault_id1).node(1,node_id1)-f(fault_id2).node(1,node_id2))^2 ...
                           +(f(fault_id1).node(2,node_id1)-f(fault_id2).node(2,node_id2))^2);
            
                % If distance < min(ds)*01
                if (dist<0.01*min(min(f(fault_id1).ds),min(f(fault_id2).ds))) 
                    'dist'
                    dist
                    %
                    % Update value of nb_crossing
                    nb_crossing = nb_crossing+1;
                    %
                end
            
            end 
        end 
    end 
end 
%
crossing = zeros(nb_crossing,4);
crossing_id =0;
for fault_id1 =1:nb_fault
    for fault_id2=fault_id1+1:nb_fault
        %
        % for each pair of node
        for node_id1=1:f(fault_id1).nb_element+1
            for node_id2=1:f(fault_id2).nb_element+1
                % Calculate the distance between the two nodes
                dist = sqrt((f(fault_id1).node(1,node_id1)-f(fault_id2).node(1,node_id2))^2 ... 
                           +(f(fault_id1).node(2,node_id1)-f(fault_id2).node(2,node_id2))^2);
            
                % If distance < min(ds)*01
                if (dist<0.01*min(min(f(fault_id1).ds),min(f(fault_id2).ds))) 
                       crossing_id = crossing_id + 1;
                       crossing(crossing_id,1) = fault_id1;
                       crossing(crossing_id,2) = fault_id2;
                       crossing(crossing_id,3) = node_id1;
                       crossing(crossing_id,4) = node_id2;
                end
            
            end 
        end 
    end 
end 

'nb_crossing'
nb_crossing

%% Add injection point 
crossing = [10 10 f(10).index_injection f(10).index_injection; crossing];
nb_crossing = nb_crossing+1;
%% Build all crossing 


for crossing_id1=1:nb_crossing
    % Add the faults to the crossing
    Crossing(crossing_id1).CrossToFault = [crossing(crossing_id1,1), crossing(crossing_id1,2)];
    % Initicializing 
    Crossing(crossing_id1).CrossToCross = [];

    for crossing_id2=1:nb_crossing
        if crossing_id1==crossing_id2
            continue
        end

        if crossing(crossing_id1,1)==crossing(crossing_id2,1)|  ...
           crossing(crossing_id1,1)==crossing(crossing_id2,2)|  ...
           crossing(crossing_id1,2)==crossing(crossing_id2,1)|  ...
           crossing(crossing_id1,2)==crossing(crossing_id2,2)

            Crossing(crossing_id1).CrossToCross = [Crossing(crossing_id1).CrossToCross, crossing_id2];
            Crossing(crossing_id1).CrossToFault = unique([Crossing(crossing_id1).CrossToFault, crossing(crossing_id2,1), crossing(crossing_id2,2)]);
        end
    end
end


%% Create map
map = ones(nb_crossing)*Inf;
for i =1:nb_crossing
    map(i,i) = 0;
    for j=Crossing(i).CrossToCross
        
        if crossing(i,1) == crossing(j,1)
            map(i,j) = abs(crossing(i,3)-crossing(j,3));

        elseif crossing(i,1) == crossing(j,2)
            map(i,j) = abs(crossing(i,3)-crossing(j,4));

        elseif crossing(i,2) == crossing(j,1)
            map(i,j) = abs(crossing(i,4)-crossing(j,3));

        elseif crossing(i,2) == crossing(j,2)
            map(i,j) = abs(crossing(i,4)-crossing(j,4));

        end
    end
end

distances = dijkstra(map,1);

%% Build fault to node
f(crossing(crossing_id,1)).crossing = [];
f(crossing(crossing_id,1)).crossing_node = [];
for crossing_id = 1:nb_crossing
    f(crossing(crossing_id,1)).crossing = [f(crossing(crossing_id,1)).crossing, crossing_id];
    f(crossing(crossing_id,1)).crossing_node = [f(crossing(crossing_id,1)).crossing_node, crossing(crossing_id,3)];

    f(crossing(crossing_id,2)).crossing = [f(crossing(crossing_id,2)).crossing, crossing_id];
    f(crossing(crossing_id,2)).crossing_node = [f(crossing(crossing_id,2)).crossing_node, crossing(crossing_id,4)];
end

%% For each fault compute distance to injection

for fault_id = 1:nb_fault
    temp = zeros(length(f(fault_id).crossing),f(fault_id).nb_element);

    for i=1:length(f(fault_id).crossing)
        lala  = [fliplr(1:f(fault_id).crossing_node(i)-1) 1:f(fault_id).nb_element+1-f(fault_id).crossing_node(i)];
     

        temp(i,:) = double(lala)-.5+distances(f(fault_id).crossing(i));
    end

    f(fault_id).dist_inj = 0.1*min(temp,[],1);
end
f(1).dist_inj






%%
% Build vector of distance from injection
DistInj = [];
for i=1:nb_fault
    DistInj = [DistInj f(i).dist_inj];
end
DistInj = sort(unique(DistInj));
Pdist_inj = zeros(length(f(1).time),length(DistInj));
Vdist_inj = zeros(length(f(1).time),length(DistInj));

for time_id=1:length(f(1).time)
    for fault_id=1:nb_fault
        for el_id = 1:f(fault_id).nb_element
            [x,x_id] = min(abs(DistInj-f(fault_id).dist_inj(el_id)));
            if x <0.01
                Pdist_inj(time_id,x_id) = max([f(fault_id).P(time_id,el_id),Pdist_inj(time_id,x_id)]);
                Vdist_inj(time_id,x_id) = max([f(fault_id).V(time_id,el_id),Vdist_inj(time_id,x_id)]);
            end

        end
    end
end


%% Distance from injection (timeslot)
close all
figure('Position',[1 1 2500 1200])
subplot(2,3,1)
% time id
time_array = [500 600 650 700 772 900];
for i = 1:6
    time_id = time_array(i)
    subplot(2,3,i)
    yyaxis left
    for fault_id = 1:nb_fault
        semilogy(f(fault_id).dist_inj,f(fault_id).V(time_id,:),'b.','MarkerSize',5)
        hold on
    end
    % semilogy(DistInj,Vdist_inj(time_id,:),'k-','LineWidth',2)


    yyaxis right
    for fault_id = 1:nb_fault
        plot(f(fault_id).dist_inj,f(fault_id).P(time_id,:)/1e6,'r.','MarkerSize',5)
        hold on
    end
    % plot(DistInj,Pdist_inj(time_id,:)/1e6,'k-','LineWidth',2)

    set(gcf,'color','w')
    set(gca,'Fontsize',20)
    
    xlim([0 max(DistInj)])
    xlabel('Distance from injection (m)')

    yyaxis left
    ylabel('Slip Velocity (m/s)')
    ylim([1e-15 1e-5])

    yyaxis right
    ylabel('Pressure (MPa)')
    ylim([0 3])


    % Plot the faults
    yyaxis right
    for fault_id = 1:nb_fault   
        temp = 2.1*ones(size(f(fault_id).dist_inj))+fault_id*0.08;
        plot(f(fault_id).dist_inj,temp,'k-','LineWidth',2)

         
        x_temp = [min(f(fault_id).dist_inj) max(f(fault_id).dist_inj)];
        if x_temp(2)<max(DistInj)/2
            text(x_temp(2)+1,temp(1),['Fault ' num2str(fault_id)])
        else
            text(x_temp(1)-3.2,temp(1),['Fault ' num2str(fault_id)])
        end
    end



    ax = gca;
    ax.YAxis(1).Color = 'b';
    ax.YAxis(2).Color = 'r';

    title(seconds2duration(f(1).time(time_id)),'FontSize',17)
end

export_fig('Distance_from_injection.pdf')


%% Distance from injection (timeslot, maximum)
close all
figure('Position',[1 1 2500 1200])
subplot(2,3,1)
% time id
time_array = [500 600 650 700 772 900];
for i = 1:6
    time_id = time_array(i)
    subplot(2,3,i)
    yyaxis left
 
    semilogy(DistInj,Vdist_inj(time_id,:),'b-','LineWidth',2)


    yyaxis right
  
    plot(DistInj,Pdist_inj(time_id,:)/1e6,'r-','LineWidth',2)
    hold on




    set(gcf,'color','w')
    set(gca,'Fontsize',20)

    xlim([0 max(DistInj)])
    xlabel('Distance from injection (m)')

    yyaxis left
    ylabel('Slip Velocity (m/s)')
    ylim([1e-15 1e-5])

    yyaxis right
    ylabel('Pressure (MPa)')
    ylim([0 3])


       % Plot the faults
    yyaxis right
    for fault_id = 1:nb_fault   
        temp = 2.1*ones(size(f(fault_id).dist_inj))+fault_id*0.08;
        plot(f(fault_id).dist_inj,temp,'k-','LineWidth',2)

         
        x_temp = [min(f(fault_id).dist_inj) max(f(fault_id).dist_inj)];
        if x_temp(2)<max(DistInj)/2
            text(x_temp(2)+1,temp(1),['Fault ' num2str(fault_id)])
        else
            text(x_temp(1)-3.2,temp(1),['Fault ' num2str(fault_id)])
        end
    end


    ax = gca;
    ax.YAxis(1).Color = 'b';
    ax.YAxis(2).Color = 'r';



    



    title(seconds2duration(f(1).time(time_id)),'FontSize',17)
end

export_fig('Distance_from_injection_simpler.pdf')

%% Distance from injection, single figure

close all
figure('Position',[1 1 1500 1200])

% Fontsize for fault number
font_size = 20;


% time id
time_array = [500 600 650 700 772 900];
for i = 1:6
    time_id = time_array(i)
    yyaxis left
    
    if i==length(time_array)
        pl1 = semilogy(DistInj,Vdist_inj(time_id,:),'b--','LineWidth',3);
    else
        pl1 = semilogy(DistInj,Vdist_inj(time_id,:),'b-','LineWidth',3);
    end
    pl1.Color = [pl1.Color 0.8*(i/length(time_array))+0.05];
    hold on

    yyaxis right
  

    if i==length(time_array)
        pl2 = plot(DistInj,Pdist_inj(time_id,:)/1e6,'r--','LineWidth',3);
    else
        pl2 = plot(DistInj,Pdist_inj(time_id,:)/1e6,'r-','LineWidth',3);
    end
    pl2.Color = [pl2.Color 0.8*(i/length(time_array))+0.05];
    hold on




    set(gcf,'color','w')
    set(gca,'Fontsize',20)

    xlim([0 max(DistInj)])
    xlabel('Distance from injection (m)')

    yyaxis left
    ylabel('Slip Velocity (m/s)')
    ylim([1e-15 1e-5])

    yyaxis right
    ylabel('Pressure (MPa)')
    ylim([0 3])


       % Plot the faults
    yyaxis right
    for fault_id = 1:nb_fault   
        temp = 2.1*ones(size(f(fault_id).dist_inj))+fault_id*0.08;
        plot(f(fault_id).dist_inj,temp,'k-','LineWidth',2)

         
        x_temp = [min(f(fault_id).dist_inj) max(f(fault_id).dist_inj)];
        if x_temp(2)<max(DistInj)/2
            text(x_temp(2)+1,temp(1),['Fault ' num2str(fault_id)],'FontSize',font_size)
        else
            text(x_temp(1)-2,temp(1),['Fault ' num2str(fault_id)],'FontSize',font_size)
        end
    end


    ax = gca;
    ax.YAxis(1).Color = 'b';
    ax.YAxis(2).Color = 'r';



    



    % title(seconds2duration(f(1).time(time_id)),'FontSize',17)
end

export_fig('Distance_from_injection_single.pdf')
export_fig('Distance_from_injection_single.png')

%% Distance from injection (timeslot, colorcoded by fault)
close all
figure('Position',[1 1 1400 1000])

figure_letter = ['a)' 'b)' 'c)' 'd'];

% subplot(2,2,1)
time_sp = [86400*365.25*10-5*60,86400*365.25*10+5*60,86400*365.25*10+30*60,86400*365.25*10+12*60*60];

% Find closest id to time
time_array= zeros(size(time_sp));
for i=1:length(time_sp)
[~,time_array(i)] = min(abs(f(1).time-time_sp(i)));
time_array(i) = time_array(i)+1;
end



color = colororder('gem12');
for i = 1:4
    figure('Position',[1 1 1000 600])

    time_id = time_array(i)
    % subplot(2,2,i)
    yyaxis left
    for fault_id = 1:nb_fault
        pl1(fault_id) = semilogy(f(fault_id).dist_inj,f(fault_id).V(time_id,:),'.','MarkerSize',8,'Color',color(fault_id,:));
        hold on
    end
    % semilogy(DistInj,Vdist_inj(time_id,:),'k-','LineWidth',2)


    % yyaxis right
    % for fault_id = 1:nb_fault
    %     pl2(fault_id) = plot(f(fault_id).dist_inj,f(fault_id).P(time_id,:)/1e6,'.','MarkerSize',5);
    %     hold on
    % end
    % plot(DistInj,Pdist_inj(time_id,:)/1e6,'k-','LineWidth',2)

    set(gcf,'color','w')
    set(gca,'Fontsize',25)
    
    xlim([0 max(DistInj)])
    xlabel('Distance from injection (m)')

    yyaxis left
    ylabel('Slip Velocity (m/s)')
    ylim([1e-15 1e-5])

    yyaxis right
    plot(DistInj,Pdist_inj(time_id,:)/1e6,'k-','LineWidth',3)
    ylabel('Pressure (MPa)','Rotation',-90)
    ylim([0 2])


    % Plot the faults
    yyaxis right
    for fault_id = 1:nb_fault   
        temp = 1.4*ones(size(f(fault_id).dist_inj))+fault_id*0.05;
        plot(f(fault_id).dist_inj,temp,'-','LineWidth',2,'Color',color(fault_id,:))

         
        x_temp = [min(f(fault_id).dist_inj) max(f(fault_id).dist_inj)];
       
        size_font = 20
        if fault_id==7 
            text(x_temp(2)+1,temp(1),['Fault ' num2str(fault_id)],'FontSize',size_font)
        elseif fault_id==9
            text(x_temp(1)-2.2,temp(1),['Fault ' num2str(fault_id)],'FontSize',size_font)
        elseif x_temp(2)<max(DistInj)/2
            text(x_temp(2)+1,temp(1),['Fault ' num2str(fault_id)],'FontSize',size_font)
        else 
            text(x_temp(1)-3.2,temp(1),['Fault ' num2str(fault_id)],'FontSize',size_font)
        end
    end



    ax = gca;
    ax.YAxis(1).Color = 'k';
    ax.YAxis(2).Color = 'k';
    % ax.YAxis(2).Visible = 'off'
    if f(1).time(time_id)-86400*365.25*10>0
        title([seconds2duration_special(f(1).time(time_id)-86400*365.25*10) ' after injection'],'FontSize',25)
    else
        title([seconds2duration_special(abs(f(1).time(time_id)-86400*365.25*10)) ' before injection'],'FontSize',25)
    end


    export_fig(['Distance_from_injection_fault_colorcoded' num2str(i) '.pdf'])

end

export_fig('Distance_from_injection_fault_colorcoded.pdf')



%% Calculate velocity of the slow slip
figure
% Pressure
k = 1e-15;
beta = 1e-8;
phi = 0.1;
T = 10000.0;
Pmax = 2e6;
eta = 1e-3;
Q0 = Pmax*sqrt(k*beta*phi)/sqrt(eta*T)
t_max = 10000;



t = 0:1:t_max;
dist = sqrt(k/(beta*eta*phi)*t);



[~,lala] = min(abs(Vdist_inj-1e-11),[],2);

plot(f(1).time-86400*365.25*10,DistInj(lala),'LineWidth',2)
xlim([0 t_max])
hold on

plot(t,dist,'LineWidth',2)

k = 2e-14
dist = sqrt(k/(beta*eta*phi)*t);
plot(t,dist,'LineWidth',2)
plot([2000 2000],[0 15],'k','LineWidth',2)


xlabel('Time since injection (s)')
ylabel('Distance (m) of Slip velocity front (10^{-11} m/s)')

 set(gcf,'color','w')
    set(gca,'Fontsize',18)

    export_fig('scaling_dist_time_diffusion.pdf')


%% Geometry
offset = 0.3;
figure('Position',[1 1 1000 1000])
for fault_id = 1:length(input)
    hold on
    plot(f(fault_id).node(1,:),f(fault_id).node(2,:),'k-','LineWidth',4,'Color',color(fault_id,:))
end
for fault_id = 1:length(input)

    % Plot fault number
    x_mean = (f(fault_id).node(1,1)+f(fault_id).node(1,end))/2;
    y_mean = (f(fault_id).node(2,1)+f(fault_id).node(2,end))/2;
    y_dist = max(f(fault_id).node(2,:))-min(f(fault_id).node(2,:));
    x_dist = max(f(fault_id).node(1,:))-min(f(fault_id).node(1,:));
    angle = 180/pi*atan(y_dist/x_dist);
    n1 = -y_dist/sqrt(x_dist^2+y_dist^2);
    n2 = x_dist/sqrt(x_dist^2+y_dist^2);
    text(x_mean+offset*n1,y_mean+offset*n2,['Fault ' num2str(fault_id)],'Rotation',angle,'FontSize',17,'HorizontalAlignment','center')
end

% Plot injection point
plot(f(10).element(1,input(10).index_injection),f(10).element(2,input(10).index_injection),'ok','MarkerFaceColor','k','MarkerSize',10)



axis equal
xlim([-6,11])
ylim([-3 17])
xlabel('Position (m)')
ylabel('Position (m)')


set(gca,'FontSize',18)
set(gcf,'color','w');
export_fig([save_name '_geometry.pdf'])


%% Calculate velocity of the slow slip
% Create fault.maxV
for fault_id=1:length(f)
    f(fault_id).maxV = max(f(fault_id).V,[],2);
    if max(f(fault_id).maxV)>1e-11

        index_temp = find(f(fault_id).maxV-1e-11>0);
        threshold_beg = index_temp(1);
        threshold_end = index_temp(end);
        f(fault_id).t_beg = f(1).time(threshold_beg)-86400*365.25*10;
        f(fault_id).t_end = f(1).time(threshold_end)-86400*365.25*10;

    end
end




figure





% [~,lala] = min(abs(Vdist_inj-1e-11),[],2);
index_dist = int16.empty;
for time_id =1:length(f(1).time)
    lala = find(Vdist_inj(time_id,:)-1e-11>0)
    if isempty(lala)
        index_dist = [index_dist 1];
    else
        index_dist = [index_dist lala(end)];
    end
    
end
index_dist = squeeze(index_dist);


Vdist_inj_threshold = zeros(size(Vdist_inj));
Vdist_inj_threshold(lala) = Vdist_inj(lala);


hold on


% Plot rectangle
for fault_id =1:length(f)
    if max(f(fault_id).maxV)>1e-11

        % Last time at witch there is 

        % Max and min of distance from injection 
        
        max_dist = max(f(fault_id).dist_inj);
        min_dist = min(f(fault_id).dist_inj);

        rectangle('position',[f(fault_id).t_beg min_dist f(fault_id).t_end-f(fault_id).t_beg max_dist-min_dist],'FaceColor',color(fault_id,:),'FaceAlpha',0.1,'EdgeColor','k','LineWidth',0.25)
    end
end


% Plot theoretical scaling
k = input(1).permeability(1);
beta = input(1).fluid_comp+input(1).rock_comp;
phi = input(1).porosity;
T = 10000.0;

eta = input(1).dyn_viscosity;

t_max = 10000;


t = 0:1:t_max;
dist = sqrt(k/(beta*eta*phi)*t);

plot(f(1).time-86400*365.25*10,DistInj(index_dist),'r','LineWidth',3)
xlim([0 t_max])
ylim([0 10])

plot(t,dist,'b','LineWidth',2)

k = 2e-14
dist = sqrt(k/(beta*eta*phi)*t);
plot(t,dist,'b','LineWidth',2)
plot([2000 2000],[0 15],'k--','LineWidth',2)
text(2250,2,'End of injection','FontSize',18,'Rotation',90)

text(2750,8.2,'\alpha=2x10^{-2} m^2/s','FontSize',18,'Rotation',42,'Color','b')
text(7500,3.2,'\alpha=10^{-3} m^2/s','FontSize',18,'Rotation',8,'Color','b')



% Plot text
for fault_id =1:length(f)
    if fault_id == 7 ||fault_id==8
        continue
    end

    if max(f(fault_id).maxV)>1e-11
        % Max and min of distance from injection 
        max_dist = max(f(fault_id).dist_inj);
        min_dist = min(f(fault_id).dist_inj);
        text(f(fault_id).t_beg+2400,min_dist+0.3,['Fault ' num2str(fault_id)],'FontSize',18)
    end
end

max_dist = max(f(7).dist_inj);
min_dist = min(f(7).dist_inj);
text(f(7).t_beg+450,min_dist+1.5,['Fault ' num2str(7)],'FontSize',18,'Rotation',90)

max_dist = max(f(8).dist_inj);
min_dist = min(f(8).dist_inj);
text(f(8).t_beg+200,min_dist+0.3,['Fault ' num2str(8)],'FontSize',18)



xlabel('Time since injection (s)')
ylabel({'Distance from injection (m)'; 'of slip velocity front (10^{-11} m/s)'})

 set(gcf,'color','w')
    set(gca,'Fontsize',18)

    export_fig('-m3','scaling_dist_time_diffusion_improved.png')



%% Find shortest distance from injection (Dijkstra algorithm)
function distances = dijkstra(map,startingpoint)

% Number of distances
N = size(map,1);

% Initialize distance
distances(1:N) = Inf;
distances(startingpoint) = 0;

% Make a visited array for each node
visited(startingpoint) = 1;
visited(1:N) = 0;

while isinf(sum(distances))

    candidates(1:N) = Inf;

    for index = 1:N
        if visited(index) == 0
            candidates(index) = distances(index);
        end
    end

    % Compute the distance
    [currentDistance, currentPoint] = min(candidates);

    for index=1:N
        newDistance = currentDistance + map(currentPoint, index);
        if newDistance<distances(index)
            distances(index) = newDistance ;
        end
    end

    % Mark the current 
    visited(currentPoint) = 1;

end


end


%%
nu = (input(1).cp^2-2*input(1).cs^2)/(0.5*(input(1).cp^2-input(1).cs^2))

%%

alpha = 1e-16/(input(1).porosity*input(1).dyn_viscosity*(input(1).fluid_comp+input(1).rock_comp))
alpha = 1e-15/(input(1).porosity*input(1).dyn_viscosity*(input(1).fluid_comp+input(1).rock_comp))

%%
Lb = -input(1).mu*input(1).Dc(1)*2*(1-input(1).cs^2/input(1).cp^2)/(input(1).sigmaN(1)*input(1).b(1))

%%
max(max(f(10).P))