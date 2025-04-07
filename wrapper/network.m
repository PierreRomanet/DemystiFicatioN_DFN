%% Load data
clear all
clc
close all
addpath '/Users/pierre/Dropbox/PapersinPreparation/Romanet2023/matlab'
addpath '/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_v08/wrapper'




% Name of the file
save_name = 'network_rate_strengthening';
save_name = 'network_rate_strengthening_variable_Lk0.0005';

version = 'v08';

% save_name = 'network_rate_strengthening_Q8e-7_10000s';
% version = 'v05'

% Path of all the simulation
direc1 = ['/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_' version '/problems/' save_name '/'];



% Base level
base_level = 1e-20;
max_level = 1e0;
skip = 10;
unit = 1.0; % meter
unit_P = 1e6;





% Load data
f1 = loadandprocessdata(direc1,'V','P','theta');
input = load_input(direc1);



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

time_sp = [86400*365.25*10-5*60,86400*365.25*10+5*60,86400*365.25*10+30*60,86400*365.25*10+12*60*60];

% Find closest id to time
time_index= zeros(size(time_sp));
for i=1:length(time_sp)
[~,time_index(i)] = min(abs(f(1).time-time_sp(i)));
time_index(i) = time_index(i)+1;
end


%%

%% Make video


% Plot data 
% Close all figures
close all
figure('position', [0 0 1400 900])


% create the video writer with 1 fps
writerObj = VideoWriter(['myVideo_' save_name '_' version] ,'MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
% open the video writer
open(writerObj)


for i=2000:skip:length(f1.time)
    clf



hold on
   
    
    
    
    xlim([-6 10])
    ylim([-3 18])


     % Improve z axis
      % clb.Label.String = 'Slip velocity (m/s)';

    set(gca,'zscale','log')
    set(gca,'ZMinorTick','off')
    zlim([base_level  max_level])

    % Make surface 
    xl = xlim;
    yl = ylim;
    zl = zlim;
    [X,Y] = meshgrid([xl(1) xl(end)],[yl(1) yl(end)]);
    Z = zl(1)*ones(size(X));
    surf(X,Y,Z,'FaceColor',[246/256,215/256,176/256],'FaceAlpha',0.3,'EdgeColor','k','linewidth',2)
    


    for fault_id=1:length(input)
    % Make plot fault
     % Make plot velocity
    h = plot3(f(fault_id).element(1,:)/unit,f(fault_id).element(2,:)/unit,squeeze(f(fault_id).V(i,:)),'linewidth',3);
    % Link fault to velocity
    plot3([f(fault_id).element(1,1)/unit f(fault_id).element(1,1)/unit],[f(fault_id).element(2,1)/unit f(fault_id).element(2,1)/unit],[base_level*2  f(fault_id).V(i,1)],'k','linewidth',2)
    plot3([f(fault_id).element(1,end)/unit f(fault_id).element(1,end)/unit],[f(fault_id).element(2,end)/unit f(fault_id).element(2,end)/unit],[base_level*2  f(fault_id).V(i,end)],'k','linewidth',2)




     surface([f(fault_id).element(1,:)/unit;f(fault_id).element(1,:)/unit],...
        [f(fault_id).element(2,:)/unit;f(fault_id).element(2,:)/unit], ...
        [ones(size(f(fault_id).element(2,:)))*base_level*2;ones(size(f(fault_id).element(2,:)))*base_level*2],...
        [f(fault_id).P(i,:)/unit_P;f(fault_id).P(i,:)/unit_P],...
        'facecol','no',...
        'edgecol','interp',...
        'linewidth',12);


    end
    

    % Plot time
    text(-5,0,1000,seconds2duration(f1.time(i)),'FontSize',20)

    maxVel = max(f1.V(i,:));
    % Plot maximum velocity
    text(-5,0,10,{'Maximum slip rate :', [num2str(maxVel,'%2.5e') ' m/s']},'FontSize',20)


    xlabel('X position (km)')
    ylabel('Y position (km)')
    zlabel('Slip rate (m/s)')

    cb = colorbar;
    colormap(jet)
    cb.Label.String = 'Pressure (MPa)';
    cb.Label.FontSize = 20;
    cb.LineWidth = 2;
    cb.Location = 'eastoutside';
    a =  cb.Position;
    set(cb,'Position',[a(1)+0.06 a(2) a(3) a(4)*0.6])
   
    % cb.Position = [.95 0.1 0.02 0.5];
    % clim([0 1])

    set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

    view(20,35)
  



    writeVideo(writerObj, getframe(gcf));

end
   % export_fig([save_name '_fig_poster5.png'])

close(writerObj);



%% Make figure paper


% Plot data 
% Close all figures
close all







figure_letter = ['a.' 'b.' 'c.' 'd.'];
% tiledlayout(2,2)

for i = 1:length(time_index)
    figure('position', [0 0 1000 600])





    hold on




    xlim([-6 10])
    ylim([-3 18])


    % Improve z axis
    % clb.Label.String = 'Slip velocity (m/s)';

    set(gca,'zscale','log')
    set(gca,'ZMinorTick','off')
    zlim([base_level  max_level])

    % Make surface
    xl = xlim;
    yl = ylim;
    zl = zlim;
    [X,Y] = meshgrid([xl(1) xl(end)],[yl(1) yl(end)]);
    Z = zl(1)*ones(size(X));
    surf(X,Y,Z,'FaceColor',[246/256,215/256,176/256],'FaceAlpha',0.3,'EdgeColor','k','linewidth',2)



    for fault_id=1:length(input)
        % Make plot fault
        % Make plot velocity
        h = plot3(f(fault_id).element(1,:)/unit,f(fault_id).element(2,:)/unit,squeeze(f(fault_id).V(time_index(i),:)),'linewidth',3);
        % Link fault to velocity
        plot3([f(fault_id).element(1,1)/unit f(fault_id).element(1,1)/unit],[f(fault_id).element(2,1)/unit f(fault_id).element(2,1)/unit],[base_level*2  f(fault_id).V(time_index(i),1)],'k','linewidth',2)
        plot3([f(fault_id).element(1,end)/unit f(fault_id).element(1,end)/unit],[f(fault_id).element(2,end)/unit f(fault_id).element(2,end)/unit],[base_level*2  f(fault_id).V(time_index(i),end)],'k','linewidth',2)




        surface([f(fault_id).element(1,:)/unit;f(fault_id).element(1,:)/unit],...
            [f(fault_id).element(2,:)/unit;f(fault_id).element(2,:)/unit], ...
            [ones(size(f(fault_id).element(2,:)))*base_level*2;ones(size(f(fault_id).element(2,:)))*base_level*2],...
            [f(fault_id).P(time_index(i),:)/unit_P;f(fault_id).P(time_index(i),:)/unit_P],...
            'facecol','no',...
            'edgecol','interp',...
            'linewidth',8);


    end


    % Plot time
    % text(-5,0,1000,seconds2duration_special(f1.time(time_index(i))),'FontSize',20)
% Plot time 
    if f(1).time(time_index(i))-86400*365.25*10>0
        text(-5,0,1000,[seconds2duration_special(f(1).time(time_index(i))-86400*365.25*10) ' after injection'],'FontSize',25)
    else
        text(-5,0,1000,[seconds2duration_special(abs(f(1).time(time_index(i))-86400*365.25*10)) ' before injection'],'FontSize',25)
    end



    maxVel = max(f1.V(time_index(i),:));
    % Plot maximum velocity
    text(-5,0,1e-2,{'Maximum slip rate :', [num2str(maxVel,'%2.5e') ' m/s']},'FontSize',25)



    xlabel('X position (m)')
    ylabel('Y position (m)')
    zlabel('Slip rate (m/s)')
    
    cb = colorbar;
    colormap(jet)
    cb.Label.String = 'Pressure (MPa)';
    if i==1
        caxis([0 0.1])
    end
    cb.Label.FontSize = 25;
    cb.LineWidth = 2;
    cb.Location = 'eastoutside';
    a =  cb.Position;
    set(cb,'Position',[a(1) a(2) a(3) a(4)])

    % cb.Position = [.95 0.1 0.02 0.5];
    % clim([0 1])
        % text(-7.5,0,100000,figure_letter(i),'FontSize',30)

    set(gca,'FontSize',25)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

    view(20,35)

    
export_fig('-m3',[save_name num2str(i) '_fig_paper.pdf'])

cb.Visible = 'off'
export_fig('-m3',[save_name num2str(i) '_fig_paper.png'])



end


export_fig('-m3',[save_name '_fig_paper.pdf'])
export_fig('-m3',[save_name '_fig_paper.png'])




%% PLot max V evolution
close all
figure(2)

for fault_id = 1:length(input)
    semilogy(f1.time/(86400*365.25),max(f(fault_id).V,[],2),'k','linewidth',1)
 hold on

end

% semilogy(f1.time/(86400*365.25),max(f(1).V,[],2),'r','linewidth',2)
% xlim([0 100])


set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

ylabel('Max slip rate (m/s)')
xlabel('Time (year)')

   export_fig([save_name '_maxV.pdf'])

% figure(4)
% 
% 
% for fault_id = 1:length(input)
%     semilogy(f1.time/(86400*365.25),max(f(fault_id).V,[],2),'k','linewidth',1)
%  hold on
% 
% end
% semilogy(f1.time/(86400*365.25),max(f(1).V,[],2),'r','linewidth',2)
% 
% 
% set(gca,'FontSize',18)
%     set(gca,'linewidth',2)
%     set(gcf,'color','w');
% xlim([10 f1.time(end)/(86400*365.25)])
% ylabel('Max slip rate (m/s)')
% xlabel('Time (year)')
% 
% 
%    export_fig([save_name '_maxV2.pdf'])
% %%
% close all
% 
% figure(4)
% for fault_id = 1:length(input)
%     semilogy(f1.time(end:end-200),max(f(fault_id).V(end:end-200,:),[],2),'k','linewidth',1)
%  hold on
% 
% end
% semilogy(f1.time(end:end-200),max(f(1).V(end:end-200,:),[],2),'r','linewidth',2)
% 
% 
% set(gca,'FontSize',18)
%     set(gca,'linewidth',2)
%     set(gcf,'color','w');
% 
% ylabel('Max slip rate (m/s)')
% xlabel('Time (year)')
% % xlim([f1.time(end-200) f1.time(end)])

   % export_fig([save_name '_maxV2.pdf'])




%% PLot max P evolution
close all
figure(3)
plot(f1.time-10*365.25*86400-10,max(f(10).P/1e6,[],2),'linewidth',2)
 xlim([0 10000])
hold on



% Pressure
k = 1e-15;
beta = 1e-8;
phi = 0.1;
T = 10000.0;
Pmax = 2e6;
eta = 1e-3;
Q0 = Pmax*sqrt(k*beta*phi)/sqrt(eta*T)

% figure()
time = 0:T*3;
plot(time,Q0 *sqrt(eta*time)/sqrt(k*beta*phi)/1e6)
% xlabel('Time (s)')
% ylabel('Pressure (MPa)')

set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

ylabel('Max Pressure (MPa)')
xlabel('Time (s)')


%% Plot geometry alone
close all
offset = 0.7;

for fault_id = 1:length(input)
    hold on
    plot(f(fault_id).node(1,:),f(fault_id).node(2,:),'k-','LineWidth',2)

    % Plot fault number
    x_mean = (f(fault_id).node(1,1)+f(fault_id).node(1,end))/2;
    y_mean = (f(fault_id).node(2,1)+f(fault_id).node(2,end))/2;
    y_dist = max(f(fault_id).node(2,:))-min(f(fault_id).node(2,:));
    x_dist = max(f(fault_id).node(1,:))-min(f(fault_id).node(1,:));
    angle = 180/pi*atan(y_dist/x_dist);
    n1 = -y_dist/sqrt(x_dist^2+y_dist^2);
    n2 = x_dist/sqrt(x_dist^2+y_dist^2);
    text(x_mean+offset*n1,y_mean+offset*n2,['Fault ' num2str(fault_id)],'Rotation',angle,'FontSize',15,'HorizontalAlignment','center')
end

% Plot injection point
plot(f(10).element(1,input(10).index_injection),f(10).element(2,input(10).index_injection),'or','MarkerFaceColor','r','MarkerSize',10)



axis equal
xlim([-6,11])
ylim([-3 17])
xlabel('Position (m)')
ylabel('Position (m)')


set(gca,'FontSize',18)
set(gcf,'color','w');
export_fig([save_name '_geometry.pdf'])

%% Injected volume vs total moment
% Calculate size of injection
id_inj = input(10).index_injection;
L_injection = sqrt((input(10).node(1,id_inj+1)-input(10).node(1,id_inj))^2+(input(10).node(2,id_inj+1)-input(10).node(2,id_inj))^2);

% Find id time beg and end injection
[~,id_beg] = min(abs(f(1).time-input(10).t_injection_beg));
[~,id_end] = min(abs(f(1).time-input(10).t_injection_end));

% Calculate volume of injection
V = zeros(1,length(f(1).time));
V(id_beg:id_end) = cumtrapz(f(10).time(id_beg:id_end),input(10).Q*L_injection*ones(size(f(10).time(id_beg:id_end))));
V(id_end:end) = V(id_end);



% Calculate slip and ds
for fault_id = 1:length(f)
f(fault_id).slip = cumtrapz(f(fault_id).time,f(fault_id).V);
f(fault_id).ds = sqrt((f(fault_id).node(1,2:end)-f(fault_id).node(1,1:end-1)).^2 ...
                     +(f(fault_id).node(2,2:end)-f(fault_id).node(2,1:end-1)).^2);

% Calculate moment
f(fault_id).M0 = zeros(1,length(f(1).time));
for t_id =1:length(f(1).time)
    f(fault_id).M0(t_id) = sum(input(fault_id).mu*f(fault_id).slip(t_id,:).*f(fault_id).ds);
end
end




% Calculate moment 
offset_end = 0;
figure
fault_name = []
for fault_id=1:length(f)
plot(V(id_beg:id_end+offset_end)*1e3,(f(fault_id).M0(id_beg:id_end+offset_end)-f(fault_id).M0(id_beg))/1e6,'LineWidth',2)
hold on
fault_name = [fault_name;['Fault ' num2str(fault_id,'%.2i')]]
end





xlabel('Injected Volume (L)')
ylabel('M_0 (MPa.m^2)')


set(gca,'FontSize',18)
set(gcf,'color','w');

legend(fault_name,'Location','West')

export_fig([save_name '_injection_moment.pdf'])


%% Moment vs time
color = colororder('gem12');

% Calculate cumulative moment
cumulative_moment = zeros(1,length(f(1).time));
for fault_id = 1:length(f)
    cumulative_moment = cumulative_moment + f(fault_id).M0;
end



figure('Position',[1 1 1000 1000])
fault_name = []
yyaxis left
for fault_id=1:length(f)
plot(f(1).time/(86400*365.25),f(fault_id).M0/1e6,'-','LineWidth',2,'Color',color(fault_id,:))
hold on
fault_name = [fault_name;['Fault ' num2str(fault_id,'%.2i')]]
end
ylabel('Individual moment M_0 (MPa.m^2)')
yyaxis right
plot(f(1).time/(86400*365.25),cumulative_moment/1e6,'k','LineWidth',4)
ylabel('Total moment M_0 (MPa.m^2)')




xlabel('Time (year)')

ax = gca;
ax.YAxis(1).Color = 'k';
ax.YAxis(2).Color = 'k';


set(gca,'FontSize',18)
set(gcf,'color','w');

legend(fault_name,'Location','West')

export_fig([save_name '_moment_time.pdf'])


%% Moment vs time zoom

figure('Position',[1 1 1000 1000])
for fault_id=1:length(f)
plot(f(1).time-(10*86400*365.25),(f(fault_id).M0-f(fault_id).M0(id_beg))/1e6,'-','LineWidth',2,'Color',color(fault_id,:))
hold on
end

plot(f(1).time-10*(86400*365.25),(cumulative_moment-cumulative_moment(id_beg))/1e6,'k','LineWidth',4)
ylabel('Moment M_0 (MN)')
xlim([0 10000])




xlabel('Time (s)')




set(gca,'FontSize',18)
set(gcf,'color','w');

legend(fault_name,'Location','East')
export_fig([save_name '_moment_time_zoom.pdf'])



%% Time vs hydraulic energy
% Calculate injection history
f(10).Q = zeros(length(f(1).time),1);
for time_id=1:length(f(1).time)
    if f(10).time(time_id)>=input(10).t_injection_beg&& f(10).time(time_id)<input(10).t_injection_end
        f(10).Q(time_id) = input(10).Q;
    else
        f(10).Q(time_id) = 0.0;
    end
end

hy_energy = cumtrapz(f(10).time,f(10).P(:,input(10).index_injection).*f(10).Q*f(10).ds(input(10).index_injection));



figure('Position',[1 1 1000 1000])
yyaxis left

plot(f(1).time-10*365.25*86400,hy_energy,'LineWidth',2)
ylabel('Hydraulic energy (J)')

yyaxis right
plot(f(1).time-10*365.25*86400,f(10).P(:,input(10).index_injection)/1e6,'LineWidth',2)


xlabel('Time (s)')
ylabel('Pressure at injection point (MPa)')


set(gca,'FontSize',18)
set(gcf,'color','w');
xlim([0 10000])
title('Hydraulic energy')

export_fig([save_name '_hy_energy_time.pdf'])


%% Cumulative moment vs hydraulic energy




figure('Position',[1 1 1000 1000])

plot(hy_energy,cumulative_moment,'LineWidth',2)
ylabel('Total moment (MN)')
xlabel('Hydraulic energy (J)')

set(gca,'FontSize',18)
set(gcf,'color','w');
% xlim([0 10000])
title('Hydraulic energy')

export_fig([save_name '_hy_energy_cumulative.pdf'])


%% Zoom

figure('Position',[1 1 1000 1000])

plot(hy_energy(id_beg:id_end+10),(cumulative_moment(id_beg:id_end+10)-cumulative_moment(id_beg))/1e6,'LineWidth',2)
ylabel('Total moment (MN)')
xlabel('Hydraulic energy (J)')

set(gca,'FontSize',18)
set(gcf,'color','w');
% xlim([0 10000])

export_fig([save_name '_hy_energy_cumulative_zoom.pdf'])


%%
input(10).Q
V = input(10).Q*(input(10).t_injection_end-input(10).t_injection_beglat)




%%
plot(f(8).P(time_index(end),:))
hold on 
plot(f(9).P(time_index(end),:))
plot(f(7).P(time_index(end),:))
plot(f(6).P(time_index(end),:))
% plot(f(9).P(time_index(end),:))
% plot(f(9).P(time_index(end),:))
% plot(f(9).P(time_index(end),:))
% plot(f(9).P(time_index(end),:))

%%
clc
% For big fault
Tt = -input(1).sigmaN(1)*(input(1).f0(1) + input(1).a(1) * log(f(1).V(1,1)/input(1).V0(1))) ...
                                  + input(1).b(1) * log(f(1).theta(1,1)*input(1).V0(1)/input(1).Dc(1))



input(2).sigmaN(1)*input(2).f0(1)


%%  Quantification of \tau

T = (-input(1).sigmaN(1)*input(1).f0(1)-Tt)/(input(1).f0(1)*max(max(f(10).P)))