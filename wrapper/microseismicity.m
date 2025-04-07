%% Load data
clear all
clc
close all
addpath '/Users/pierre/Dropbox/PapersinPreparation/Romanet2025_tripleJunction/matlab'



% Name of the file
save_name = 'LSBB_with_800faults_noLoading_D7000';% Most important simulation%

% Path of all the original simulation of the paper
%direc0 = ['/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_v04/problems/LSBB_with_800faults_noLoading1_paper/'];



% Path of all the simulation
direc1 = ['/Users/pierre/Dropbox/OnGoingResearch/Fluid_and_earthquakes/Softwares/DemystiFicatioN_v07/problems/' save_name '/'];




% Base level
base_level = 1e-15;
max_level = 1e0;
skip = 20;
unit = 1.0; % meter
unit_P = 1e6;

% Load data
f1 = loadandprocessdata(direc1,'V','P','slip','theta','Sn');
input = load_input(direc1);
% input0 = load_input(direc0);
%
% Reorganise data according to each fault
id_beg = 1;
for fault_id = 1:length(input)

    % Set id_end
    id_end = id_beg + length(input(fault_id).a)-1;

    % Create structure
    f(fault_id).V = f1.V(:,id_beg:id_end);
    f(fault_id).theta = f1.theta(:,id_beg:id_end);
    f(fault_id).P = f1.P(:,id_beg:id_end);
    f(fault_id).Sn = f1.Sn(:,id_beg:id_end);
    f(fault_id).node = input(fault_id).node;
    f(fault_id).element = (f(fault_id).node(:,2:end)+f(fault_id).node(:,1:end-1))/2;
    f(1).time = f1.time;

    % Calculate element
    f(fault_id).element =  (f(fault_id).node(:,2:end)+f(fault_id).node(:,1:end-1))/2;



    id_beg = id_end + 2;

end


figure
plot(f(1).V(end,:))


%%
% Calculate Coulomb stress change due to the main fault
close all
% Create element 
x=-10.1:0.2:10;
y =-4.1:0.2:4;
[X,Y] = meshgrid(x,y);

element = zeros(2,size(X,1)*size(X,2));
element(1,:) = reshape(X,1,[]);
element(2,:) = reshape(Y,1,[]);


% Calculate slip 
slip = cumtrapz(f(1).time,f(1).V);

% Calculate stresses
scaling = 32.04e9/pi*(1-3.464e3/5000.0);
[S11_total,S22_total,S12_total,S11_gradient,S22_gradient,S12_gradient]=calculate_modeII(f(1).node,element,scaling);

% Calculate Coulomb change
coulomb = (S12_total+0.6*S22_total)*slip(2000,:).';
coulomb = reshape(coulomb,length(y),length(x));

% pcolor
% pcolor(X,Y,coulomb)








%%

time = 86400*365.25*10+86400*2;
[~,time_id] = min(abs(f(1).time-time)); 

% Plot data 
% Close all figures
close all
figure('position', [0 0 1400 900])


% create the video writer with 1 fps
writerObj = VideoWriter(['myVideo_' save_name] ,'MPEG-4');
writerObj.FrameRate = 10;
writerObj.Quality = 100;
% open the video writer
open(writerObj)

fault_index_array = zeros(length(input),1);
for i=1:skip:time_id%length(f1.time)% [500 1000 1500 10000]
    clf

    ax1 =axes;



    % Make plot velocity
    h = plot3(f(1).element(1,:)/unit,f(1).element(2,:)/unit,squeeze(f(1).V(i,:)),'linewidth',3);
    hold on


    % Make surface 
    % Calculate Coulomb change
    coulomb = (S12_total+0.6*S22_total)*slip(i,:).';
    coulomb = reshape(coulomb,length(y),length(x));
    zl = zlim;
    Z = base_level*ones(size(X));
    surf(ax1,X,Y,Z,coulomb/1e6);

    colormap(redblue(100))
    clim([-1 1])

    cb1 = colorbar;
    cb1.Label.String = 'Coulomb stress change (MPa)';
    cb1.Label.FontSize = 20;
    cb1.LineWidth = 2;
    cb1.Location = 'northoutside';
    a =  cb1.Position;
    set(cb1,'Position',[a(1)+a(3)*0.25 a(2)*0.9 a(3)*0.5 a(4)*0.6])
        set(ax1,'FontSize',18)

    shading interp
    % colorbar
    % Make plot fault

    % Link fault to velocity
    plot3([f(1).element(1,1)/unit f(1).element(1,1)/unit],[f(1).element(2,1)/unit f(1).element(2,1)/unit],[base_level*2  f(1).V(i,1)],'k','linewidth',2)
    plot3([f(1).element(1,end)/unit f(1).element(1,end)/unit],[f(1).element(2,end)/unit f(1).element(2,end)/unit],[base_level*2  f(1).V(i,end)],'k','linewidth',2)


    % Plot all the smaller faults
    % Make plot velocity
    for fault_id = 2:length(input)


        plot3(ax1,f(fault_id).node(1,:)/unit,f(fault_id).node(2,:)/unit,base_level*2*ones(size(f(fault_id).node(2,:))),'k','linewidth',3);


        if fault_index_array(fault_id) == 1
                plot3(f(fault_id).node(1,2)/unit,f(fault_id).node(2,2)/unit,base_level*2.5,'k','linewidth',1,'Marker','pentagram','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',15);
        end

        % Make a start in case of earthquake
        if max(f(fault_id).V(i,:))>1e-3
             plot3(f(fault_id).node(1,2)/unit,f(fault_id).node(2,2)/unit,base_level*2.5,'k','linewidth',3,'Marker','pentagram','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',25);
             fault_index_array(fault_id) = 1;
        end




        %h = plot3(f(fault_id).element(1,:)/unit,f(fault_id).element(2,:)/unit,squeeze(f(fault_id).V(i,:)),'linewidth',3);

        % Link fault to velocity
        %plot3([f(fault_id).element(1,1)/unit f(fault_id).element(1,1)/unit],[f(fault_id).element(2,1)/unit f(fault_id).element(2,1)/unit],[base_level*2  f(fault_id).V(i,1)],'k','linewidth',2)
        %plot3([f(fault_id).element(1,end)/unit f(fault_id).element(1,end)/unit],[f(fault_id).element(2,end)/unit f(fault_id).element(2,end)/unit],[base_level*2  f(fault_id).V(i,end)],'k','linewidth',2)
    end



    % plot pressure
    % surface([f(1).element(1,:)/unit;f(1).element(1,:)/unit],...
    %     [f(1).element(2,:)/unit;f(1).element(2,:)/unit], ...
    %     [ones(size(f(1).element(2,:)))*base_level*2;ones(size(f(1).element(2,:)))*base_level*2],...
    %     [f(1).Sn(i,:)/unit_P;f(1).Sn(i,:)/unit_P],...
    %     'facecol','no',...
    %     'edgecol','interp',...
    %     'linewidth',12);






    % Plot time
    text(-12,0,1e-5,seconds2duration(f1.time(i)),'FontSize',20)

    maxVel = max(f(1).V(i,:));
    % Plot maximum velocity
    text(-12,0,1e-7,{'Maximum slip rate :', [num2str(maxVel,'%2.5e') ' m/s']},'FontSize',20)



    xlim([-10 10])
    ylim([-4 4])
    zlim([base_level  max_level])


     % Improve z axis
      % clb.Label.String = 'Slip velocity (m/s)';

    set(gca,'zscale','log')
    set(gca,'ZMinorTick','off')




    xlabel('X position (m)')
    ylabel('Y position (m)')
    zlabel('Slip rate (m/s)')



    set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

    ax2 = axes
    surface(ax2,[f(1).element(1,:)/unit;f(1).element(1,:)/unit],...
        [f(1).element(2,:)/unit;f(1).element(2,:)/unit], ...
        [ones(size(f(1).element(2,:)))*base_level*2;ones(size(f(1).element(2,:)))*base_level*2],...
        [f(1).P(i,:)/unit_P;f(1).P(i,:)/unit_P],...
        'facecol','no',...
        'edgecol','interp',...
        'linewidth',12);

    hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];

 colormap(ax2,jet)
        cb = colorbar(ax2);
    cb.Label.String = 'Pressure (MPa)';
    cb.Label.FontSize = 20;
    cb.LineWidth = 2;
    cb.Location = 'northoutside';
    a =  cb.Position;
    set(cb,'Position',[a(1)+a(3)*0.25 a(2)*0.8 a(3)*0.5 a(4)*0.6])
        set(ax2,'FontSize',18)

    % cb.Position = [.95 0.1 0.02 0.5];
    clim([0 1])

        view(20,35)


    writeVideo(writerObj, getframe(get(groot,'CurrentFigure')));

end
   % export_fig([save_name '_fig_poster5.png'])

close(writerObj);




%% PLot max V evolution
% close all
figure('Position',[1 1 1500 1000])

for fault_id = 1:length(input)
    semilogy(f1.time/(86400*365.25),max(f(fault_id).V,[],2),'k','linewidth',1)
 hold on

end

semilogy(f1.time/(86400*365.25),max(f(1).V,[],2),'r','linewidth',2)
% xlim([0 100])


set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

ylabel('Max slip rate (m/s)')
xlabel('Time (year)')

export_fig([save_name '_maxV.png'])
xlim([10 f1.time(end)/(86400*365.25)])
 export_fig([save_name '_ZoommaxV.png'])


figure('Position',[1 1 1500 1000])

% for fault_id = 1:length(input)
%     semilogy(f1.time/(86400*365.25),max(f(fault_id).V,[],2),'k','linewidth',1)
%  hold on
% 
% end

semilogy(f1.time/(86400*365.25),max(f(1).V,[],2),'r','linewidth',2)
% xlim([0 100])


set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

ylabel('Max slip rate (m/s)')
xlabel('Time (year)')

% export_fig([save_name '_maxV.png'])
xlim([10 f1.time(end)/(86400*365.25)])
 export_fig([save_name '_ZoommaxV.pdf'])




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
figure(3)
plot(f1.time/(86400*365.25),max(f(1).P/1e6,[],2),'linewidth',2)
 



set(gca,'FontSize',18)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

ylabel('Max Pressure (MPa)')
xlabel('Time (year)')


   export_fig([save_name '_maxV.pdf'])


%% Pressure
k = 1e-15;
beta = 1e-8;
phi = 0.1;
T = 1000.0;
Pmax = 2e6;
eta = 1e-3;
Q0 = Pmax*sqrt(k*beta*phi)/sqrt(eta*T)

figure()
time = 0:T*3;
plot(time,Q0 *sqrt(eta*time)/sqrt(k*beta*phi)/1e6)
xlabel('Time (s)')
ylabel('Pressure (MPa)')


%% Calculate Coulomb stress change due to the main fault
close all
% Create element 
x=-10.1:0.2:10;
y =-4.1:0.2:4;
[X,Y] = meshgrid(x,y);

element = zeros(2,size(X,1)*size(X,2));
element(1,:) = reshape(X,1,[]);
element(2,:) = reshape(Y,1,[]);


% Calculate slip 
slip = cumtrapz(f(1).time,f(1).V);

% Calculate stresses
scaling = 32.04e9/pi*(1-3.464e3/5000.0);
[S11_total,S22_total,S12_total,S11_gradient,S22_gradient,S12_gradient]=calculate_modeII(f(1).node,element,scaling);

% Calculate Coulomb change
coulomb = (S12_total+0.6*S22_total)*slip(2000,:).';
coulomb = reshape(coulomb,length(y),length(x));

% pcolor
pcolor(X,Y,coulomb)
colormap(redblue(100))
clim([-1e6 1e6])
shading interp
colorbar


%%
figure
% plot(f(1).node(1,:),f(1).node(2,:))
% plot(element(1,:),element(2,:),'o')
plot(f(1).element(1,:),slip(2000,:))


%% Find distance of microseismicity
% Build a catalog of earthquakes
V_threshold = 1e-3;
catalog = [] ;
EQ_on_going = 0;
nb_EQ = 0;
fault_ID_EQ_ongoing = 1;
for time_id =1:length(f(1).time)
    for fault_id=2:length(f)
        if (max(f(fault_id).V(time_id,:))>=V_threshold)&& EQ_on_going ==0
            nb_EQ = nb_EQ + 1;
            EQ_on_going = 1;
            fault_ID_EQ_ongoing = fault_id;
            temp = zeros(1,3);
            temp(1,1) = fault_id;
            temp(1,2) = time_id;
        end





    end
    if (max(f(fault_ID_EQ_ongoing).V(time_id,:))<=V_threshold)&& EQ_on_going ==1
        EQ_on_going = 0;
        temp(1,3) = time_id;
        catalog = vertcat(catalog,temp);
    end
end

% Calculate distance from earthquake
dist = zeros(size(catalog,1),1);
% Loop over the catalog
for eq_id=1:size(catalog,1)

    dist1 = sqrt((f(catalog(eq_id,1)).element(1,1)-f(1).element(1,1))^2 + ...
                 (f(catalog(eq_id,1)).element(2,1)-f(1).element(2,1))^2);
    dist2 = sqrt((f(catalog(eq_id,1)).element(1,1)-f(1).element(1,end))^2 + ...
                 (f(catalog(eq_id,1)).element(2,1)-f(1).element(2,end))^2);
    dist(eq_id) = min(dist1,dist2);

end
%%
figure()% Pressure
k = 8e-17%input(1).permeability(1);
beta = input(1).fluid_comp+input(1).rock_comp;
phi = input(1).porosity;
T = 10000.0;
Pmax = 2e6;
eta = input(1).dyn_viscosity;
t_max = 180000;



t = 0:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);





% 
plot(t/60,dist_th,'b','linewidth',2)
hold on
k = input(1).permeability(1);
t = 0:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);
plot(t/60,dist_th,'b','linewidth',2)


text(1500,2.95,'k = 8\times 10^{-17} m^2','Rotation',18,'FontSize',18,'color','b')
text(300,3.7,'k = 1\times 10^{-15} m^2','Rotation',65,'FontSize',18,'color','b')
plot([input(1).t_injection_end-86400*10*365.25 input(1).t_injection_end-86400*10*365.25]/60,[0 6],'k--','LineWidth',1)
text(150,3.7,'End of injection','Rotation',90,'FontSize',18,'color','k')

scatter((f(1).time(catalog(:,2))-86400*10*365.25)/60,dist,"pentagram",'MarkerEdgeColor','k','MarkerFaceColor','r','SizeData',100)
ylim([0 6])
xlim([0 (f(1).time(end)-86400*10*365.25)/60])
% xlim([0 10000])
% xlim([0 100000])
xlim([0 3000])
set(gca,'FontSize',18)
set(gcf,'color','w');
xlabel('Time (min)')
ylabel('Distance from main fault edges (m)')
export_fig([save_name '_diffusion.pdf'])



%% Make snapshop 
% Plot data 
% Close all figures
figure('position', [0 0 1400 900])
%[104 130 350 1250 2577 15000]


time_sp = [86400*365.25*10-5*60,86400*365.25*10+20*60,86400*365.25*10+30*60,input(1).t_injection_end,86400*365.25*10+4*60*60,86400*365.25*10+24*60*60];

% Find closest id to time
time_index= zeros(size(time_sp));
for i=1:length(time_sp)
[~,time_index(i)] = min(abs(f(1).time-time_sp(i)));
time_index(i) = time_index(i)+1;
end



for i=time_index
    clf
    ax1 =axes;
    


    % Make plot velocity
    h = plot3(f(1).element(1,:)/unit,f(1).element(2,:)/unit,squeeze(f(1).V(i,:)),'linewidth',3);
    hold on
   

    % Make surface 
    % Calculate Coulomb change
    coulomb = (S12_total+0.6*S22_total)*slip(i,:).';
    coulomb = reshape(coulomb,length(y),length(x));
    zl = zlim;
    Z = base_level*ones(size(X));
    surf(ax1,X,Y,Z,coulomb/1e6);

    colormap(redblue(100))
    clim([-1 1])

    cb1 = colorbar;
    cb1.Label.String = 'Coulomb stress change (MPa)';
    cb1.Label.FontSize = 22;
    cb1.LineWidth = 2;
    % cb1.Location = 'northoutside';
    cb1.Visible = 'off'
    a =  cb1.Position;
    % set(cb1,'Position',[a(1)+a(3)*0.25 a(2)*0.9 a(3)*0.5 a(4)*0.6])
        set(ax1,'FontSize',22)

    shading interp
    % colorbar
    % Make plot fault

    % Link fault to velocity
    plot3([f(1).element(1,1)/unit f(1).element(1,1)/unit],[f(1).element(2,1)/unit f(1).element(2,1)/unit],[base_level*2  f(1).V(i,1)],'k','linewidth',2)
    plot3([f(1).element(1,end)/unit f(1).element(1,end)/unit],[f(1).element(2,end)/unit f(1).element(2,end)/unit],[base_level*2  f(1).V(i,end)],'k','linewidth',2)


    % Plot all the smaller faults
    % Make plot velocity
    for fault_id = 2:length(input)


        plot3(ax1,f(fault_id).node(1,:)/unit,f(fault_id).node(2,:)/unit,base_level*2*ones(size(f(fault_id).node(2,:))),'k','linewidth',3);

    end

    % Plot EQ
    for EQ_id=1:size(catalog,1)
        if i>=catalog(EQ_id,2)&& i<catalog(EQ_id,2)
            plot3(f(catalog(EQ_id,1)).node(1,2)/unit,f(catalog(EQ_id,1)).node(2,2)/unit,base_level*2.5,'k','linewidth',3,'Marker','pentagram','MarkerFaceColor','y','MarkerEdgeColor','k','MarkerSize',25);
        elseif i>=catalog(EQ_id,2)
            plot3(f(catalog(EQ_id,1)).node(1,2)/unit,f(catalog(EQ_id,1)).node(2,2)/unit,base_level*2.5,'k','linewidth',1,'Marker','pentagram','MarkerFaceColor','r','MarkerEdgeColor','k','MarkerSize',15);
        end

    end



        %h = plot3(f(fault_id).element(1,:)/unit,f(fault_id).element(2,:)/unit,squeeze(f(fault_id).V(i,:)),'linewidth',3);

        % Link fault to velocity
        %plot3([f(fault_id).element(1,1)/unit f(fault_id).element(1,1)/unit],[f(fault_id).element(2,1)/unit f(fault_id).element(2,1)/unit],[base_level*2  f(fault_id).V(i,1)],'k','linewidth',2)
        %plot3([f(fault_id).element(1,end)/unit f(fault_id).element(1,end)/unit],[f(fault_id).element(2,end)/unit f(fault_id).element(2,end)/unit],[base_level*2  f(fault_id).V(i,end)],'k','linewidth',2)

    fontsize = 30;

    % Plot time 
    if f(1).time(i)-86400*365.25*10>0
        text(-12,0,1e-5,[seconds2duration_special(f(1).time(i)-86400*365.25*10) ' after injection'],'FontSize',fontsize)
    else
        text(-12,0,1e-5,[seconds2duration_special(abs(f(1).time(i)-86400*365.25*10)) ' before injection'],'FontSize',fontsize)
    end
    % text(-12,0,1e-5,seconds2duration_special(f1.time(i)),'FontSize',20)

    maxVel = max(f(1).V(i,:));
    % Plot maximum velocity
    text(-12,0,1e-7,{'Maximum slip rate :', [num2str(maxVel,'%2.5e') ' m/s']},'FontSize',fontsize)


     
    xlim([-10 10])
    ylim([-4 4])
    zlim([base_level  max_level])


     % Improve z axis
      % clb.Label.String = 'Slip velocity (m/s)';

    set(gca,'zscale','log')
    set(gca,'ZMinorTick','off')




    xlabel('X position (m)')
    ylabel('Y position (m)')
    zlabel('Slip rate (m/s)')

   

    set(gca,'FontSize',fontsize)
    set(gca,'linewidth',2)
    set(gcf,'color','w');
    
    ax2 = axes
    surface(ax2,[f(1).element(1,:)/unit;f(1).element(1,:)/unit],...
        [f(1).element(2,:)/unit;f(1).element(2,:)/unit], ...
        [ones(size(f(1).element(2,:)))*base_level*2;ones(size(f(1).element(2,:)))*base_level*2],...
        [f(1).P(i,:)/unit_P;f(1).P(i,:)/unit_P],...
        'facecol','no',...
        'edgecol','interp',...
        'linewidth',12);

    hLink = linkprop([ax1,ax2],{'XLim','YLim','ZLim','CameraUpVector','CameraPosition','CameraTarget'});
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];

 colormap(ax2,jet)
        cb = colorbar(ax2);
    cb.Label.String = 'Pressure (MPa)';
    cb.Label.FontSize = 22;
    cb.LineWidth = 2;
    % cb.Location = 'northoutside';
    cb.Visible = 'off'
    a =  cb.Position;
    % set(cb,'Position',[a(1)+a(3)*0.25 a(2)*0.8 a(3)*0.5 a(4)*0.6])
        set(ax2,'FontSize',22)

    % cb.Position = [.95 0.1 0.02 0.5];
    clim([0 1])

        view(20,35)
pause(1)
 export_fig([save_name '_' num2str(i) '_fig.png'])


end





%%
figure


ax1 =axes;
 colormap(redblue(100))
    clim([-1 1])

    cb1 = colorbar;
    cb1.Label.String = 'Coulomb stress change (MPa)';
    cb1.Label.FontSize = 22;
    cb1.LineWidth = 2;
    cb1.Location = 'northoutside';
    % cb1.Visible = 'off'
    a =  cb1.Position;
    set(cb,'Position',[a(1) a(2)*1.5 a(3) a(4)])
        set(ax1,'FontSize',22)

    export_fig([save_name '_colorbar1.pdf'])


        figure
    ax2 = axes

   

  cb = colorbar(ax2);
colormap(jet)
    cb.Label.String = 'Pressure (MPa)';
    cb.Label.FontSize = 22;
    cb.LineWidth = 2;
    cb.Location = 'northoutside';
    % cb.Visible = 'off'
    a =  cb.Position;
    set(cb,'Position',[a(1) a(2)*0.8 a(3) a(4)])
        set(ax2,'FontSize',22)

    % cb.Position = [.95 0.1 0.02 0.5];
    clim([0 1]) 
    
    export_fig([save_name '_colorbar2.pdf'])


    %% Geometry enhanced

    figure('Position',[1 1 1000 1000])
    plot(f(1).node(1,:),f(1).node(2,:),'LineWidth',2)


ylim([-0.25 0.25])
xlim([-5 5])


    set(gca,'FontSize',22)
    set(gca,'linewidth',2)
    set(gcf,'color','w');

export_fig([save_name '_geometryMainFault.pdf'])


%% Calculate volume and flow rate
f(1).ds = sqrt((f(1).node(1,2:end)-f(1).node(1,1:end-1)).^2+(f(1).node(2,2:end)-f(1).node(2,1:end-1)).^2);
Q = input(1).Q
V = Q * (input(1).t_injection_end-input(1).t_injection_beg)


%%
figure('Position',[1 1 1500 600])% Pressure
tiledlayout(1,2)
size_font = 25
nexttile
k = 1e-16%input(1).permeability(1);
beta = input(1).fluid_comp+input(1).rock_comp;
phi = input(1).porosity;
T = 10000.0;
Pmax = 2e6;
eta = input(1).dyn_viscosity;
t_max = 180000;



t = 0:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);






plot(t/60,dist_th,'b','linewidth',2)
hold on
k = input(1).permeability(1);
t = 0:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);
plot(t/60,dist_th,'b','linewidth',2)


text(1500,3.25,'k = 1\times 10^{-16} m^2','Rotation',24,'FontSize',size_font,'color','b')
text(195,2.8,'k = 1\times 10^{-15} m^2','Rotation',75,'FontSize',size_font,'color','b')
plot([input(1).t_injection_end-86400*10*365.25 input(1).t_injection_end-86400*10*365.25]/60,[0 6],'k--','LineWidth',1)

scatter((f(1).time(catalog(:,2))-86400*10*365.25)/60,dist,"pentagram",'MarkerEdgeColor','k','MarkerFaceColor','r','SizeData',175)
text(150,1.1,'End of injection','Rotation',90,'FontSize',size_font,'color','k')

ylim([0 4.5])
set(gca,'FontSize',size_font)
set(gcf,'color','w');
xlabel('Time (min)')
ylabel('Distance from main fault edges (m)')

nexttile

k = 10e-17%input(1).permeability(1);
beta = input(1).fluid_comp+input(1).rock_comp;
phi = input(1).porosity;
T = 10000.0;
Pmax = 2e6;
eta = input(1).dyn_viscosity;
t_max = 180000;



t = 1:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);






plot(log10(t),dist_th,'b','linewidth',2)
hold on
k = input(1).permeability(1);
dist_th = sqrt(k/(beta*eta*phi)*t);
plot(log10(t),dist_th,'b','linewidth',2)


text(log10(1e4),1.2,'k = 1\times 10^{-16} m^2','Rotation',34,'FontSize',size_font,'color','b')
text(log10(1e4),3.,'k = 1\times 10^{-15} m^2','Rotation',55,'FontSize',size_font,'color','b')
plot(log10([input(1).t_injection_end-86400*10*365.25 input(1).t_injection_end-86400*10*365.25]),[0 6],'k--','LineWidth',1)
text(log10(100*56),2.75,'End of injection','Rotation',90,'FontSize',size_font,'color','k')

scatter(log10(f(1).time(catalog(:,2))-86400*10*365.25),dist,"pentagram",'MarkerEdgeColor','k','MarkerFaceColor','r','SizeData',175)
ylim([0 4.5])
set(gca,'FontSize',size_font)
set(gcf,'color','w');

xl1 = [];
xl = xticklabels();
for i=1:length(xl)
    xl{i} = ['10^' xl{i}];
end
xticklabels(xl);
xlim([log10(1000) log10(2e5)])

xlabel('Time (seconds)')
ylabel('Distance from main fault edges (m)')
export_fig([save_name '_diffusion_log.pdf'])

%%

figure('Position',[1 1 800 600])% Pressure
size_font = 25

k = 1e-16%input(1).permeability(1);
beta = input(1).fluid_comp+input(1).rock_comp;
phi = input(1).porosity;
T = 10000.0;
Pmax = 2e6;
eta = input(1).dyn_viscosity;
t_max = 180000;



t = 0:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);


plot(t/60,dist_th,'b','linewidth',2)
hold on
k = input(1).permeability(1);
t = 0:1:t_max;
dist_th = sqrt(k/(beta*eta*phi)*t);
plot(t/60,dist_th,'b','linewidth',2)


text(1500,3.25,'\alpha = 10^{-4} m^2/s','Rotation',25,'FontSize',size_font,'color','b')
text(220,3.1,'\alpha = 10^{-3} m^2/s','Rotation',75,'FontSize',size_font,'color','b')
plot([input(1).t_injection_end-86400*10*365.25 input(1).t_injection_end-86400*10*365.25]/60,[0 6],'k--','LineWidth',1)

scatter((f(1).time(catalog(:,2))-86400*10*365.25)/60,dist,"pentagram",'MarkerEdgeColor','k','MarkerFaceColor','r','SizeData',175)
text(150,1.1,'End of injection','Rotation',90,'FontSize',size_font,'color','k')

ylim([0 4.5])
set(gca,'FontSize',size_font)
set(gcf,'color','w');
xlabel('Time (min)')
ylabel('Distance from main fault edges (m)')
xlim([0 3000])
export_fig([save_name '_diffusion2.pdf'])

%% Velocity of the center of the fault

figure
nb_half = floor(length(f(1).V(1,:))/2);
loglog(f(1).time-86400*10*365.25,mean(f(1).V(:,:),2),'linewidth',2)
xlabel('Time (s)')
ylabel('Mean velocity main fault (m/s)')
set(gca,'FontSize',size_font)
set(gcf,'color','w');
export_fig([save_name '_V_main_fault.pdf'])


% xlim([86400*10*365.25 86400*10*365.25+10000])

%% Calculate nucleation lengthscale
% Lb
Lb = -input(1).mu*input(1).Dc(1)*2*(1-input(1).cs^2/input(1).cp^2)/(input(1).sigmaN(1)*input(1).b(1))
input(1).node(2,2)
% Lnuc

mean(f(1).ds)*length(f(1).ds)

% Gridsize

%% Calculate slip 
f(1).slip = cumtrapz(f(1).time,f(1).V);

close all
plot((f(1).time-86400*10*365.25)/60,max(f(1).slip/1e-3,[],2),'linewidth',2)
xlim([0 200000/60])
xlabel('Time (min)')
ylabel('Maximum of slip (mm)')
set(gca,'FontSize',size_font)
set(gcf,'color','w');
export_fig([save_name '_max_slip.pdf'])


%% Calculate shear and normal traction on each fault
% Calculate shear traction on each fault
for i=1:length(f)
    f(i).Coulomb = -input(i).f0
end


%% Diffusion
alpha = 1e-16/(input(1).porosity*input(1).dyn_viscosity*(input(1).fluid_comp+input(1).rock_comp))
alpha = 1e-15/(input(1).porosity*input(1).dyn_viscosity*(input(1).fluid_comp+input(1).rock_comp))

%% Friction at which small faults are triggered
% Calculate friction for each fault
for fault_id = 2:length(input)

    % Create structure
    f(fault_id).friction = input(fault_id).f0 + ...
                           input(fault_id).a*log(f(fault_id).V/input(fault_id).V0) + ...
                           input(fault_id).b*log(f(fault_id).theta*input(fault_id).V0/input(fault_id).Dc);


end

figure('position',[1 1 1000 1000])
% For each EQ
for EQ_id=1:size(catalog,1)
    plot(f(1).time(catalog(EQ_id,2))-86400*10*365.25,f(catalog(EQ_id,1)).friction(catalog(EQ_id,2)),'ro','MarkerFaceColor','red')
    hold on


end
xlabel('Time (s)')
ylabel('Friction f')
title('Friction at the onset of micro-earthquake (threhold slip velocity 10^{-3} m/s)')
xlim([0 200000])

set(gca,'FontSize',20)
set(gcf,'color','w');
export_fig([save_name '_friction.pdf'])


%% Calculate theta evolution for each subfault

figure
for i=2:100%length(f)
plot(f(1).time -86400*10*365.25,f(i).theta(:,1))
hold on

end

xlim([0 f(1).time(end)-86400*10*365.25])

%% Time step
close all

semilogy(diff(f(1).time))




%%
for i=1:1:length(f(1).time)
plot(f(1).element(1,:),f(1).P(i,:).')
title(seconds2duration(f(1).time(i)))
set(gca,'FontSize',20)
pause(0.05)


end

%% Problem with pressure evolution

[~,time_id] = min(abs(f(1).time-input(1).t_injection_end));


figure('position',[1 1 1000 700])
plot(f(1).element(1,:),f(1).P(time_id:300:end,:)','Linewidth',2)


xlabel('Position (m)')
ylabel('Pressure (MPa)')
title('Pressure evolution after stopping the injection')
set(gca,'FontSize',20)
set(gcf,'Color','W')
export_fig('problem_LSBB_paper.png','-m3')


%%


plot(f(1).time,max(f(1).P,[],2))
xlim([10*365.25*86400 10*365.25*86400+90000])




%%
clc
% For big fault
Tt = -input(1).sigmaN(1)*(input(1).f0(1) + input(1).a(1) * log(f(1).V(1,1)/input(1).V0(1))) ...
                                  + input(1).b(1) * log(f(1).theta(1,1)*input(1).V0(1)/input(1).Dc(1))

% For small faults
Tt = -input(2).sigmaN(1)*(input(2).f0(1) + input(2).a(1) * log(f(2).V(1,1)/input(2).V0(1))) ...
                                  + input(2).b(1) * log(f(2).theta(1,1)*input(2).V0(1)/input(2).Dc(1))


input(2).sigmaN(1)*input(2).f0(1)



%% Max V just before injection
clf
[~,id] = min(abs(f(1).time-86400*10*365.25));

id = id -1
max(f(1).V(id,:))




semilogy(f(1).time,max(f(1).V,[],2))
hold on
semilogy(f(1).time(id),max(f(1).V(id,:)),'ro')


%%
% For big fault
Tt = -input(1).sigmaN(1)*(input(1).f0(1) + input(1).a(1) * log(f(1).V(1,1)/input(1).V0(1))) ...
                                  + input(1).b(1) * log(f(1).theta(1,1)*input(1).V0(1)/input(1).Dc(1))



input(2).sigmaN(1)*input(2).f0(1)


%%  Quantification of \tau

T = (-input(1).sigmaN(1)*input(1).f0(1)-Tt)/(input(1).f0(1)*max(max(f(1).P)))