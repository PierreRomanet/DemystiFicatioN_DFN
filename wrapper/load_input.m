function f = load_input(direc)
% loadandprocessdata  reads the data from the given directory from the output 
% of FastCycles 2D
%
%   f = loadandprocessdata(direct) reads the velocity V of size time_step x
%   space_points
%
%   f = loadandprocessdata(direc,'St') reads the elastic shear traction of 
%   size time_step x space_points.
%
%   The different option are 
%   the velocity: 'V'
%   the state variable: 'theta'
%   the elastic shear traction: 'St'
%   the elastic normal traction: 'Sn'
%
%   f = loadandprocessdata(direc,'V','St') reads the velocity and the 
%   elastic shear traction, both of size time_step x space_points 

% Read all files name in the directory



% Get the number of faults
temp = h5info([direc 'config']);
% Initialisation counting number of faults
nb_fault = 0;
for group_id = 1:length(temp.Groups)
    if contains(temp.Groups(group_id).Name,'fault')
        nb_fault = nb_fault+1;
    end
end


% For each fault read the parameters
for fault_id=1:nb_fault

    % Create name of the fault 
    fault_name = ['/fault' num2str(fault_id)];

    % Friction
    f(fault_id).a = h5read([direc 'config'],[fault_name '/a']);
    f(fault_id).b = h5read([direc 'config'],[fault_name '/b']);
    f(fault_id).Dc = h5read([direc 'config'],[fault_name '/Dc']);
    f(fault_id).f0 = h5read([direc 'config'],[fault_name '/f0']);
    f(fault_id).V0 = h5read([direc 'config'],[fault_name '/V0']);
    f(fault_id).mu =  h5read([direc 'config'],'/material_loading/mu');


    % Geometry
    f(fault_id).node = h5read([direc 'config'],[fault_name '/node']);

    % Fluid injection
    f(fault_id).nb_source = h5read([direc 'config'],[fault_name '/nb_source']);
    if     f(fault_id).nb_source >=1 
            f(fault_id).index_injection = h5read([direc 'config'],[fault_name '/index_injection']);
            f(fault_id).t_injection_beg =  h5read([direc 'config'],[fault_name '/t_injection_beg']);
            f(fault_id).t_injection_end =  h5read([direc 'config'],[fault_name '/t_injection_end']);
            f(fault_id).Q =  h5read([direc 'config'],[fault_name '/Q']);
            
    end
    
    % Hydraulic model
    f(fault_id).porosity =  h5read([direc 'config'],'/material_loading/porosity');
    f(fault_id).rock_comp =  h5read([direc 'config'],'/material_loading/rock_comp');
    f(fault_id).fluid_comp =  h5read([direc 'config'],'/material_loading/fluid_comp');
    f(fault_id).dyn_viscosity =  h5read([direc 'config'],'/material_loading/dyn_viscosity');
        f(fault_id).fluid_density =  h5read([direc 'config'],'/material_loading/fluid_density');

        f(fault_id).cp =  h5read([direc 'config'],'/material_loading/cp');
    f(fault_id).cs =  h5read([direc 'config'],'/material_loading/cs');

    f(fault_id).permeability =  h5read([direc 'config'],[fault_name '/permeability']);

    % 
    f(fault_id).sigmaN =  h5read([direc 'config'],[fault_name '/sigmaN']);



end