function [time_all,xg_all,h_all,parameters] = GroundingLine_FlowlineModel(accum,A_glen,plotting)
%GroundingLine_FlowlineModel is a flowline model that calculates the
%time-dependent grounding line positions, thickness and velocity 
%for a marine ice sheet with a Weertman sliding law and bed topography from
%Schoof, JGR, 2007.
%The model accept accumulation rate (in meters/sec)
%and the Glen's  law rate factor as inputs and outputs the time step and
%time_dependent grounding line position and thickness. 
%If the plotting input is set to 1, the model will plot the flowline 
%thickness and velocity during calculation.
%This flowline model uses an operator-split integration that calculates 
%ice velocity, thickness & grounding line position using an implicit method.
%The details can be found in Robel et al., JGR, 2014

%This code has been prepared for public distribution by Alex Robel
%You may contact Alex at robel@caltech.edu for more information on this and
%other simple marine ice sheet models.

%% Parameters
%%Set Grid Resolution
parameters.grid.n_nodes = 2000;      %Horizontal Resolution
parameters.grid.n2_nodes = 40;
parameters.grid.gz_nodes = 400;      %Horizontal Resolution in grounding zone
parameters.grid.sigma_gz = 0.97;
parameters.Dupont_G = 0;          %lateral shear stress

parameters.year = 3600*24*365;                     %length of a year in seconds
parameters.tfinal = 10e3.*parameters.year;          %total time of integration
parameters.nsteps = 5e3;                           %number of time steps

%%Time step parameters
parameters.dtau = parameters.tfinal/parameters.nsteps; %length of time steps
parameters.dtau_max = parameters.dtau;

%%Newton method Parameters
parameters.HS_sensitivity = pi*parameters.year;     %sensitivity of the HS function (as this gets larger, theta approaches the actual HS function)
parameters.uverbose = 1;
parameters.iteration_threshold = 1e-3;
parameters.hiter_max=1e3;
parameters.uiter_max=5e2;
parameters.titer_max=4e1;
parameters.CFL=50;

%%Grid Parameters
parameters.grid.n_elements = parameters.grid.n_nodes-1;           %number of finite elements (= n_nodes-1 in 1-D), h and N have length n_elements
parameters.grid.sigma_node = [linspace(0,0.97,parameters.grid.n_nodes-parameters.grid.gz_nodes),linspace(0.97+(.03/parameters.grid.gz_nodes),1,parameters.grid.gz_nodes)]'; %node positions scaled to (0,1) with refinement near GL
parameters.grid.sigma_element =...
    (parameters.grid.sigma_node(1:parameters.grid.n_nodes-1)+...
    parameters.grid.sigma_node(2:parameters.grid.n_nodes))/2;     %element centres scaled to (0,1)

%%Glen's Law parameters
parameters.B_Glen = (A_glen^(-1/3)).* ones(parameters.grid.n_elements,1);                     %B in Glen's law (vertically averaged if necessary)
parameters.n_Glen = 3;

%%Physical parameters
parameters.rho = 900;  %917                                 %ice density
parameters.rho_w = 1000;  %1028                               %water density
parameters.g = 9.81;                                    %acceleration due to gravity
parameters.D_eps = 1e-10;                               %strain rate regularizer
parameters.u_eps = 1e-9;                %velocity regularizer
parameters.u_in = 0./parameters.year; 

parameters.accumrate = accum;
parameters.buttress = 0;

%%Sliding Law Parameters
parameters.frictionlaw = 'Weertman';

parameters.C_schoof = 7.624e6;      %See Schoof (2007)
parameters.m_schoof = 1/3;          %See Schoof (2007)

parameters.B_shear = 0;
parameters.width_shear = 1e3;

parameters.float = 1;
n_plots = 1e2;                        %number of plots to make

%% Pre-allocate storage
h_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
xg_all = nan*ones(1,round(parameters.nsteps/5));
time_all = nan*ones(1,round(parameters.nsteps/5));
% ub_all = nan*ones(parameters.grid.n_nodes,round(parameters.nsteps/5));

%% Initialize variables
disp('Initializing Variables...')
load ConvergedWeertman_SchoofBed_ac0.3_Ag1e-25_N2000_400.mat   %load converged solution from earlier run

% Start use these lines if you want to make a new initial condition
% x_g = 1500e3;
% h_g = -(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
% h = 3500 - (3500-h_g).*parameters.grid.sigma_element.^3;

parameters.uverbose = 1;
[u,error,uiter] = velocity_solve_fix(100.*parameters.grid.sigma_node./parameters.year,h,x_g,parameters); %initialize velocity


%% Run Integration
warning('off')

%%Pre-allocate storage
h_all = nan*ones(parameters.grid.n_elements,round(parameters.nsteps/5));
xg_all = nan*ones(1,round(parameters.nsteps/5));
time_all = nan*ones(1,round(parameters.nsteps/5));
ub_all = nan*ones(parameters.grid.n_nodes,round(parameters.nsteps/5));


%%Run integration
parameters.uverbose = 0;
time = 0;
t=0;
q=0;

while(parameters.year*time <= parameters.tfinal)
    
    t=t+1;
    
    %calculate velocity, thickness and GL position 
    %using a time step that is reduced if either of
    %the solvers reaches max iterations
    parameters.hx_g_old = [h;x_g];
    u_old = u;
    
    hiter=parameters.hiter_max;
    parameters.dtau = parameters.dtau_max; %start with max time step length
    courant = parameters.CFL;
    
    while(hiter==parameters.hiter_max || courant >= parameters.CFL)
        
        %calculate thickness and GL position (with velocity from previous
        %time step)
        [h,x_g,hiter] = thickness_wGL_solve_fix(parameters.hx_g_old(1:end-1),parameters.hx_g_old(end),u_old,parameters);
        
        %calculate velocity (with current thickness and GL pos)
        [u,error,uiter] = velocity_solve_fix(u_old,h,x_g,parameters);
      
        %Calculate CFL number
        u_eff = 0.5.*(u(1:end-1)+u(2:end)) - parameters.grid.sigma_element.*(x_g-parameters.hx_g_old(end))./parameters.dtau;
        courant = max(abs(u_eff./(x_g.*diff(parameters.grid.sigma_node)./parameters.dtau)));
        
        h_g = -(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
        q_sch=(h_g^((parameters.n_Glen+3+parameters.m_schoof)/(1+parameters.m_schoof)))*(((parameters.B_Glen(1)^(-3)).*(parameters.rho.*parameters.g)^(parameters.n_Glen+1) * (1-(parameters.rho/parameters.rho_w))^parameters.n_Glen / ((4^3)*parameters.C_schoof))^(1/(parameters.m_schoof+1))).*(1-parameters.buttress)^(parameters.n_Glen/(1+parameters.m_schoof));
        q_num=u(end).*h_g;
        100*(q_num-q_sch)./q_sch

        %If necessary, adapt time step length
        if(hiter==parameters.hiter_max || courant >= parameters.CFL)
            parameters.dtau = parameters.dtau/2; 
        end
        
    end
    
    %plot solution on occasional time steps
    if((mod(t,parameters.nsteps/n_plots)==1 || t==1))
        if(plotting);FlowlinePlot;end
    end
    
    %save solutions
    q=q+1;
    h_all(:,q) = h;
    xg_all(q) = x_g;
    time_all(q) = time;     

    
    time = time + (parameters.dtau./parameters.year);
    if(sum(isnan(u))>0);disp('Integration FAIL');break;end
    if(parameters.dtau <= (2^(-20)).*parameters.dtau_max);disp('Integration FAIL');break;end
    
end
end