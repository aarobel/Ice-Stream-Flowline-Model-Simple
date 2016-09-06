h_g=-(parameters.rho_w/parameters.rho).*Base(x_g,parameters);
bed_elev = Base(x_g.*parameters.grid.sigma_node,parameters);

x_element=repmat(x_g*parameters.grid.sigma_element,[1 parameters.grid.n2_nodes]);

if(t~=1)
    set(p1,'color','k')
    set(p2,'color','k')

end
    
%%

figure(1);set(1,'units','normalized','position',[0 0.1 0.5 0.75]);

subplot(3,1,1)
area(linspace(0,1500,1000),Base(linspace(0,1500e3,1000),parameters),-1000,'FaceColor',[0.7,0.5,0.2]);hold on;
% plot(x_g.*parameters.grid.sigma_node/1000,bed_elev,'k','linewidth',2);hold on
p2 = plot([x_g.*parameters.grid.sigma_node;x_g]/1000,[bed_elev+[h;h_g];bed_elev(end)],'b','linewidth',3);hold on
plot(x_g/1000,bed_elev(end),'b.','markersize',8)
ylim([-1000 6000])
xlabel('x (km)','fontsize',26);
ylabel('Ice Height (m)','fontsize',26)
set(gca,'fontsize',26)
title(['Time is ' int2str(time) ' years'],'fontsize',26)

subplot(3,1,2)
p1 = plot(parameters.grid.sigma_node,u*parameters.year,'r','linewidth',3);hold on
xlabel('Position with ice sheet from ice divide to GL','Interpreter','LaTeX','fontsize',26);
ylim([min(u*parameters.year) max([u*parameters.year;10])]);
ylabel('Velocity (m/yr)','fontsize',26)
set(gca,'fontsize',26)

subplot(3,1,3)
plot(time./1000,x_g./1e3,'k.','markersize',30);hold on
xlabel('time (kyr)','Interpreter','LaTeX','fontsize',26);
ylabel('GL Position (km)','fontsize',26)
set(gca,'fontsize',26)
drawnow