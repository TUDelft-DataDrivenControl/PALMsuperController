clear;clc;close all;

wf = '6turb';

addpath(genpath('mexcdf'));
addpath(strcat(wf));
filename = strcat(wf,'_m01.nc');

Animate    = 10;

time       = double(nc_varget(filename,'time'));
u          = double(nc_varget(filename,'u'));
v          = double(nc_varget(filename,'v'));
%w          = double(nc_varget(filename,'w'));
x          = double(nc_varget(filename,'x'));
y          = double(nc_varget(filename,'y'));
xu         = double(nc_varget(filename,'xu'));
yv         = double(nc_varget(filename,'yv'));
zw_3d      = double(nc_varget(filename,'zw_3d'));
zu_3d      = double(nc_varget(filename,'zu_3d'));

dx         = diff(x);
dy         = diff(y);
dxu        = diff(xu);
dyv        = diff(yv);

switch lower(filename)
                      
    case lower('6turb_m01.nc')
        Cry = [1091.0, 1469.0, 1091.0, 1469.0, 1091.0, 1469.0]-yv(1);
        Crx = [6000.0, 6000.0, 6630.0, 6630.0, 7260.0, 7260.0]-xu(1);
        
        %M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power]
        M1  = dlmread('Data_PALM\6turb_adm\6turb_adm_turbine_parameters01.txt','',1,0);
        M2  = dlmread('Data_PALM\6turb_adm\6turb_adm_turbine_parameters02.txt','',1,0);
        M3  = dlmread('Data_PALM\6turb_adm\6turb_adm_turbine_parameters03.txt','',1,0);
        M4  = dlmread('Data_PALM\6turb_adm\6turb_adm_turbine_parameters04.txt','',1,0);
        M5  = dlmread('Data_PALM\6turb_adm\6turb_adm_turbine_parameters05.txt','',1,0);
        M6  = dlmread('Data_PALM\6turb_adm\6turb_adm_turbine_parameters06.txt','',1,0);
        
        t     = M1(:,1);
        Ur    = [M1(:,2)';M2(:,2)';M3(:,2)';M4(:,2)';M5(:,2)';M6(:,2)'];
        Uinf  = [M1(:,3)';M2(:,3)';M3(:,3)';M4(:,3)';M5(:,3)';M6(:,3)'];
        phi   = [M1(:,6)';M2(:,6)';M3(:,6)';M4(:,6)';M5(:,6)';M6(:,6)'];
        P     = [M1(:,8)';M2(:,8)';M3(:,8)';M4(:,8)';M5(:,8)';M6(:,8)'];
        CT    = [M1(:,4)';M2(:,4)';M3(:,4)';M4(:,4)';M5(:,4)';M6(:,4)'];
        F     = [M1(:,7)';M2(:,7)';M3(:,7)';M4(:,7)';M5(:,7)';M6(:,7)'];
        a     = [M1(:,5)';M2(:,5)';M3(:,5)';M4(:,5)';M5(:,5)';M6(:,5)'];
        nz    = 1;
        Dr    = 120;
        u_Inf = u(1,nz,1,1);        
end

N     = size(Crx,2);
x     = x-x(1);
y     = y-y(1);
xu    = xu-xu(1);
yv    = yv-yv(1);

for i = 1:length(Crx)
    [~,xline(i,1)]     = min(abs(x-Crx(i)));
    [ML_prim, L_prim ] = min(abs(y- (Cry(i)-Dr/2)));
    [MR_prim, R_prim ] = min(abs(y- (Cry(i)+Dr/2)));
    yline{i}           = L_prim:1: R_prim;
end


%% Plot flow fields
scrsz = get(0,'ScreenSize');
if Animate>0
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

st = 1;
turb_coords = (.5*Dr*exp(1i*phi(:,st)*pi/180)).';  % Yaw angles
uk          = squeeze(u(st,nz,:,:)).';
vk          = squeeze(v(st,nz,:,:)).';

if Animate > 0
    subplot(2,3,1)
    [~,p1] = contourf(y,xu,uk,'Linecolor','none');  colormap(hot);
    caxis([min(min(min(u(:,nz,:,:)))) max(max(max(u(:,nz,:,:))))+.01]);
    hold all; cb1 = colorbar; axis equal; axis tight;
    title('$u$ [m/s]','interpreter','latex');
    t1 = text(0,x(end)+200,['Time ', num2str(time(st),'%.1f'), 's'],'interpreter','latex');
    uTurbs = plot([Cry-real(turb_coords); Cry+real(turb_coords)],...
        [Crx-imag(turb_coords); Crx+imag(turb_coords)],'k','linewidth',1);
    hold off
    
    subplot(2,3,2)
    [~,p2] = contourf(yv,x,vk,'Linecolor','none');  colormap(hot);
    caxis([min(min(min(v(:,nz,:,:)))) max(max(max(v(:,nz,:,:))))+.01]);
    hold all; cb2 = colorbar; axis equal; axis tight;
    title('$v$ [m/s]','interpreter','latex');
    vTurbs = plot([Cry-real(turb_coords); Cry+real(turb_coords)],...
        [Crx-imag(turb_coords); Crx+imag(turb_coords)],'k','linewidth',1);
    hold off
    
    subplot(2,3,3)
    p3 = plot(t(st),P(:,st));
    title('$P$ [W]','interpreter','latex');
    axis([0,t(end) 0 max(max(P(:,st:end)))+10^5]); grid;
    
    subplot(2,3,4)
    p4 = plot(y,uk(floor((Crx(1)+Dr)/dx(1)),:),'r' );
    ylabel('$u$','interpreter','latex');xlabel('$y$','interpreter','latex');
    title('wake cross-section','interpreter','latex');
    axis tight;axis([0,y(end) 2 12]); grid;
    
    subplot(2,3,5)
    p5 = plot(x,mean(uk(:, yline{1}),2));
    title('$U_c$ [m/s]','interpreter','latex'); grid;
    axis([0,x(end) 0 max(max(mean(uk(:, yline{1}),2)))+1]);
    
    subplot(2,3,6)
    if exist('CT','var')
        p6 = plot(t(st),CT(:,st));
        ylabel('$C_T$','interpreter','latex')
        title('thrust coeficient []','interpreter','latex');
        axis([0,t(end) 0 max(max(CT(:,st:end)))+.2]);
    else
        p6 = plot(t(st),Ft(:,st));
        ylabel('$F_t$','interpreter','latex')
        title('thrust force [N]','interpreter','latex');
        axis([0,t(end) 0 max(max(Ft(:,st:end)))+10^5]);
    end
    grid
    drawnow
    
    for k=st:Animate:min(length(t),length(time))
        turb_coords  = (.5*Dr*exp(1i*phi(:,k)*pi/180)).';  % Yaw angles
        uk          = squeeze(u(k,nz,:,:)).';
        vk          = squeeze(v(k,nz,:,:)).';
        
        set(p1,'Zdata',uk);
        set(p2,'Zdata',vk);
        
        set(cb1,'Limits',[min(min(uk)) max(max(uk))+.01]);
        set(cb2,'Limits',[min(min(vk)) max(max(vk))]);
        
        set(t1,'String',['Time ', num2str(time(k),'%.1f'), 's']);
        for ll=1:N
            set(uTurbs(ll),'XData',[Cry(ll)-real(turb_coords(ll)) Cry(ll)+real(turb_coords(ll))])
            set(uTurbs(ll),'YData',[Crx(ll)-imag(turb_coords(ll)) Crx(ll)+imag(turb_coords(ll))])
            
            set(vTurbs(ll),'XData',[Cry(ll)-real(turb_coords(ll)) Cry(ll)+real(turb_coords(ll))])
            set(vTurbs(ll),'YData',[Crx(ll)-imag(turb_coords(ll)) Crx(ll)+imag(turb_coords(ll))])
            
            set(p3(ll),'XData',t(st:k))
            
            
            set(p3(ll),'YData',P(ll,st:k))
            
            if exist('ywakef','var')
                set(p6(1),'XData',t(st:k))
                set(p6(1),'YData',ywakef(1,st:k));
                set(p7(1),'XData',t(st:k))
                set(p7(1),'YData',ywake_r(1,st:k));
            elseif exist('CT','var')
                set(p6(ll),'XData',t(st:k))
                set(p6(ll),'YData',CT(ll,st:k))
            else
                set(p6(ll),'XData',t(st:k))
                set(p6(ll),'YData',Ft(ll,st:k))
            end
        end
        
        set(p5,'YData',mean(uk(:, yline{1}),2))
        set(p4,'YData',uk(floor((Crx(1)+Dr)/dx(1)),:))
        drawnow
        
    end
end

   