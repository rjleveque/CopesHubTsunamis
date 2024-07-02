clear

%colormap for max tsunami amplitude plot
load maxcmap.mat

% Read bathy/topo data in most format
fname='/Users/yongwei/Documents/work/2020_CoPe_proposal/Audrey_sources_testing/model_grids/csz_gebco2021_15s.most';  % Name of file in MOST format
fid1=fopen(fname,'r'); 
[A]=fscanf(fid1,'%f',[1 2]);
[lon_topo]=fscanf(fid1,'%f',[A(1)]);
[lat_topo]=fscanf(fid1,'%f',[A(2)]);
[d1]=fscanf(fid1,'%f',[A(1) A(2)]);
fclose(fid1);
d1=-d1';
nxC=length(lon_topo);
nyC=length(lat_topo);
set(gcf,'position',[100 100 nxC*cosd(40)*3*0.3 nyC*0.3]);
    
% Plot max tsunami amplitude for ft
%% Read max tsunami amplitude for ft in netcdf format producde by MOST
maxh=ncread('random_mur13_deep_maxh.nc','MAX_HEIGHT');
maxh=maxh';
%% Plot max tsunami amplitude for ft
ax1=axes('Position',[0.05 0.05 0.27 0.9]);
pcolor(lon_topo-360,lat_topo,maxh/100); shading interp;
caxis([-4 20]);
colormap(ax1,maxcmap);
colorbar('position',[0.28 0.08 0.015 0.25]);
hold on
contour(lon_topo-360,lat_topo,d1,[0 0],'k','LineWidth',2);
xlabel('Lon');
ylabel('Lat');
set(gca,'FontSize',16,'linewidth',2);
title(['ft: Max Tsunami' newline 'Surface Elevation (m)']);
clear maxh

%% Read max tsunami amplitude for nosub in netcdf format producde by MOST
maxh=ncread('random_mur13_deep_nosub_maxh.nc','MAX_HEIGHT');
maxh=maxh';
%% Plot max tsunami amplitude for nosub
ax2=axes('Position',[0.37 0.05 0.27 0.9]); 
pcolor(lon_topo-360,lat_topo,maxh/100); shading interp;
caxis([-4 20]);
colormap(ax2,maxcmap);
colorbar('position',[0.6 0.08 0.015 0.25]);
hold on
contour(lon_topo-360,lat_topo,d1,[0 0],'k','LineWidth',2);
xlabel('Lon');
set(gca,'FontSize',16,'linewidth',2);
title(['nosub: Max Tsunami' newline 'Surface Elevation (m)']);
clear maxh

%% Read max tsunami amplitude for nosub-homogeneous in netcdf format producde by MOST
maxh=ncread('random-mur13-deep-nosub-homogeneous-halfspace_maxh.nc','MAX_HEIGHT');
maxh=maxh';
%% Plot max tsunami amplitude for nosub-homogeneous
ax3=axes('Position',[0.69 0.05 0.27 0.9]); 
pcolor(lon_topo-360,lat_topo,maxh/100); shading interp;
caxis([-4 20]);
colormap(ax3,maxcmap);
colorbar('position',[0.92 0.08 0.015 0.25]);
hold on
contour(lon_topo-360,lat_topo,d1,[0 0],'k','LineWidth',2);
xlabel('Lon');
set(gca,'FontSize',16,'linewidth',2);
title(['nosub-homogeneous: Max Tsunami' newline 'Surface Elevation (m)']);
