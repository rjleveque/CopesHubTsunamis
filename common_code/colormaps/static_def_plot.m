clear

%colormap for deformation plot
load defcmap.mat

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
   
% Plot deformation for ft
%% Read ft deformation file in most format
fname='../random-mur13-deep_def_15s.most';  % Name of file in MOST format
fid1=fopen(fname,'r'); 
[A]=fscanf(fid1,'%f',[1 2]);
[lon_topo]=fscanf(fid1,'%f',[A(1)]);
[lat_topo]=fscanf(fid1,'%f',[A(2)]);
[def1]=fscanf(fid1,'%f',[A(1) A(2)]);
fclose(fid1);
def1=-def1';
%% plot ft deformation
ax1=axes('Position',[0.05 0.05 0.27 0.9]);
pcolor(lon_topo-360,lat_topo,def1); shading interp;
caxis([-4 4]);
colormap(ax1,defcmap);
colorbar('position',[0.28 0.08 0.015 0.25]);
hold on
contour(lon_topo-360,lat_topo,d1,[0 0],'k','LineWidth',2);
xlabel('Lon');
ylabel('Lat');
set(gca,'FontSize',16,'linewidth',2);
title(['ft: vertical' newline 'displacement (m)']);

% Plot deformation for nosub
%% Read nosub deformation file in most format
fname='../random-mur13-deep-nosub_def_15s.most';  % Name of file in MOST format
fid1=fopen(fname,'r'); 
[A]=fscanf(fid1,'%f',[1 2]);
[lon_topo]=fscanf(fid1,'%f',[A(1)]);
[lat_topo]=fscanf(fid1,'%f',[A(2)]);
[def2]=fscanf(fid1,'%f',[A(1) A(2)]);
fclose(fid1);
def2=-def2';
%% plot nosub deformation
ax2=axes('Position',[0.37 0.05 0.27 0.9]); 
pcolor(lon_topo-360,lat_topo,def2); shading interp;
caxis([-4 4]);
colormap(ax2,defcmap);
colorbar('position',[0.6 0.08 0.015 0.25]);
hold on
contour(lon_topo-360,lat_topo,d1,[0 0],'k','LineWidth',2);
xlabel('Lon');
        %ylabel('Lat');
set(gca,'FontSize',16,'linewidth',2);
title(['nosub: vertical' newline 'displacement (m)']);

% Plot deformation for nosub-nomogeneous
%% Read nosub-nomogeneous deformation file in most format
fname='../random-mur13-deep-nosub-homogeneous-halfspace-noQ_def_15s.most';  % Name of file in MOST format
fid1=fopen(fname,'r'); 
[A]=fscanf(fid1,'%f',[1 2]);
[lon_topo]=fscanf(fid1,'%f',[A(1)]);
[lat_topo]=fscanf(fid1,'%f',[A(2)]);
[def3]=fscanf(fid1,'%f',[A(1) A(2)]);
fclose(fid1);
def3=-def3';

%% plot nosub-homogeneous deformation
ax3=axes('Position',[0.69 0.05 0.27 0.9]); 
pcolor(lon_topo-360,lat_topo,def3); shading interp;
caxis([-4 4]);
colormap(ax3,defcmap);
colorbar('position',[0.92 0.08 0.015 0.25]);
hold on
contour(lon_topo-360,lat_topo,d1,[0 0],'k','LineWidth',2);
xlabel('Lon');
        %ylabel('Lat');
set(gca,'FontSize',16,'linewidth',2);
title(['nosub-homogeneous: vertical' newline 'displacement (m)']);
        