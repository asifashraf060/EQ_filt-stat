%script to filter out EQ data from USGS catalog based on an interface
%asif Feb26, 2023
clear all, close all

%% INPUT

% path to the USGS CSV file
the_EQs    = '/Users/asifashraf/OneDrive - University Of Oregon/UNL_research/2nd Manuscript_2D-models/EQ_filter/M2_Y1900.csv';
% path to the matlab structure of the interface
the_INT    = '/Users/asifashraf/OneDrive - University Of Oregon/UNL_research/2nd Manuscript_2D-models/EQ_filter/MC_Slab.mat';
% filter value for depth error
dpEr_fltr  = 2;
% filter value for location error
locEr_fltr = 2;
% filter value for magnitude
mg_fltr    = 2;
% path to the save folder
out_dir    = '/Users/asifashraf/OneDrive - University Of Oregon/UNL_research/2nd Manuscript_2D-models/EQ_filter/';
% Boundary of latitude values; start from lowest latitude
lat_lim = [43 44; 44 46; 46 47; 47 48; 48 49];
% number of cluster within those boundary in those exact orders
clstr_n = [1 2 2 3 2];

%% CALCULATION
EQ_tb    = readtable(the_EQs);
EQ_lt    = table2array(EQ_tb(:,2));
EQ_ln    = table2array(EQ_tb(:,3));
EQ_dp    = table2array(EQ_tb(:,4));
EQ_mg    = table2array(EQ_tb(:,5));
EQ_locEr = table2array(EQ_tb(:,16));
EQ_dpEr  = table2array(EQ_tb(:,17));


%% FILTER EARTHQUAKES BASED ON ERROR VALUES
disp('Filering for input parameters ...')
%find indices that contain nans for location and depth error
nan_inds  = unique(vertcat((find(isnan(EQ_locEr) == 1)), (find(isnan(EQ_dpEr) == 1))));
%find indices that contain higher values than filter values
fltr_inds = unique(vertcat( (find(EQ_locEr > locEr_fltr)),...
                                (find(EQ_dpEr > dpEr_fltr)),...
                                    (find(EQ_mg < mg_fltr))));
%make an array of indices by combining fltr_inds and nan_inds
rmv_inds  = unique(vertcat(nan_inds, fltr_inds));
%indices that are alive after the filtering
add_inds  = setdiff([1:length(EQ_lt)], rmv_inds);
%filter out all the arrays of input earthquakes
EQ_lt_f   = EQ_lt(add_inds); EQ_ln_f   = EQ_ln(add_inds);
EQ_dp_f   = EQ_dp(add_inds); EQ_mg_f   = EQ_mg(add_inds);


%% FILTER EARTHQUAKES BASED ON INPUT INTERFACE

disp('Filering for input interface ...')
INT    = load(the_INT);

%use places that doesn't have any nans
ind_ln  = find(isnan(INT.MC_Lon) == 0);
ind_lt  = find(isnan(INT.MC_Lat) == 0);
ind_dp  = find(isnan(INT.MC_depth) == 0);
ind_use = intersect(ind_ln, ind_lt);
ind_use = intersect(ind_use, ind_dp);

INT_ln = INT.MC_Lon(ind_use); %check the fieldname before input
INT_lt = INT.MC_Lat(ind_use);
INT_dp = INT.MC_depth(ind_use);
%perform scattered Interpolation on both Interface and Earthquake depths
disp('Performing scattered interpolation for the inported interface ...')
F = scatteredInterpolant(INT_ln, INT_lt, INT_dp);
[INT_ln_grd, INT_lt_grd] = meshgrid(INT_ln, INT_lt);
INT_dp_grd = F(INT_ln_grd, INT_lt_grd);

%loop through every filtered EQ
clear EQ_LAT; clear EQ_LON; clear EQ_DEPTH; clear EQ_MAG
for i = 1:length(EQ_lt_f)
    lt     = EQ_lt_f(i);
    ln     = EQ_ln_f(i);
    dp     = EQ_dp_f(i);
    mg     = EQ_mg_f(i);
    lt_ind_INT = find(abs(lt - INT_lt) == min(abs(lt - INT_lt)));
    ln_ind_INT = find(abs(ln - INT_ln) == min(abs(ln - INT_ln)));
    int_dp     = abs(INT_dp_grd(lt_ind_INT,ln_ind_INT));
    if dp<int_dp
        EQ_LAT(:,i)     = nan;
        EQ_LON(:,i)     = nan;
        EQ_DEPTH(:,i)   = nan;
        EQ_MAG(:,i)     = nan;
    end
    if dp>int_dp
        EQ_LAT(:,i)     = lt;
        EQ_LON(:,i)     = ln;
        EQ_DEPTH(:,i)   = dp;
        EQ_MAG(:,i)     = mg;
    end
end
disp('Earthquakes are filtered based on the Input interface ...')

% plot the filtered EQ in 3-D
figure(1), clf
plot3(INT_ln, INT_lt, (INT_dp)*(-1), '.')
hold on
plot3(EQ_LON, EQ_LAT, EQ_DEPTH, 'or')
grid on
xlim([-127 -122])
ylim([42 52])
xlabel('Longitude (degree decimal)')
ylabel('Latitude (degree decimal)')
zlabel('Depth (km)')
title('red=EQ||blue=interface')
set(gca, 'ZDir', 'reverse', 'FontSize', 16)

% make all the nan values zero and save it in a text file
EQ_LAT   = EQ_LAT(find(isnan(EQ_LAT)==0));
EQ_LON   = EQ_LON(find(isnan(EQ_LON)==0));
EQ_DEPTH = EQ_DEPTH(find(isnan(EQ_DEPTH)==0));
EQ_MAG   = EQ_MAG(find(isnan(EQ_MAG)==0));

fid = fopen(append(out_dir, 'EQ_filtered_M', string(mg_fltr), '_dpEr', string(dpEr_fltr), '_locEr', string(locEr_fltr)  ,'.txt'), 'w');
for i = 1:length(EQ_LAT)
    fprintf(fid, '%f %f %f %f\n', EQ_LAT(i), EQ_LON(i), EQ_DEPTH(i), EQ_MAG(i));
end
fclose(fid)
disp('Filtered earthquakes are written in the text file')

%% Perform statistical analysis over the EQ locations
disp('Performing k-means clustering...')

% plot the statistical analysis
figure(9913), clf

%plot the arrow of subduction direction
hold on

theta = 68;     %angle of the arrow, N68E is the subduction direction from DeMats et al. (2010)
r     = .5;      %magnitude of the arrow
x     = -125;   %x position to plot
y     = 45;     %y poisiton to plot
u     = r * cosd(90-theta); % convert polar (theta, r) to cartesian
v     = r * sind(90-theta);

h = quiver(x, y, u, v, 'k');

% loop through each latitude boundaries and apply K-means clustering
for ii = 1:length(lat_lim)
    
    % specifying the limits
    lim = lat_lim(ii,:);
    
    % find indices for earthquakes within that boundary
    ind_latLim = find(EQ_LAT>=min(lim) & EQ_LAT<=max(lim));

    % assign lat and lon of those eqs to data variable
    data = [];
    data(:,1) = EQ_LAT(ind_latLim);
    data(:,2) = EQ_LON(ind_latLim);
    
    % perform K-means clustering
    [idx, C] = kmeans(data, clstr_n(ii));

    % loop to go through each cluster and...
    % plot those clustered eqs and fit a best fit line through those eqs
    for i = 1:clstr_n(ii)

        xx = data(:,2); yy = data(:,1);

        % find the indices for values within a single group
        ind_gr = find(idx == i);

        % assign those clustered eqs into new varibles x and y
        x = xx(ind_gr); y = yy(ind_gr); x = reshape(x, [], 1);

        % Fit a linear regression model
        mdl = fitlm(x, y);
        yfit = predict(mdl, x);

        % plot the best fit line
        figure(9913)
        hold on
        plot(x, y, '.', 'MarkerSize', 15)
        hold on
        plot(x, yfit, '-b', 'LineWidth', 2)
        
        % save these clustered positions in a text file
        fid = fopen(append(out_dir, 'lat_bndr', string(ii) ,'_eq_cluster_', string(i),'.txt'), 'w');
        
        for j = 1:length(x)
            fprintf(fid, '%f %f %f\n', x(j), y(j), yfit(j));
        end
        fclose(fid)
    end
end

