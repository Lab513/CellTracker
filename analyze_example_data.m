% example scrit how to use the HoghTracker class
clear all;
clear classes;


% Example 1: measuring cytosolic fluorescent expression
f_path_trans     = 'example_data/cytosolic_expression/img_trans_%04d.tif';
f_path_expr      = 'example_data/cytosolic_expression/img_yfp_%04d.tif';

% initialize HoughTracker class
% argument for the constructor is a name-string
% output images are saved in the folder "name-string"_tracker
ht = HoughTracker('pics_tracker_expr');

for i=1:3:25
    % rad the trans image
    img_trans = uint16( imread( sprintf( f_path_trans, i ) ) );
    % if a fluorescent image has been acquired, read it
    if exist( sprintf( f_path_expr, i ) )
        img_expr  = uint16( imread( sprintf( f_path_expr, i ) ) ); 
    else
        img_expr = 0;
    end
    % do another time point
    ht.add_timepoint( img_trans, img_expr, 0, 0 );
end

% get the data and plot it
% the first row is the time (int time-points)
time = ht.data_int(1,:)';
% all apart from the first fow is data
expr_data = ht.data_int(2:end,:)';
figure
plot(time, expr_data);









% Eample 2: measuring nuclear localization of a transcription factor
f_path_trans     = 'example_data/nuclear_localization/img_trans_%04d.tif';
f_path_expr      = 'example_data/nuclear_localization/img_yfp_%04d.tif';
f_path_nucleus   = 'example_data/nuclear_localization/img_cfp_%04d.tif';
f_path_tf        = 'example_data/nuclear_localization/img_rfp_%04d.tif';

ht = HoughTracker('pics_tracker_loc', 'cells_radii', [6, 30]);


for i=61:70
    img_trans   = uint16( imread( sprintf( f_path_trans, i ) ) );
    img_nucleus = uint16( imread( sprintf( f_path_nucleus, i ) ) );
    img_tf      = uint16( imread( sprintf( f_path_tf, i ) ) );
    
    % for some time point we have also expression data
    if exist( sprintf( f_path_expr, i ) )
        img_expr = uint16( imread( sprintf( f_path_expr, i ) ) );
    else
        img_expr = 0;
    end
    % do another time point
    ht.add_timepoint( img_trans, img_expr, img_tf, img_nucleus ); 
end

% plot the data
figure
time_int = ht.data_int(1,:);
data_int = ht.data_int(2:end,:)';
plot(time_int, data_int)
title('expression')
figure
time_coloc = ht.data_coloc(1,:);
data_coloc = ht.data_coloc(2:end,:)';
plot(time_coloc, data_coloc)
title('co-localization')

% plot also cell-radius
figure
plot(time_coloc, ht.data_size')
title('cell-size')

