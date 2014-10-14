%% Code by Jannis Uhlendorf 2012

classdef HoughTracker < handle
  
  properties
    cells_radii    = [6, 20]; % pixels
    old_circen     = [];
    old_circen2    = [];

    mapping        = [];
    %id_counter     = 0;
    data_int       = [];
    data_coloc     = [];
    data_size      = [];
    
    movie_tracking = [];
    movie_nuc      = [];
    tp             =  0;
    base_name      = 0;
    lost_cells     = [];
    lost_cells_mapped = [];
  end

  methods

    function obj = HoughTracker( base_name )
      obj.base_name = base_name;
      %obj.movie_tracking = avifile( [ base_name, '_tracking.avi'] );
      %obj.movie_nuc      = avifile( [ base_name, '_nucleus.avi'] );
    end

    function delete( obj )
      %obj.movie_tracking = close( obj.movie_tracking );
      %obj.movie_nuc      = close( obj.movie_nuc );
    end
    
    function x = add_timepoint( obj, img_trans, img_expr, img_hog, img_nuclear )
      [accum, circen, cirrad] = obj.hough_transform( img_trans );
      obj.tp = obj.tp+1;
      if ~isempty(obj.old_circen)
	% tracking
	% append cells that have been lost in the round before
	old = [ obj.old_circen; obj.old_circen2(obj.lost_cells,:) ]; 
	m=obj.tracking( old, circen, cirrad );
	%m = zeros( length(old), length(cirrad) );

	% separate result in normal and refound cells
	m_new = m(1:end-length(obj.lost_cells),:);
	m_old = m(end-length(obj.lost_cells)+1:end,:);
	
	% determine which cells have been lost (lost for original order, lost mapped for mapped order)
	lost = find( all(m_new==0,2) );
	[lost_mapped,dummy] = find(obj.mapping(:,lost));
	% old code, the one above does the same
	%v = zeros( size(obj.mapping,2), 1 );
	%v(lost)=1:length(lost);
	%l_ind = obj.mapping*v;
	%lost_mapped=zeros(length(lost),1);
	%for i=1:length(lost)
	%  lost_mapped(i) = find(l_ind==i);
	%end

	obj.mapping = obj.mapping * m_new; % update permutation matrix

	% modify mapping to bring in lost cells again
	if ~isempty(obj.lost_cells)
	  for o=find(any(m_old,2))'
	    obj.mapping( obj.lost_cells_mapped(o), : ) = m_old( o, : );
	  end
        end
	for new=find(any(obj.mapping,1)==0) % append row for new objects
	  obj.mapping = [obj.mapping; 1:size(obj.mapping,2)==new];
	end
	obj.lost_cells = lost;
	obj.lost_cells_mapped = lost_mapped;
	obj.old_circen2 = obj.old_circen;
      else
	obj.mapping = eye(size(circen,1)); % first permutation matrix is identity
      end
      obj.old_circen = circen;

      ids = (obj.mapping')*(1:size(obj.mapping,1))';
      ids(ids==0)=[];
      obj.plot_cells( img_trans, circen, cirrad, ids );
      
      [size_x,size_y]=size(img_trans);
      [X,Y] = meshgrid(1:size_y, 1:size_x);

      cirrad = cirrad;
      mask = zeros(size_x,size_y);
      for ii = 1 : size(circen, 1)
	this_mask = (X-circen(ii,1)).^2 + (Y-circen(ii,2)).^2 < cirrad(ii)^2;	  
	mask = mask + ii*this_mask;
	mask = min( mask, ii ); % if cells are overlapping ...
      end

      cellsize = obj.mapping * cirrad;
      cellsize = [obj.tp; cellsize];
      if ~isempty(obj.data_size) & size(obj.data_size,1)<size(cellsize,1)
	num_new = size(cellsize,1)-size(obj.data_size,1);
	obj.data_size = [obj.data_size; NaN*ones(num_new,size(obj.data_size,2))];
	obj.data_int  = [obj.data_int; NaN*ones(num_new,size(obj.data_int,2))];
	obj.data_coloc  = [obj.data_coloc; NaN*ones(num_new,size(obj.data_coloc,2))];
      end
      obj.data_size = [obj.data_size cellsize];

      if img_expr~=0
	expr = obj.get_intensities( img_trans, img_expr, 0, mask, cirrad );
	expr(isnan(expr))=0;
	expr_new = obj.mapping * expr';
	expr_new(expr_new==0)=NaN;
	x = expr_new;
	expr_new = [obj.tp; expr_new];
	obj.data_int  = [obj.data_int expr_new];
      end
      if (img_hog~=0) & (img_nuclear~=0)
	coloc = obj.colocalization( img_trans, img_hog, img_nuclear, mask, cirrad );
	coloc(isnan(coloc))=0;
	coloc_new = obj.mapping * coloc';
	coloc_new(coloc_new==0)=NaN;
	coloc_new = [obj.tp; coloc_new];
	obj.data_coloc = [obj.data_coloc coloc_new];
      end
      %x = expr_new;
    end

    function x = get_intensities( obj, img_trans, img_gfp, img_rfp, mask, cirrad )
      x=zeros(1,max(mask(:)));
      for i=1:max(mask(:))
	this_mask = (mask==i);
	crop = img_gfp( this_mask );
	x(i) = sum(crop(:)) / numel(crop);
	%x(i) = sum(crop(:)) / cirrad(i)^2;
      end
    end

    function x = colocalization( obj, img_trans, img_gfp, img_rfp, mask, cirrad )
      x=zeros(1,max(mask(:)));
      t         = graythresh( img_rfp );
      t_img_rfp = im2bw( img_rfp, t*1.5 );
      
      nuc_gfp   = img_gfp.*uint16(t_img_rfp);
      cyt_gfp   = img_gfp - nuc_gfp;

      % debug .. write nuclear images to disk
      h = figure('visible', 'off');
      imshow( cyt_gfp, [min(cyt_gfp(:)), max(cyt_gfp(:))] );
      saveas( h, [ obj.base_name sprintf( '_tracker/img_nuclear_%04d.tif', obj.tp) ], 'tif' );

      for i=1:max(mask(:))
	this_mask = (mask==i);
	crop_gfp  = img_gfp( this_mask );
	crop_nuc_gfp = nuc_gfp( this_mask );
	crop_cyt_gfp = cyt_gfp( this_mask );

	x_in  = sum(crop_nuc_gfp)/sum( crop_nuc_gfp~=0 );
	x_out = sum(crop_cyt_gfp)/sum( crop_cyt_gfp~=0 );
	x(i) = x_in/x_out;
      end
    end

    function [accum, circen, cirrad] = hough_transform( obj, img_trans )
      %[accum, circen, cirrad] = CircularHough_Grd(double(img_trans), obj.cells_radii, 20, 15, 1)
      %[accum, circen, cirrad] = CircularHough_Grd(double(img_trans), obj.cells_radii, 20, 15, 1);
      [accum, circen, cirrad] = CircularHough_Grd(double(img_trans), obj.cells_radii, 50, 14, 1);
      if any(cirrad <= 0)
        inds = find(cirrad>0);
        cirrad = cirrad(inds);
        circen = circen(inds,:);
      end
					 
      % remove multiple detections
      M_x = circen(:,1)*ones(1,size(circen,1));
      M_y = circen(:,2)*ones(1,size(circen,1));
      identical = (abs(M_x - M_x') + abs(M_y - M_y'))<.1;
      identical = triu(identical, 1 );
      [x,y] = find( identical );
      circen(y,:)= [];
      cirrad(y)  = [];

    end

    function x = tracking( obj, old_centers, new_centers, cirrad )

      function dm = distance_matrix(a,b)
        % custom function to compute distance matrix between two 1d-coordinate vectors
        am = a*ones(1,length(b));
        bm = b*ones(1,length(a));
        dm=(am-bm').^2;
        %dm = log(abs(am-bm'));
      end
      
      D = distance_matrix(old_centers(:,1),new_centers(:,1)) + distance_matrix(old_centers(:,2),new_centers(:,2));
      % make the distance between two cells which lie in the same radius smaller than zero
      P = [];
      for i=1:size(old_centers,1)
        P = [ P; sqrt(D(i,:))<cirrad' ];
      end
      M = P.*D;
      m = max( M(:) );
      D = D - m - .001;

      %x=obj.linear_constraint_mapping( D );
      x=obj.linear_constraint_mapping_lp_solve( D );
    end

    function mapping = linear_constraint_mapping( obj, D )
      % solve tracking from distance matrix D using binary integer programming
      
      % define the problem as c'*b != max while Ax<b and A_eq*x=b_eq
      c = reshape(D',1,[]);
      [dim_1, dim_2] = size(D);
      dd = dim_1*dim_2;
      
      % construct contraints matrix
      A=[];
      b=ones(dim_1+dim_2,1);
      for i=1:dim_1     % constraints preventing multiple mappings for one row-object
        a = zeros(1,dd);
        a(((i-1)*dim_2)+1:i*dim_2)=1;
        A=[A; a];
      end
      for i=1:dim_2     % contraints prevening multiple mappings for one column-object
        x = zeros(1,dim_2);
        x(i)=1;
        a = repmat(x, 1, dim_1 );
        A = [A; a];
      end

      % make sure all mappeable objects are mapped
      A_eq = ones(1,dd);
      b_eq = min(dim_1, dim_2);

      x_0 = reshape( eye(dim_1,dim_2)', 1, [])';

      % do the actual mapping
      options = optimset( 'TolXInteger', 0 );
      %opt = bintprog( c, A, b, A_eq, b_eq, x_0 ); %, options);
      opt = bintprog( c, A, b ); %, A_eq, b_eq, x_0 ); %, options); 
      m=reshape(opt',dim_2,dim_1)';
      [a,b]=find(m);
      mapping=zeros(1,max(b));
      mapping(b)=a;
      mapping = m;
    end

    function mapping = linear_constraint_mapping_lp_solve( obj, D )
        c = reshape(D',1,[]);
        [dim_1, dim_2] = size(D);
        dd = dim_1*dim_2;
        
        lp=mxlpsolve('make_lp', 0, length(c));
        mxlpsolve('set_verbose', lp, 3);
        mxlpsolve('set_obj_fn', lp, c);
        %mxlpsolve('add_constraint', lp, [0, 78.26, 0, 2.9], 2, 92.3);
        
        for i=1:dim_1
            a = zeros(1,dd);
            a(((i-1)*dim_2)+1:i*dim_2)=1;
            mxlpsolve('add_constraint', lp, a, 1, 1 );
        end
        for i=1:dim_2
            x = zeros(1,dim_2);
            x(i)=1;
            a = repmat(x, 1, dim_1 );
            mxlpsolve('add_constraint', lp, a, 1, 1 );
        end
        %mxlpsolve('add_constraint', lp, ones(1,dd), 3, min(dim_1,dim_2) );
        for i=1:length(c)
            %mxlpsovle('set_binary', lp, i, 'TRUE');
            mxlpsolve('set_binary', lp, i, 1 );
        end
        
        mxlpsolve('solve',lp);
        opt=mxlpsolve('get_variables', lp );
        m=reshape(opt',dim_2,dim_1)';
        [a,b]=find(m);
        mapping=zeros(1,max(b));
        mapping(b)=a;
        mapping = m;
    end


    function pic = plot_cells( obj, img_trans, circen, cirrad, ids )
      h = figure('visible', 'off');
      imshow( img_trans, [min(img_trans(:)), max(img_trans(:))] );
      hold on;
      %plot(circen(:,1),circen(:,2),'r+', 'linewidth', 2)
      for i=1:size(circen,1)
        rectangle('Position',[circen(i,1) - cirrad(i), circen(i,2) - cirrad(i), 2*cirrad(i), 2*cirrad(i)], 'Curvature', [1,1], 'edgecolor', 'b', 'linewidth', 1.5);
	text( circen(i,1)-5, circen(i,2)-5, num2str(ids(i)), 'Color', 'red' );
	%text( circen(i,1)-5, circen(i,2)+25, num2str(floor(circen(i,1))), 'Color', 'red', 'fontsize', 8 );
	%text( circen(i,1)-5, circen(i,2)+35, num2str(floor(circen(i,2))), 'Color', 'red', 'fontsize', 8 );
      end
      hold off;

      saveas( h, [ obj.base_name sprintf( '_tracker/img_tracking_%04d.tif', obj.tp) ], 'tif' );
    end
    




  end
end
