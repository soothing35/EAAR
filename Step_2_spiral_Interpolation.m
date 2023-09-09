function interpolation_file = Spiral_Interpolation(picking_file,pix_size,particle_size,round)
 % picking_file is the name of file saving picking result of a spiral got from spiral_picker_gui
 % such as 0206_1
 % round is such as 1, 2, 3 etc
    if nargin<2
        pix_size = 0.106;  %0.065;  %1.3 %nm/pixel
        particle_size = 2.9;  %3.2;%nm
        round = 1;
    elseif nargin<3
        particle_size = 2.9;  %3.2;%nm
        round = 1;
    elseif nargin<4
        round = 1;
    end
     
    m_pth = fileparts(which('Spiral_Interpolation.m'));   
    [parentdir,~,~] = fileparts(m_pth);
    picking_file_path =[parentdir '/Picking-result/'];   %picking file path
    
    if round==1
        a = load([picking_file_path picking_file '.mat']);
    else
        pth = which([picking_file '.mat']); 
        a = load(pth);
    end

    %% determine the step size in spacing data
    step_size = ceil(particle_size/pix_size); %pixel
   
    
    %% Beizer Interpolation
    k =  a.Points_position;
    [cp1,cp2] = findControlPoints(k); 
    [Inter_points] = displayBezier(k,cp1,cp2);
    
    
    %% calculate the arclength of the curve
    px = Inter_points(:,1);
    py = Inter_points(:,2);
    [arclen,seglen] = arclength(px,py);

    %% calculate the curvature of the spicific points 
    num_InterPoints = size(Inter_points,1);
    Curv = zeros(num_InterPoints,1);
    %tic
    for j = 1:num_InterPoints-2
        x1 = Inter_points(j,1);
        y1 = Inter_points(j,2);
    
        x2 = Inter_points(j+1,1);
        y2 = Inter_points(j+1,2);
    
        x3 = Inter_points(j+2,1);
        y3 = Inter_points(j+2,2);
   
        Curv(j+1) = 2*((x2-x1)*(y3-y2)-(y2-y1)*(x3-x2)) / sqrt( ...
    ((x2-x1)^2+(y2-y1)^2)*((x3-x2)^2+(y3-y2)^2)*((x1-x3)^2+(y1-y3)^2));
    end
    %toc
 
    %% calculate the slope( tangent line)
    %approximate the derivative (slope) at each point with the following second order differencing
    t = Inter_points(:,1)'; p = Inter_points(:,2)';
    n = size(Inter_points,1);
    td = [t(3),t(1:n-1)]; tu = [t(2:n),t(n-2)];
    pd = [p(3),p(1:n-1)]; pu = [p(2:n),p(n-2)];
    dpdt = ((pu-p)./(tu-t).*(t-td)+(p-pd)./(t-td).*(tu-t))./(tu-td);
    slope = dpdt';
    
    %% calculate the spacing data
    point_num =  ceil(arclen/step_size);
    [linspace_x,linspace_y] = linspacearc(Inter_points(:,1),Inter_points(:,2),point_num);
    linspace_points = [linspace_x' linspace_y'];
    
    
    %% calculate the slope for spacing data( tangent line)
    %approximate the derivative (slope) at each point with the following second order differencing
    t = linspace_points(:,1)'; p = linspace_points(:,2)';
    n = size(linspace_points,1);
    td = [t(3),t(1:n-1)]; tu = [t(2:n),t(n-2)];
    pd = [p(3),p(1:n-1)]; pu = [p(2:n),p(n-2)];
    dpdt = ((pu-p)./(tu-t).*(t-td)+(p-pd)./(t-td).*(tu-t))./(tu-td);
    linspace_slope = dpdt';
    
    
    if round ==1
        cd (m_pth);
        cd ../; 
        pp=pwd;
        if ~exist('Results','dir')    %% Results folder is in the upper level of .m file
            mkdir('Results');
        end
        %cd([pp,'\Results']);
         cd([pp,'/Results']);
        if ~exist(picking_file,'dir') 
            mkdir(picking_file);
        end
        cd(picking_file);         
        save([picking_file '_Interpolation_' num2str(round) '.mat'],'Inter_points','arclen','seglen','Curv','slope','linspace_points','step_size','linspace_slope');
        interpolation_file = ([picking_file '_Interpolation_' num2str(round)]);
        path (path,[pp '/' 'Results']);
        path (path, [pp '/' 'Results' '/' picking_file]);
    else
        pos = strfind(picking_file,'_');
        picking_file = picking_file(1:pos(size(pos,2)-1)-1);
        
        cd (m_pth);
        cd ../; pp=pwd;
        if ~exist('Results','dir') 
            mkdir('Results');
        end
        cd('Results');

        if ~exist(picking_file,'dir') 
            mkdir(picking_file);
        end
        cd(picking_file);         
        save([picking_file '_Interpolation_' num2str(round) '.mat'],'Inter_points','arclen','seglen','Curv','slope','linspace_points','step_size','linspace_slope');
        interpolation_file = ([picking_file '_Interpolation_' num2str(round)]);    
    end
    fclose('all'); 
    
    cd (m_pth);

end
    
