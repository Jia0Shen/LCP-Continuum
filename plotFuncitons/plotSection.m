function plotSection(g_array,radius,color)

    %Author: Caleb Rucker - May 10 2011
    %This function plots a tubular surface along a curve to represent a section
    %of a concentric-tube robot.
    %Inputs:
    %g_array = an N x 16 matrix where N is the number of points along the
    %curve. Each row  should be a reshaped 4x4 homogeneous transformation matrix
    %(g_array(i,:)=reshape(ith_g_matrix,1,16);
    %radius = desired radius of the plotted tube
    %color = desired color of the tube (eg. color = [1 0 0])
    %After calling this function, it will probably be a good idea to add
    %lighting with the command - camlight headlight, as shown in the example
    %file


    %Get size info
    [m,n]=size(g_array);

    %Set the number of axial rectangles to draw on the tube surface (sort of a
    %"resolution")
    radial_pts = 5; 

    %First draw a filled in circle at the base of the section 
    %*************************************************************************
    g_0=reshape(g_array(1,:),4,4); %First frame

    %make a parameter going from 0 to 2pi for calculating circle points
    tcirc = [0:2*pi/(radial_pts-1):2*pi];

    %Make points on a circle in the local x-y plane
    basecirc = [radius*sin(tcirc);radius*cos(tcirc);zeros(1,length(tcirc));ones(1,length(tcirc))];    

    %Transform these points into the world frame
    basecirc_rot = g_0*basecirc;

    %Pick out x, y, and z components of these points
    xedge = basecirc_rot(1,:);      
    yedge = basecirc_rot(2,:);
    zedge = basecirc_rot(3,:);

    %Draw patches to fill in the circle
    patch(xedge,yedge,zedge,color,'EdgeAlpha',0)
    material shiny
    hold on
    %*************************************************************************

    %Loop to draw each cylindrical segment
    %*************************************************************************
    for i=1:m-1
        g=reshape(g_array(i,:),4,4); %The current frame
        
        g_ahead=reshape(g_array(i+1,:),4,4); %The next frame
        basecirc_rot = g*basecirc; %Current frame circle points
        basecirc_rot_ahead = g_ahead*basecirc; %Next frame circle points     

        
        %Loop to draw each square patch for this segment
        for j=1:radial_pts-1
        xedge = [basecirc_rot(1,j),basecirc_rot(1,j+1),basecirc_rot_ahead(1,j+1),basecirc_rot_ahead(1,j)];      
        yedge = [basecirc_rot(2,j),basecirc_rot(2,j+1),basecirc_rot_ahead(2,j+1),basecirc_rot_ahead(2,j)];
        zedge = [basecirc_rot(3,j),basecirc_rot(3,j+1),basecirc_rot_ahead(3,j+1),basecirc_rot_ahead(3,j)];
        patch(xedge,yedge,zedge,color,'EdgeAlpha',0)
        end

    end
    %**************************************************************************

    %Last, draw a filled in circle at the end of the section
    %*************************************************************************
    g_end=reshape(g_array(end,:),4,4);%Last frame

    %Transform points into the world frame
    basecirc_rot = g_end*basecirc;

    %Get edge points
    xedge = basecirc_rot(1,:);      
    yedge = basecirc_rot(2,:);
    zedge = basecirc_rot(3,:);
    
  
    

    %Draw patches to fill in the circle
    patch(xedge,yedge,zedge,color,'EdgeAlpha',1)
    material shiny
    
    

end
