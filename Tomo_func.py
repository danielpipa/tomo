import numpy
import scipy.ndimage


#-----------------------------------------------------------------
#Shift operator
#Receives layers and shifts them according with angles and altitudes
#-----------------------------------------------------------------
def S(layersXY,sizes, alt_p,angles):
    N_layers = len(alt_p);   #number of layers
    sizeG = sizes[0];    #size of Ground Layer (and wfs)
    for i in range(1,N_layers):  #skip ground layer
        x_shift = numpy.tan(angles[0]*numpy.pi/648000)*alt_p[i];
        y_shift = numpy.tan(angles[1]*numpy.pi/648000)*alt_p[i];
        #zero padding to size of layer
        temp = layersXY[i];
        tempR = scipy.ndimage.interpolation.shift(temp.real,[x_shift,y_shift])    #shift real
        tempI = scipy.ndimage.interpolation.shift(temp.imag,[x_shift,y_shift])    #shift imaginary
        temp = tempR+1j*tempI;
        crop = int((len(temp)-sizeG)/2);
        if crop != 0:
            layersXY[i] = temp[crop:-crop,crop:-crop];

    return layersXY;

#-----------------------------------------------------------------
#Shift operator
#Receives layers and shifts them according with angles and altitudes
#-----------------------------------------------------------------
def St(layersXY,sizes, alt_p,angles):
    N_layers = len(alt_p);   #number of layers
    sizeG = sizes[0];    #size of Ground Layer (and wfs)
    for i in range(1,N_layers):  #skip ground layer
        x_shift = numpy.tan(angles[0]*numpy.pi/648000)*alt_p[i];
        y_shift = numpy.tan(angles[1]*numpy.pi/648000)*alt_p[i];
        #zero padding to size of layer
        temp = numpy.pad(layersXY[i],(int((sizes[i]-sizeG)/2), int((sizes[i]-sizeG)/2)),'constant');
        tempR = scipy.ndimage.interpolation.shift(temp.real,[-x_shift,-y_shift])    #shift real
        tempI = scipy.ndimage.interpolation.shift(temp.imag,[-x_shift,-y_shift])    #shift imaginary
        layersXY[i] = tempR+1j*tempI;

    return layersXY;

#-----------------------------------------------------------------
#Projection Operator
#Sums all layers
#-----------------------------------------------------------------
def P(layersXY):
    N_layers = len(layersXY);   #number of layers
    y = layersXY[0];
    for i in range(1,N_layers):
        y = y+layersXY[i];
    return y;


#-----------------------------------------------------------------
#Transpose Projection Operator
#Copy data for all layers
#-----------------------------------------------------------------
def Pt(data,sizes):
    N_layers = len(sizes);   #number of layers
    y = numpy.array(range(N_layers), dtype=object)
    for i in range(0,N_layers):
        y[i] = data;
    return y;

#-----------------------------------------------------------------
# Pupil function
#---Cria circulo da pupila de tamanho size e raio radius
#-----------------------------------------------------------------
def HH(data, sizes, alt_p, WFS_angles, t):

    if(t):
        Tn_wfs = len(WFS_angles);    #number of WFS
        n_subaps = sizes[0]**2;
        WFS = numpy.array(range(Tn_wfs), dtype=object);
        for wfs_n in range(0,Tn_wfs):
            WFS[wfs_n] = numpy.reshape(data[wfs_n*n_subaps:wfs_n*n_subaps+n_subaps],[sizes[0],sizes[0]]);

        Tn_layers = len(sizes);
        out = numpy.zeros(numpy.sum(sizes**2))
        for wfs_n in range(0,Tn_wfs):
            temp = St(Pt(WFS[wfs_n],sizes),sizes,alt_p,WFS_angles[wfs_n,:]);
            l_n=0;
            temp1 = numpy.reshape(temp[l_n],[-1,1]);
            for l_n in range(1,Tn_layers):
                temp1 = numpy.append(temp1, numpy.reshape(temp[l_n],[-1,1]));
            out = out + temp1;
    else:
        Tn_layers = len(sizes);     #Number of layers
        layersXY = numpy.array(range(Tn_layers), dtype=object)
        layersXY[0] = numpy.reshape(data[0:sizes[0]**2],[sizes[0],sizes[0]]);
        for l_n in range(1,Tn_layers):
            layersXY[l_n] = numpy.reshape(data[l_n*sizes[l_n-1]**2:l_n*sizes[l_n-1]**2+sizes[l_n]**2],[sizes[l_n],sizes[l_n]]);
        Tn_wfs = len(WFS_angles);   #Number of WFS

        wfs_n = 0;
        out = P(S(layersXY,sizes,alt_p,WFS_angles[wfs_n,:]));  #%output is the WFSs slopes
        for wfs_n in range(1,Tn_wfs):
            temp = P(S(layersXY,sizes,alt_p,WFS_angles[wfs_n,:]));
            out = numpy.append(out,temp);

    return out;
