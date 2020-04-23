import numpy as np
import move as mv
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import stats
import stats2 as ps
import copy
import os
import sys
import math
import bisect

def algoritmo_test(chain,d):
    g = chain
    theta = []
    angles = []
    radius = []
    for k in range(0,len(g)-3):
        g = chain.copy()
        point_a = np.array([g[k][0],g[k][1],g[k][2]])
        point_b = np.array([g[k+1][0],g[k+1][1],g[k+1][2]])
        point_c = np.array([g[k+2][0],g[k+2][1],g[k+2][2]])
        point_d = np.array([g[k+3][0],g[k+3][1],g[k+3][2]])
        
        ag,bg,cg = mv.f_correction(point_c,point_d,d)
        point_d[0] = ag
        point_d[1] = bg
        point_d[2] = cg
    
        ah,bh,ch,dh,eh,fh = mv.b2_correction(point_c,point_b,point_a,d)
        point_b[0] = ah
        point_b[1] = bh
        point_b[2] = ch
        point_a[0] = dh
        point_a[1] = eh
        point_a[2] = fh

        ag,bg,cg = mv.f_correction(point_b,point_a,d)
        point_a[0] = ag
        point_a[1] = bg
        point_a[2] = cg
        
        point_d = mv.traslation(point_d,-1*point_c)
        point_b = mv.traslation(point_b,-1*point_c)
        point_a = mv.traslation(point_a,-1*point_c)
        point_c = mv.traslation(point_c,-1*point_c)

        rot_matrix_x = mv.rotation_matrix_x(point_b,0,center=True)
        point_a = np.dot(rot_matrix_x,point_a)
        point_b = np.dot(rot_matrix_x,point_b)
        point_c = np.dot(rot_matrix_x,point_c)
        point_d = np.dot(rot_matrix_x,point_d)

        rot_matrix_y = mv.rotation_matrix_y(point_b,0,center=True)
        point_a = np.dot(rot_matrix_y,point_a)
        point_b = np.dot(rot_matrix_y,point_b)
        point_c = np.dot(rot_matrix_y,point_c)
        point_d = np.dot(rot_matrix_y,point_d)


        rot_matrix_z = mv.rotation_matrix_z(point_a,0,center=True)
        point_a = np.dot(rot_matrix_z,point_a)
        point_b = np.dot(rot_matrix_z,point_b)
        point_c = np.dot(rot_matrix_z,point_c)
        point_d = np.dot(rot_matrix_z,point_d)
        
        a,b,c = mv.cart2sph(point_d[0], point_d[1], point_d[2])
        radius.append(a)
        theta.append(b)
        angles.append(c)
        #angles.append(cart2sph(point_d[0], point_d[1], point_d[2])[2])
    
    return radius,theta,angles

def calculate_end_to_end223(data,L,w=0,get_chain=False):
    tmp_chain = []
    for i,v in enumerate(data):
        
        x = v[0]
        y = v[1]
        z = v[2]
        if i > 0 :
            a = int(np.abs(tmp_chain[-1][0]-x)/L)
            b = int(np.abs(tmp_chain[-1][1]-y)/L)
            c = int(np.abs(tmp_chain[-1][2]-z)/L)
            for m in range(0,a+1):
                if tmp_chain[-1][0]-x > L/2:
                    x = x + L
                if tmp_chain[-1][0]-x < -L/2:
                    x = x - L
            for m in range(0,b+1):
                if tmp_chain[-1][1]-y > L/2:
                    y = y + L
                if tmp_chain[-1][1]-y < -L/2:
                    y = y - L
            for m in range(0,c+1):
                if tmp_chain[-1][2]-z > L/2:
                    z = z + L
                if tmp_chain[-1][2]-z < -L/2:
                    z = z - L
            tmp_chain.append([x,y,z])
        else:
            tmp_chain.append([x,y,z])
    if w==0:
        ee = np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[-1]))**2
    else:
        ee = np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[w]))**2
    if get_chain == False:
        return(ee)
    else:
        return(tmp_chain)

def calculate_end_to_end_optimize(data,L,w,true_chain):
    tmp_chain = copy.deepcopy(true_chain)
    x = data[0]
    y = data[1]
    z = data[2]
    a = int(np.abs(tmp_chain[-1][0]-x)/L)
    b = int(np.abs(tmp_chain[-1][1]-y)/L)
    c = int(np.abs(tmp_chain[-1][2]-z)/L)
    for m in range(0,a+1):
        if tmp_chain[-1][0]-x > L/2:
            x = x + L
        if tmp_chain[-1][0]-x < -L/2:
            x = x - L
    for m in range(0,b+1):
        if tmp_chain[-1][1]-y > L/2:
            y = y + L
        if tmp_chain[-1][1]-y < -L/2:
            y = y - L
    for m in range(0,c+1):
        if tmp_chain[-1][2]-z > L/2:
            z = z + L
        if tmp_chain[-1][2]-z < -L/2:
            z = z - L
    tmp_chain.append([x,y,z])
    return(np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[-1]))**2, tmp_chain[-1] )


def calculate_end_to_end(data,L,w=0,get_chain=False):
    tmp_chain = []
    for i,v in enumerate(data):
        x = v[0]
        y = v[1]
        z = v[2]
        if i > 0 :
            a = int(np.abs(tmp_chain[-1][0]-x)/L)
            b = int(np.abs(tmp_chain[-1][1]-y)/L)
            c = int(np.abs(tmp_chain[-1][2]-z)/L)
            for m in range(0,a+1):
                if tmp_chain[-1][0]-x > L/2:
                    x = x + L
                if tmp_chain[-1][0]-x < -L/2:
                    x = x - L
            for m in range(0,b+1):
                if tmp_chain[-1][1]-y > L/2:
                    y = y + L
                if tmp_chain[-1][1]-y < -L/2:
                    y = y - L
            for m in range(0,c+1):
                if tmp_chain[-1][2]-z > L/2:
                    z = z + L
                if tmp_chain[-1][2]-z < -L/2:
                    z = z - L
            tmp_chain.append([x,y,z])
        else:
            tmp_chain.append([x,y,z])
        if w!=0 and w == i-1:
            break
    if w==0:
        ee = np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[-1]))**2
    else:
        ee = np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[w]))**2
        #print(ee)
    if get_chain == False:
        return(ee)
    else:
        return(tmp_chain)
               

def give_me_a_random_number(phi):
    hist, bins = np.histogram(phi, bins=300)
    bin_midpoints = bins[:-1] + np.diff(bins)/2
    cdf = np.cumsum(hist)
    cdf = cdf / cdf[-1]
    values = np.random.rand(1)
    value_bins = np.searchsorted(cdf, values)
    random_from_cdf = bin_midpoints[value_bins]
    return(random_from_cdf)
        
def calculate_end_to_end_mean(data,L,get_chain=False):
    ee = 0
    for k in data:
        tmp_chain = []
        for i,v in enumerate(k):
            x = v[0]
            y = v[1]
            z = v[2]
            if i > 0 :
                a = int(np.abs(tmp_chain[-1][0]-x)/L)
                b = int(np.abs(tmp_chain[-1][1]-y)/L)
                c = int(np.abs(tmp_chain[-1][2]-z)/L)
                for m in range(0,a+1):
                    if tmp_chain[-1][0]-x > L/2:
                        x = x + L
                    if tmp_chain[-1][0]-x < -L/2:
                        x = x - L
                for m in range(0,b+1):
                    if tmp_chain[-1][1]-y > L/2:
                        y = y + L
                    if tmp_chain[-1][1]-y < -L/2:
                        y = y - L
                for m in range(0,c+1):
                    if tmp_chain[-1][2]-z > L/2:
                        z = z + L
                    if tmp_chain[-1][2]-z < -L/2:
                        z = z - L
                tmp_chain.append([x,y,z])
            else:
                tmp_chain.append([x,y,z])
        ee += np.linalg.norm(np.array(tmp_chain[0])-np.array(tmp_chain[-1]))**2
    if get_chain == False:
        return(ee/len(data))
    else:
        return(tmp_chain)

def sigmoide(x,a):
    #return(a*x/np.sqrt(1+x**2))
    return(math.tanh(x*a))

def get_label(x,y,z,rcut2,L):
    ratio = math.ceil(L/rcut2)
    x = math.floor(x/rcut2)
    y = math.floor(y/rcut2)
    z = math.floor(z/rcut2)
    if x == ratio:
        x = x - 1
    if y == ratio:
        y = y - 1
    if z == ratio:
        z = z - 1
    return(z*ratio**2 + y*ratio + x)

def same_algorithm(N,M,phi,ee,tol,name,pbc=True):
    #N Numero de Cadenas
    #M Numero de moleculas
    Data = np.zeros((N,M,3))
    label = np.ones((N,M))*-1
    lo = 1.54
    Tho = 1.911
    ro = 0.761
    Mm = 14.02658
    Mol = 6.02214129*10**23
    L = ((N*M*Mm/(ro*Mol))**(1/3))*1e8
    print(L)
    L = float('%.4f'%(L))
    rcut = L/10
    rcut2 =  1 
    ratio = math.ceil(L/rcut2)
    valores = []
    
    
    j = 0
    #print("rcut",rcut)
    print("L:",L)
    while (j <= N-1):
        print("Cadena: " + str(j))
        ### Primera Molecula ###
        i=0
        s = True
        Data[j][0][0],Data[j][0][1],Data[j][0][2] = np.random.uniform(0,1)*L,np.random.uniform(0,1)*L,np.random.uniform(0,1)*L
        while( Data[j][0][0] >=L or Data[j][0][1]>=L or Data[j][0][2]>=L or Data[j][0][0] <= 0 or Data[j][0][1]<=0 or Data[j][0][2]<=0 or s):
            Data[j][0][0],Data[j][0][1],Data[j][0][2] = np.random.uniform(0,1)*L,np.random.uniform(0,1)*L,np.random.uniform(0,1)*L
            s = False
            x = math.floor(Data[j][0][0]/rcut2)
            y = math.floor(Data[j][0][1]/rcut2)
            z = math.floor(Data[j][0][2]/rcut2)
            for cnt1 in range(-2,3):
                for cnt2 in range(-2,3):
                    for cnt3 in range(-2,3):
                        x_a=x_b=y_a=y_b=z_a=z_b = False

                        tmp_x = x + cnt3
                        tmp_y = y + cnt2
                        tmp_z = z + cnt1

                        if tmp_x == ratio:
                            x_a = True
                            tmp_x = 0
                        if tmp_x == ratio+1:
                            x_a = True
                            tmp_x = 1
                        elif tmp_x == -2:
                            x_b = True
                            tmp_x = ratio-2
                        elif tmp_x == -1:
                            x_b = True
                            tmp_x = ratio-1
                        if tmp_y == ratio:
                            y_a = True
                            tmp_y = 0
                        elif tmp_y == ratio+1:
                            y_a = True
                            tmp_y = 1
                        elif tmp_y == -2:
                            y_b = True
                            tmp_y = ratio-2
                        elif tmp_y == -1:
                            y_b = True
                            tmp_y = ratio-1
                        if tmp_z == ratio:
                            z_a = True
                            tmp_z = 0
                        elif tmp_z == ratio+1:
                            z_a = True
                            tmp_z = 1
                        elif tmp_z == -2:
                            z_b = True
                            tmp_z = ratio-2
                        elif tmp_z == -1:
                            z_b = True
                            tmp_z = ratio-1

                        t_label = math.ceil(tmp_z*ratio**2 + tmp_y*ratio + tmp_x)
                        t_label = np.where(label==t_label)
                        points_in_space = np.array(Data[t_label])
                        if x_a == True:
                            points_in_space[:,0] = points_in_space[:,0] + L
                        if x_b == True:
                            points_in_space[:,0] = points_in_space[:,0] - L
                        if y_a == True:
                            points_in_space[:,1] = points_in_space[:,1] + L
                        if y_b == True:
                            points_in_space[:,1] = points_in_space[:,1] - L
                        if z_a == True:
                            points_in_space[:,2] = points_in_space[:,2] + L
                        if z_b == True:
                            points_in_space[:,2] = points_in_space[:,2] - L
                        if np.sum(np.sqrt(np.sum((points_in_space-np.array([Data[j][i][0],Data[j][i][1],Data[j][i][2]]))**2,axis=1))>rcut2) < len(points_in_space):
                            print(":c")
                            s = True  
        label[j,0] = get_label(Data[j][0][0],Data[j][0][1],Data[j][0][2],rcut2,L)
        ### Segunda Molecula ###
        i=1
        #print("Segunda Molecula de la cadena: " + str(j))
        Th = np.random.uniform(0,1)*np.pi
        #Ph = give_me_a_random_number(phi)
        Ph = np.random.uniform(0,2*np.pi,1)
        Data[j][1][0] = Data[j][0][0] + lo*np.sin(Th)*np.cos(Ph) 
        Data[j][1][1] = Data[j][0][1] + lo*np.sin(Th)*np.sin(Ph) 
        Data[j][1][2] = Data[j][0][2] + lo*np.cos(Th) 
        s = True
        while( Data[j][1][0] >=L or Data[j][1][1]>=L or Data[j][1][2]>=L or Data[j][1][0] <= 0 or Data[j][1][1]<=0 or Data[j][1][2]<=0 or s):
        #while(s):
            Th = np.random.rand(1)[0]*np.pi
            ph = give_me_a_random_number(phi)[0]
            Data[j][1][0] = Data[j][0][0] + lo*np.sin(Th)*np.cos(Ph) 
            Data[j][1][1] = Data[j][0][1] + lo*np.sin(Th)*np.sin(Ph) 
            Data[j][1][2] = Data[j][0][2] + lo*np.cos(Th)
            if pbc == True:
                if(Data[j][1][0] > L):
                    Data[j][1][0] = Data[j][1][0] - L
                elif(Data[j][1][0] < 0):
                    Data[j][1][0] = Data[j][1][0] + L
                if(Data[j][1][1] > L):
                    Data[j][1][1] = Data[j][1][1] - L
                elif(Data[j][1][1] < 0):
                    Data[j][1][1] = Data[j][1][1] + L
                if(Data[j][1][2] > L):
                    Data[j][1][2] = Data[j][1][2] - L
                elif(Data[j][1][2] < 0):
                    Data[j][1][2] = Data[j][1][2] + L
                
            s = False
            x = math.floor(Data[j][1][0]/rcut2)
            y = math.floor(Data[j][1][1]/rcut2)
            z = math.floor(Data[j][1][2]/rcut2)
            for cnt1 in range(-2,3):
                for cnt2 in range(-2,3):
                    for cnt3 in range(-2,3):
                        x_a=x_b=y_a=y_b=z_a=z_b = False

                        tmp_x = x + cnt3
                        tmp_y = y + cnt2
                        tmp_z = z + cnt1

                        if tmp_x == ratio:
                            x_a = True
                            tmp_x = 0
                        if tmp_x == ratio+1:
                            x_a = True
                            tmp_x = 1
                        elif tmp_x == -2:
                            x_b = True
                            tmp_x = ratio-2
                        elif tmp_x == -1:
                            x_b = True
                            tmp_x = ratio-1
                        if tmp_y == ratio:
                            y_a = True
                            tmp_y = 0
                        elif tmp_y == ratio+1:
                            y_a = True
                            tmp_y = 1
                        elif tmp_y == -2:
                            y_b = True
                            tmp_y = ratio-2
                        elif tmp_y == -1:
                            y_b = True
                            tmp_y = ratio-1
                        if tmp_z == ratio:
                            z_a = True
                            tmp_z = 0
                        elif tmp_z == ratio+1:
                            z_a = True
                            tmp_z = 1
                        elif tmp_z == -2:
                            z_b = True
                            tmp_z = ratio-2
                        elif tmp_z == -1:
                            z_b = True
                            tmp_z = ratio-1

                        t_label = math.ceil(tmp_z*ratio**2 + tmp_y*ratio + tmp_x)
                        t_label = np.where(label==t_label)
                        points_in_space = np.array(Data[t_label])
                        if x_a == True:
                            points_in_space[:,0] = points_in_space[:,0] + L
                        if x_b == True:
                            points_in_space[:,0] = points_in_space[:,0] - L
                        if y_a == True:
                            points_in_space[:,1] = points_in_space[:,1] + L
                        if y_b == True:
                            points_in_space[:,1] = points_in_space[:,1] - L
                        if z_a == True:
                            points_in_space[:,2] = points_in_space[:,2] + L
                        if z_b == True:
                            points_in_space[:,2] = points_in_space[:,2] - L
                        if np.sum(np.sqrt(np.sum((points_in_space-np.array([Data[j][i][0],Data[j][i][1],Data[j][i][2]]))**2,axis=1))>rcut2) < len(points_in_space):
                            print(":c")
                            s = True
        label[j,1] = get_label(Data[j][1][0],Data[j][1][1],Data[j][1][2],rcut2,L)           
        ### Tercera Molecula ###
        i=2
        #print("Tercera Molecula de la cadena: " + str(j))
        
        temporal_0 = mv.traslation(Data[j][0],-Data[j][1])
        rot_matrix_x = mv.rotation_matrix_x(temporal_0,0,center=True)
        temporal_1 = np.dot(rot_matrix_x,temporal_0)
        rot_matrix_y = mv.rotation_matrix_y(temporal_1,0,center=True)
        #r_phi = give_me_a_random_number(phi)[0]
        r_phi = np.random.uniform(0,2*np.pi,1)
        Data[j][2][0] = lo*np.sin(Tho)*np.cos(r_phi) 
        Data[j][2][1] = lo*np.sin(Tho)*np.sin(r_phi) 
        Data[j][2][2] = lo*np.cos(Tho) 
        Data[j][2] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),Data[j][2])),Data[j][1])
        s = True
        while( Data[j][2][0] >=L or Data[j][2][1]>=L or Data[j][2][2]>=L or Data[j][2][0] <= 0 or Data[j][2][1]<=0 or Data[j][2][2]<=0 or s):
            r_phi = np.random.uniform(0,2*np.pi,1)
            Data[j][2][0] = lo*np.sin(Tho)*np.cos(r_phi) 
            Data[j][2][1] = lo*np.sin(Tho)*np.sin(r_phi) 
            Data[j][2][2] = lo*np.cos(Tho) 
            Data[j][2] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),Data[j][2])),Data[j][1])
            
            if pbc == True:
                if(Data[j][2][0] > L):
                    Data[j][2][0] = Data[j][2][0] - L
                elif(Data[j][2][0] < 0):
                    Data[j][2][0] = Data[j][2][0] + L
                if(Data[j][2][1] > L):
                    Data[j][2][1] = Data[j][2][1] - L
                elif(Data[j][2][1] < 0):
                    Data[j][2][1] = Data[j][2][1] + L
                if(Data[j][2][2] > L):
                    Data[j][2][2] = Data[j][2][2] - L
                elif(Data[j][2][2] < 0):
                    Data[j][2][2] = Data[j][2][2] + L
                
            s = False
            x = math.floor(Data[j][2][0]/rcut2)
            y = math.floor(Data[j][2][1]/rcut2)
            z = math.floor(Data[j][2][2]/rcut2)
            for cnt1 in range(-2,3):
                for cnt2 in range(-2,3):
                    for cnt3 in range(-2,3):
                        x_a=x_b=y_a=y_b=z_a=z_b = False

                        tmp_x = x + cnt3
                        tmp_y = y + cnt2
                        tmp_z = z + cnt1

                                
                        if tmp_x == ratio:
                            x_a = True
                            tmp_x = 0
                        if tmp_x == ratio+1:
                            x_a = True
                            tmp_x = 1
                        elif tmp_x == -2:
                            x_b = True
                            tmp_x = ratio-2
                        elif tmp_x == -1:
                            x_b = True
                            tmp_x = ratio-1
                        if tmp_y == ratio:
                            y_a = True
                            tmp_y = 0
                        elif tmp_y == ratio+1:
                            y_a = True
                            tmp_y = 1
                        elif tmp_y == -2:
                            y_b = True
                            tmp_y = ratio-2
                        elif tmp_y == -1:
                            y_b = True
                            tmp_y = ratio-1
                        if tmp_z == ratio:
                            z_a = True
                            tmp_z = 0
                        elif tmp_z == ratio+1:
                            z_a = True
                            tmp_z = 1
                        elif tmp_z == -2:
                            z_b = True
                            tmp_z = ratio-2
                        elif tmp_z == -1:
                            z_b = True
                            tmp_z = ratio-1

                        t_label = math.ceil(tmp_z*ratio**2 + tmp_y*ratio + tmp_x)
                        t_label = np.where(label==t_label)
                        points_in_space = np.array(Data[t_label])
                        if x_a == True:
                            points_in_space[:,0] = points_in_space[:,0] + L
                        if x_b == True:
                            points_in_space[:,0] = points_in_space[:,0] - L
                        if y_a == True:
                            points_in_space[:,1] = points_in_space[:,1] + L
                        if y_b == True:
                            points_in_space[:,1] = points_in_space[:,1] - L
                        if z_a == True:
                            points_in_space[:,2] = points_in_space[:,2] + L
                        if z_b == True:
                            points_in_space[:,2] = points_in_space[:,2] - L
                        if np.sum(np.sqrt(np.sum((points_in_space-np.array([Data[j][i][0],Data[j][i][1],Data[j][i][2]]))**2,axis=1))>rcut2) < len(points_in_space):
                            print(":c")
                            s = True

        label[j,2] = get_label(Data[j][2][0],Data[j][2][1],Data[j][2][2],rcut2,L)
        ### Moleculas Restantes ###
        i = 3
        intentos = 0
        angulos_test = []
        cadena_coeficiente = []
        true_chain = None
        while(i <=M-1):
            if i%100==0:
                print(str(i)+" Molecula de la cadena: " + str(j))
            drop = 1
            emn = 1     
            point_a = copy.deepcopy(Data[j][i-3])
            point_b = copy.deepcopy(Data[j][i-2])
            point_c = copy.deepcopy(Data[j][i-1])
            
            if point_c[0]-point_b[0] > L/2:
                point_b[0] = point_b[0] + L
            elif point_c[0]-point_b[0] < -L/2:
                point_b[0] = point_b[0] - L
                
            if point_c[1]-point_b[1] > L/2:
                point_b[1] = point_b[1] + L
            elif point_c[1]-point_b[1] < -L/2:
                point_b[1] = point_b[1] - L
                
            if point_c[2]-point_b[2] > L/2:
                point_b[2] = point_b[2] + L
            elif point_c[2]-point_b[2] < -L/2:
                point_b[2] = point_b[2] - L

            if point_b[0]-point_a[0] > L/2:
                point_a[0] = point_a[0] + L
            elif point_b[0]-point_a[0] < -L/2:
                point_a[0] = point_a[0] - L
                
            if point_b[1]-point_a[1] > L/2:
                point_a[1] = point_a[1] + L
            elif point_b[1]-point_a[1] < -L/2:
                point_a[1] = point_a[1] - L
            if point_b[2]-point_a[2] > L/2:
                point_a[2] = point_a[2] + L
            elif point_b[2]-point_a[2] < -L/2:
                point_a[2] = point_a[2] - L
                
            
            temporal_0 = mv.traslation(point_b,-point_c)
            temporala_0 = mv.traslation(point_a,-point_c)
            rot_matrix_x = mv.rotation_matrix_x(temporal_0,0,center=True)
            temporal_1 = np.dot(rot_matrix_x,temporal_0)
            temporala_1 = np.dot(rot_matrix_x,temporala_0)
            rot_matrix_y = mv.rotation_matrix_y(temporal_1,0,center=True)
            temporal_2 = np.dot(rot_matrix_y,temporal_1)
            temporala_2 = np.dot(rot_matrix_y,temporala_1)
            rot_matrix_z = mv.rotation_matrix_z(temporala_2,0,center=True)
            
            r_phi = np.random.uniform(0,2*np.pi,1)[0]
            Data[j][i][0] = lo*np.sin(Tho)*np.cos(r_phi) 
            Data[j][i][1] = lo*np.sin(Tho)*np.sin(r_phi) 
            Data[j][i][2] = lo*np.cos(Tho)
            Data[j][i] = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),np.dot(np.linalg.inv(rot_matrix_z),Data[j][i]))),Data[j][i-1]) 
            s = True
            #angulos_test.append(r_phi)
            while( Data[j][i][0] >=L or Data[j][i][1]>=L or Data[j][i][2]>=L or Data[j][i][0] <= 0 or Data[j][i][1]<=0 or Data[j][i][2]<=0 or s):
                temp_10 = []
                temp_10_index = []
                temp_10_end_to_end = []
                temp_angle = []
                intentos = intentos + 1
                flag = True
                temporal_angles = []
                temp_final_coord = []
                counter = 0
                #angulos_test = 
                angulos_test.sort()
                for temporal_angles in range(0,25):
                    r_phi = np.random.uniform(0,2*np.pi,1)[0]
                    tmp_molecule = [lo*np.sin(Tho)*np.cos(r_phi),
                                    lo*np.sin(Tho)*np.sin(r_phi),
                                    lo*np.cos(Tho)]
                    tmp_molecule = mv.traslation(np.dot(np.linalg.inv(rot_matrix_x),np.dot(np.linalg.inv(rot_matrix_y),np.dot(np.linalg.inv(rot_matrix_z),tmp_molecule))),Data[j][i-1])
                    temp_angle.append(r_phi)
                    if pbc == True:
                        if tmp_molecule[0] > L:
                            tmp_molecule[0] = tmp_molecule[0] - L
                        elif tmp_molecule[0] < 0:
                            tmp_molecule[0] = tmp_molecule[0] + L           
                        if tmp_molecule[1] > L:
                            tmp_molecule[1] = tmp_molecule[1] - L
                        elif tmp_molecule[1] < 0:
                            tmp_molecule[1] = tmp_molecule[1] + L
                        if tmp_molecule[2] > L:
                            tmp_molecule[2] = tmp_molecule[2] - L
                        elif tmp_molecule[2] < 0:
                            tmp_molecule[2] = tmp_molecule[2] + L
                    temp_10.append(tmp_molecule)
                    Data[j][i][0] = tmp_molecule[0]
                    Data[j][i][1] = tmp_molecule[1]
                    Data[j][i][2] = tmp_molecule[2]
                    #angulos_test[-1] = r_phi
                    if true_chain == None:
                        bb_w = ((i+1)*ee/M - calculate_end_to_end(Data[j],L,i,False))**2
                    else:
                        eec, coord = calculate_end_to_end_optimize(Data[j][i],L,i,true_chain)
                        bb_w = ((i+1)*ee/M - eec)**2
                        temp_final_coord.append(coord)
                    new_list = angulos_test[:]
                    bb_w2 = (ps.ks_2samp(phi,new_list,r_phi)[1])
                    temp_10_end_to_end.append(bb_w)
                    temp_10_index.append(bb_w2)

                temp_10_end_to_end = np.array(temp_10_end_to_end)/((i+1)*ee/M)**2
                cadena_coeficiente = (1+temp_10_end_to_end)/np.array(temp_10_index)
                temp_10_end_to_end = np.ndarray.tolist((1-(1+temp_10_end_to_end)/np.array(temp_10_index))**2)
                if true_chain == None:
                    check = [[x,y,z,r,p] for y,x,z,r,p in sorted(zip(temp_10_end_to_end,temp_10_index,temp_angle,temp_10,cadena_coeficiente))]
                else:
                    check = [[x,y,z,r,p,u] for y,x,z,r,p,u in sorted(zip(temp_10_end_to_end,temp_10_index,temp_angle,temp_10,cadena_coeficiente,temp_final_coord))]
                for molecule_index in check:
                    s = False
                    Data[j][i][0] = molecule_index[3][0]
                    Data[j][i][1] = molecule_index[3][1]
                    Data[j][i][2] = molecule_index[3][2]
                
                    x = math.floor(Data[j][i][0]/rcut2)
                    y = math.floor(Data[j][i][1]/rcut2)
                    z = math.floor(Data[j][i][2]/rcut2)
                    
                    for cnt1 in range(-2,3):
                        for cnt2 in range(-2,3):
                            for cnt3 in range(-2,3):
                                x_a=x_b=y_a=y_b=z_a=z_b = False
                                
                                tmp_x = x + cnt3
                                tmp_y = y + cnt2
                                tmp_z = z + cnt1
                                
                                if tmp_x == ratio:
                                    x_a = True
                                    tmp_x = 0
                                if tmp_x == ratio+1:
                                    x_a = True
                                    tmp_x = 1
                                elif tmp_x == -2:
                                    x_b = True
                                    tmp_x = ratio-2
                                elif tmp_x == -1:
                                    x_b = True
                                    tmp_x = ratio-1
                                if tmp_y == ratio:
                                    y_a = True
                                    tmp_y = 0
                                elif tmp_y == ratio+1:
                                    y_a = True
                                    tmp_y = 1
                                elif tmp_y == -2:
                                    y_b = True
                                    tmp_y = ratio-2
                                elif tmp_y == -1:
                                    y_b = True
                                    tmp_y = ratio-1
                                if tmp_z == ratio:
                                    z_a = True
                                    tmp_z = 0
                                elif tmp_z == ratio+1:
                                    z_a = True
                                    tmp_z = 1
                                elif tmp_z == -2:
                                    z_b = True
                                    tmp_z = ratio-2
                                elif tmp_z == -1:
                                    z_b = True
                                    tmp_z = ratio-1
                                    
                                t_label = math.ceil(tmp_z*ratio**2 + tmp_y*ratio + tmp_x)
                                t_label = np.where(label==t_label)
                                points_in_space = np.array(Data[t_label])
                                if x_a == True:
                                    points_in_space[:,0] = points_in_space[:,0] + L
                                if x_b == True:
                                    points_in_space[:,0] = points_in_space[:,0] - L
                                if y_a == True:
                                    points_in_space[:,1] = points_in_space[:,1] + L
                                if y_b == True:
                                    points_in_space[:,1] = points_in_space[:,1] - L
                                if z_a == True:
                                    points_in_space[:,2] = points_in_space[:,2] + L
                                if z_b == True:
                                    points_in_space[:,2] = points_in_space[:,2] - L
                                v = np.sqrt(np.sum((points_in_space-np.array([Data[j][i][0],Data[j][i][1],Data[j][i][2]]))**2,axis=1))
                                #print(points_in_space)
                                #print(v)
                                if np.sum(np.sqrt(np.sum((points_in_space-np.array([Data[j][i][0],Data[j][i][1],Data[j][i][2]]))**2,axis=1))>rcut2) < len(points_in_space):
                                    s = True
                    if s == False:
                        valores.append(molecule_index[4])
                        break


                intentos = intentos + 1
                if intentos > 100:
                    emn = 0
                    drop = drop + np.random.randint(1,99)
                    break
            if emn == 1:
                #angulos_test[(i-3)] = molecule_index[2]
                angulos_test.append(molecule_index[2])
                label[j,i] = get_label(Data[j][i][0],Data[j][i][1],Data[j][i][2],rcut2,L)
                i = i + 1
                if true_chain == None:
                    true_chain = calculate_end_to_end(Data[j],L,i,True)
                else:
                    true_chain.append(molecule_index[-1])
                if drop > 0:
                    drop = drop - 1
            if emn == 0:
                emn = 1
                i = i - drop
                #angulos_test = angulos_test[0:(i-3)]
                if true_chain != None:
                    true_chain = true_chain[0:i]
                if i < 3:
                    break
            if drop > 100:
                drop = 0
                break
            intentos = 0
        if intentos < 10:
            target = calculate_end_to_end(Data[j],L)
            np.save("restart_"+name+".npy",Data)
            print(target)
            print(stats.ks_2samp(phi,angulos_test)[1])
            j= j + 1
    return(Data,label)


if __name__ == '__main__':
	nombre_salida = "data"+sys.argv[1]+".npy"
	print(nombre_salida)
	g = open("/home/pibarra/generador_cadenas/repo_memoria/generador_cadenas/c.txt").read().split("\n")
	data = {}
	for c in g:
    		temporal = c.split(" ")
    		if len(temporal) > 1:
        		data[int(temporal[0])] = {"mol":int(temporal[1]),
                		                  "x":float(temporal[2]),
                        		          "y":float(temporal[3]),
                                		  "z":float(temporal[4])}
	atoms = []
	temporal = []
	mol = 1
	for i in range(1,len(data)+1):
    		if data[i]['mol'] != mol:
        		mol = data[i]['mol']
        		atoms.append(temporal)
        		temporal = []
    		temporal.append(np.array([data[i]["x"],data[i]["y"],data[i]["z"]]))
	atoms.append(temporal)
	d = [3.2188079615971077e+00,1.4197189203783057e+02]
	r = []
	theta = []
	phi = []
	for t,k in enumerate(atoms):
    		a,b,c = algoritmo_test(k,d)
    		if np.sum(np.abs(np.array(a)) > 3.0):
        		print(t+1)
        		print(np.where(np.abs(np.array(a)) > 3.0))
        		#print(a)
    		r = r + a
    		theta = theta + b
    		phi = phi + c
	phi = np.sort(phi)
	datita,_ = same_algorithm(int(sys.argv[4]),int(sys.argv[3]),phi,int(sys.argv[2]),0.05,nombre_salida,True)  
	np.save(nombre_salida,datita)

