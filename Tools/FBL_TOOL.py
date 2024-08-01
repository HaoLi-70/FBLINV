import sys, os, struct

import numpy as np

################################################################################

def writemodel(filename, data, Rrange0=None, rhotype=0, btype=1, \
        vtype=-1):
    
    if (data['r'].ndim != 1): 
      print('wrong dimension of r')
      return {'check': False}
    if (data['theta'].ndim != 1): 
      print('wrong dimension of theta')
      return {'check': False}
    if (data['phi'].ndim != 1): 
      print('wrong dimension of phi')
      return {'check': False}
    if (data['temp'].ndim != 3): 
      print('wrong dimension of temp')
      return {'check': False}
    if (data['rho'].ndim != 3):
      print('wrong dimension of rho') 
      return {'check': False}
    if (data['Br'].ndim != 3): 
      print('wrong dimension of Br')
      return {'check': False}
    if (data['Bt'].ndim != 3): 
      print('wrong dimension of Btheta')
      return {'check': False}
    if (data['Bp'].ndim != 3): 
      print('wrong dimension of Bphi')
      return {'check': False}
    


    nr = data['r'].shape[0]
    ntheta = data['theta'].shape[0]
    nphi = data['phi'].shape[0]

    if(Rrange0 is None):
      Rrange = [0,nr]
    else:
      Rrange = Rrange0

    shapes = data['temp'].shape
    if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
      print('wrong tenp shape')
      print('nr = ', nr, shapes[2])
      print('ntheta = ', ntheta, shapes[1])
      print('nphi = ', nphi, shapes[0])
      return {'check': False}
    shapes = data['rho'].shape
    if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
      print('wrong rho shape')
      print('nr = ', nr, shapes[2])
      print('ntheta = ', ntheta, shapes[1])
      print('nphi = ', nphi, shapes[0])
      return {'check': False}
    shapes = data['Br'].shape
    if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
      print('wrong Br shape')
      print('nr = ', nr, shapes[2])
      print('ntheta = ', ntheta, shapes[1])
      print('nphi = ', nphi, shapes[0])
      return {'check': False}
    shapes = data['Bt'].shape
    if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
      print('wrong Btheta shape')
      print('nr = ', nr, shapes[2])
      print('ntheta = ', ntheta, shapes[1])
      print('nphi = ', nphi, shapes[0])
      return {'check': False}
    shapes = data['Bp'].shape
    if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
      print('wrong Bphi shape')
      print('nr = ', nr, shapes[2])
      print('ntheta = ', ntheta, shapes[1])
      print('nphi = ', nphi, shapes[0])
      return {'check': False}

    nrout = Rrange[1]-Rrange[0]
    
    dataw = {}

    dataw['r'] = np.array(data['r'][Rrange[0]:Rrange[1]])
    print('R range = ', dataw['r'][0], dataw['r'][nrout-1])
    dataw['theta'] = np.array(data['theta'][:])
    dataw['phi']= np.array(data['phi'][:])
    
    dataw['Br'] = np.array(data['Br'][:,:,Rrange[0]:Rrange[1]])
    dataw['Bt'] = np.array(data['Bt'][:,:,Rrange[0]:Rrange[1]])
    dataw['Bp'] = np.array(data['Bp'][:,:,Rrange[0]:Rrange[1]])
    dataw['temp'] = np.array(data['temp'][:,:,Rrange[0]:Rrange[1]])
    dataw['rho'] = np.array(data['rho'][:,:,Rrange[0]:Rrange[1]])

    f = open(filename,'wb')

    f.write(b'fbmd')

    f.write(struct.pack('I',nrout))

    f.write(struct.pack('I',ntheta))

    f.write(struct.pack('I',nphi))

    arr = np.array(dataw['r'],dtype=np.float64)
    f.write(arr.tobytes())

    arr = np.array(dataw['theta'],dtype=np.float64)
    f.write(arr.tobytes())

    arr = np.array(dataw['phi'],dtype=np.float64)
    f.write(arr.tobytes())

    for ir in range(nrout):
      for itheta in range(ntheta):
        arr = np.array(dataw['temp'][:,itheta,ir],dtype=np.float64)
        f.write(arr.tobytes())

    # 0: log10 nh; 1: nh; 2: log10 ne; 3: ne
    f.write(struct.pack('I',rhotype))

    if(rhotype==0 or rhotype==2):
      dataw['lgrho'] = np.log10(dataw['rho'])
      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['lgrho'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())
    elif(rhotype==1 or rhotype==3):
      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['rho'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())

    # 0: br bt bp; 1: b thetab phib 
    f.write(struct.pack('I',btype))

    if(btype==0):
      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['Br'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())

      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['Bt'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())

      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['Bp'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())
    else:

      dataw['B'] = np.sqrt(dataw['Br']*dataw['Br'] \
          +dataw['Bt']*dataw['Bt']+dataw['Bp']*dataw['Bp'])
      dataw['thetaB'] = np.arccos(dataw['Br']/dataw['B'])
      dataw['phiB'] = np.arctan2(dataw['Bp'],dataw['Bt'])

      for i in range(nphi):
        for j in range(ntheta):
          for k in range(nrout):
            if(dataw['phiB'][i,j,k]<0):
              dataw['phiB'][i,j,k] = dataw['phiB'][i,j,k] \
                  +2.0*np.pi
      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['B'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())

      for ir in range(nrout):
        for itheta in range(ntheta):
          arr = np.array(dataw['thetaB'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())

      for ir in range(nrout):
        for itheta in range(ntheta):        
          arr = np.array(dataw['phiB'][:,itheta,ir],dtype=np.float64)
          f.write(arr.tobytes())

    if vtype>=0 and "Vr" in data:

      if (data['Vr'].ndim != 3): 
        f.close() 
        print('wrong dimension of Vr')
        return {'check': True}
      if (data['Vt'].ndim != 3):
        f.close()  
        print('wrong dimension of Vtheta')
        return {'check': True}
      if (data['Vp'].ndim != 3): 
        f.close() 
        print('wrong dimension of Vphi')
        return {'check': True}

      shapes = data['Vr'].shape
      if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
        f.close() 
        print('wrong Vr shape')
        print('nr = ', nr, shapes[2])
        print('ntheta = ', ntheta, shapes[1])
        print('nphi = ', nphi, shapes[0])
        return {'check': True}
      shapes = data['Vt'].shape
      if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
        f.close() 
        print('wrong Vtheta shape')
        print('nr = ', nr, shapes[2])
        print('ntheta = ', ntheta, shapes[1])
        print('nphi = ', nphi, shapes[0])        
        return {'check': True}
      shapes = data['Vp'].shape
      if (shapes[0]!=nphi or shapes[1]!=ntheta or shapes[2]!=nr): 
        f.close() 
        print('wrong Vphi shape')
        print('nr = ', nr, shapes[2])
        print('ntheta = ', ntheta, shapes[1])
        print('nphi = ', nphi, shapes[0])        
        return {'check': True}

      # 0: Vr Vt Vp; 1: V thetaV phiV 
      f.write(struct.pack('I',vtype))

      dataw['Vr'] = np.array(data['Vr'][:,:,Rrange[0]:Rrange[1]])
      dataw['Vt'] = np.array(data['Vt'][:,:,Rrange[0]:Rrange[1]])
      dataw['Vp'] = np.array(data['Vp'][:,:,Rrange[0]:Rrange[1]])

      if(vtype==0):
        for ir in range(nrout):
          for itheta in range(ntheta):
            arr = np.array(dataw['Vr'][:,itheta,ir],dtype=np.float64)
            f.write(arr.tobytes())

        for ir in range(nrout):
          for itheta in range(ntheta):
            arr = np.array(dataw['Vt'][:,itheta,ir],dtype=np.float64)
            f.write(arr.tobytes())

        for ir in range(nrout):
          for itheta in range(ntheta):
            arr = np.array(dataw['Vp'][:,itheta,ir],dtype=np.float64)
            f.write(arr.tobytes())

      elif(vtype==1):

        dataw['V'] = np.sqrt(dataw['Vr']*dataw['Vr'] \
            +dataw['Vt']*dataw['Vt']+dataw['Vp']*dataw['Vp'])
        dataw['thetaV'] = np.arccos(dataw['Vr']/dataw['V'])
        dataw['phiV'] = np.arctan2(dataw['Vp'],dataw['Vt'])

        for i in range(nphi):
          for j in range(ntheta):
            for k in range(nrout):
              if(dataw['phiV'][i,j,k]<0):
                dataw['phiV'][i,j,k] = dataw['phiV'][i,j,k] \
                    +2.0*np.pi
        for ir in range(nrout):
          for itheta in range(ntheta):
            arr = np.array(dataw['V'][:,itheta,ir],dtype=np.float64)
            f.write(arr.tobytes())

        for ir in range(nrout):
          for itheta in range(ntheta):
            arr = np.array(dataw['thetaV'][:,itheta,ir],dtype=np.float64)
            f.write(arr.tobytes())

        for ir in range(nrout):
          for itheta in range(ntheta):        
            arr = np.array(dataw['phiV'][:,itheta,ir],dtype=np.float64)
            f.write(arr.tobytes())

    f.close() 
    
    return {'check':True}

################################################################################

def readmodel(filename):

        
    try:
      f = open(filename,'rb')
    except:
      print('can not open '+filename)
      f.close() 
      return {'check': False}
    
    try:
      bytes = f.read(4)
      print(bytes)

      bytes = f.read(4*3)

      arr = np.array(struct.unpack('i'*3,bytes))


      dim1 = arr[0]
      dim2 = arr[1]
      dim3 = arr[2]

      print(arr)


      #bytes = f.read(8)

      #print('a=',struct.unpack('d',bytes)[0])
      res = {}
      bytes = f.read(8*dim1)
      res['r'] = np.array(struct.unpack('d'*dim1,bytes))
      bytes = f.read(8*dim2)
      res['theta'] = np.array(struct.unpack('d'*dim2,bytes))
      bytes = f.read(8*dim3)
      res['phi'] = np.array(struct.unpack('d'*dim3,bytes))

      bytes = f.read(8*dim1*dim2*dim3)
      res['temp'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))

      bytes = f.read(4)
      res['rhotype'] = np.array(struct.unpack('i',bytes))
      
      bytes = f.read(8*dim1*dim2*dim3)
      res['rho'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))

      bytes = f.read(4)
      res['btype'] = np.array(struct.unpack('i',bytes))

      if(res['btype']==0):
        bytes = f.read(8*dim1*dim2*dim3)
        res['Br'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['Btheta'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['Bphi'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
      else:
        bytes = f.read(8*dim1*dim2*dim3)
        res['B'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['ThetaB'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['PhiB'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))

      bytes = f.read(4)
      res['vtype'] = np.array(struct.unpack('i',bytes))

      if(res['vtype']==0):
        bytes = f.read(8*dim1*dim2*dim3)
        res['Vr'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['Vtheta'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['Vphi'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
      elif (res['vtype']==1):
        bytes = f.read(8*dim1*dim2*dim3)
        res['V'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['ThetaV'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
  
        bytes = f.read(8*dim1*dim2*dim3)
        res['PhiV'] = np.reshape(np.array(struct.unpack('d'*dim1*dim2*dim3,bytes)),(dim1,dim2,dim3))
        
      f.close() 
      return res
    except:
      
      print("error in reading " + filename)
      f.close() 
      return {'check': False}

################################################################################

################################################################################

def readfbl(filename):
    
    try:
      f = open(filename,'rb')
    except:
      print('can not open '+filename)
      f.close() 
      return
    
    try:
        
      bytes = f.read(4)
      print(bytes)

      bytes = f.read(4*4)
      arr = np.array(struct.unpack('i'*4,bytes))
      print(arr)

      bytes = f.read(8*4)
      arr2 = np.array(struct.unpack('d'*4,bytes))
      print(arr2)

      num = arr[2]*3+arr[3]*2
      nn = num*arr[1]*arr[0]

      data=np.zeros((num,arr[1],arr[0]),dtype=np.float64)

      print(8*nn)


      for i in range(arr[0]):
        for j in range(arr[1]):
          bytes = f.read(8*num)
          data[0:num,j,i] = np.array(struct.unpack('d'*num,bytes))
          #print(data[0:num,j,i])
      f.close() 
      return data
    except:
      
      print("error writing " + filename)
      f.close() 
      return {'check': False}
    
################################################################################    

def readcoeff(filename):
    
    try:
      f = open(filename,'rb')
    except:
      print('can not open '+filename)
      f.close() 
      return {'check': False}
    
    try:
        
      f = open(filename,'rb')
      bytes = f.read(4)
      print(bytes)
        
      res = {'check': True}

      bytes = f.read(4*6)

      arr = np.array(struct.unpack('i'*6,bytes))
      res['inv'] = arr[0:4]
      res['btype'] = arr[4]
      res['rhotype'] = arr[5]
  
      print('a',arr)
      bytes = f.read(8)
      res['Rbound'] = np.array(struct.unpack('d',bytes))

      N = np.zeros(4,dtype=np.int16)
      L = np.zeros(4,dtype=np.int16)

      print(N,L)
      symble = ['T', 'nH', 'ThetaB', 'PhiB']

      sub = 20
      for i in range(4):
        if(arr[i]!=0):
          bytes = f.read(4)
          N[i] = struct.unpack('i',bytes)[0]
          res[symble[i]+'_N'] = N[i]
          bytes = f.read(4)
          L[i] = struct.unpack('i',bytes)[0]
          res[symble[i]+'_L'] = L[i]
          Num = N[i]*(L[i]+1)*(L[i]+1)
          print('coe',i,N[i],L[i],Num)

          bytes = f.read(8*Num)
          res[symble[i]] = np.reshape(np.array(struct.unpack('d'*Num, \
              bytes)),(N[i],(L[i]+1)*(L[i]+1)))
          sub = sub+8*Num*8

      f.close() 
      return res
    except:
      
      print("error reading " + filename)
      f.close() 
      return {'check': False}
    
################################################################################

def readresult(filename):
      
          
    try:
      f = open(filename,'rb')
    except:
      print('can not open '+filename)
      f.close() 
      return {'check': False}
      
    try:

      res = {'check': True}

      bytes = f.read(4)


      match bytes:
      
        case b'syth':
      
              
          bytes = f.read(4*3)
          arr = np.array(struct.unpack('i'*3,bytes))
          res['Ny'] = arr[0]
          res['Nz'] = arr[1]
          res['Nstk'] = arr[2]

          bytes = f.read(4)                
          res['Nline'] = struct.unpack('i',bytes)[0]

          if res['Nline'] > 0:
            res['Lines'] = []
            for i in range(res['Nline']):
              res['Lines'].append({})

              bytes = f.read(8)   
              res['Lines'][i]['L'] = struct.unpack('d',bytes)[0]
              res['Lines'][i]['data'] = np.zeros((res['Nstk'], res['Ny'], res['Nz']),dtype=np.float64)
      
          bytes = f.read(4)
          res['NThom'] = struct.unpack('i',bytes)[0]


          if res['NThom'] > 0:
            res['Thoms'] = []
            for i in range(res['NThom']):
              res['Thoms'].append({})
              bytes = f.read(8)                
              res['Thoms'][i]['L'] = struct.unpack('d',bytes)[0]
              res['Thoms'][i]['data'] = np.zeros((2, res['Ny'], res['Nz']),dtype=np.float64)

                  
          bytes = f.read(8*4)
          arr2 = np.array(struct.unpack('d'*4,bytes))                
          res['FOV'] =  arr2.reshape((2,2))


              
          for iy in range(res['Ny']):
            for iz in range(res['Nz']):
              for i in range(res['Nline']):
                bytes = f.read(8*res['Nstk'])
                res['Lines'][i]['data'][:,iz,iy] = np.array(struct.unpack('d'*res['Nstk'],bytes))
              for i in range(res['NThom']):
                bytes = f.read(8*2)
                res['Thoms'][i]['data'][:,iz,iy] = np.array(struct.unpack('d'*2,bytes))  

          try:
            bytes = f.read(4)            
            res['nspec'] = int(struct.unpack('i',bytes)[0])

          except:
            res['nspec'] = 0


          if res['nspec'] > 0:
            try:
              bytes = f.read(4*2)
              arr = np.array(struct.unpack('i'*2,bytes))
              res['Nys'] = arr[0]
              res['Nzs'] = arr[1]
    
 
              bytes = f.read(8*4)
              arr = np.array(struct.unpack('d'*4,bytes))
                        
              res['FOVSPEC'] =  arr.reshape((2,2))
    
   
    
              res['prof'] = []
              for i in range(res['nspec']):
                res['prof'].append({})
                        
              for ispec in range(res['nspec']):
        
                bytes = f.read(4)            
                res['prof'][ispec]['nl'] = struct.unpack('i',bytes)[0]
                res['prof'][ispec]['L'] = np.zeros(res['prof'][ispec]['nl'],np.float64)
                res['prof'][ispec]['prof'] = np.zeros((res['prof'][ispec]['nl'], \
                    res['Nstk'], res['Nzs'], res['Nys']), np.float64)
   
                bytes = f.read(8*res['prof'][ispec]['nl']) 
                arr = np.array(struct.unpack('d'*res['prof'][ispec]['nl'],bytes))
                res['prof'][ispec]['L'][:] = arr[:]
    
              for iy in range(res['Nys']):
                for iz in range(res['Nzs']):
                  for ispec in range(res['nspec']):
                    for istk in range(res['Nstk']):
                      bytes = f.read(8*res['prof'][ispec]['nl'])   
                      arr = np.array(struct.unpack('d'*res['prof'][ispec]['nl'],bytes))    
                      res['prof'][ispec]['prof'][:,istk,iz,iy] = arr[:]
              
            except:
              print('error in reading profile ')




      
              
        case b'sigg':
      
          res['type'] = 'single grid'
      
          bytes = f.read(8*3)            
          arr = np.array(struct.unpack('d'*3,bytes))
          res['Y'] = arr[0]
          res['Z'] = arr[1]
          res['X'] = arr[2]
      
          print(res['Y'], res['Z'], res['X'])
      
          bytes = f.read(4)            
          res['Nstk'] = struct.unpack('i',bytes)[0]
          
          bytes = f.read(4)                
          res['Nline'] = struct.unpack('i',bytes)[0]
          print(' line number', res['Nline'])
          if res['Nline'] > 0:
            res['Lines'] = []
            for i in range(res['Nline']):
              res['Lines'].append({})

              bytes = f.read(8)   
              res['Lines'][i]['L'] = struct.unpack('d',bytes)[0]
              res['Lines'][i]['data'] = np.zeros(res['Nstk'],dtype=np.float64)
      
          bytes = f.read(4)
          res['NThom'] = struct.unpack('i',bytes)[0]
          print(' Thomson number', res['NThom'])

          if res['NThom'] > 0:
            res['Thoms'] = []
            for i in range(res['NThom']):
              res['Thoms'].append({})
              bytes = f.read(8)                
              res['Thoms'][i]['L'] = struct.unpack('d',bytes)[0]
              res['Thoms'][i]['data'] = np.zeros(2,dtype=np.float64)
              
      
      
          print(res['Nstk'], res['Nline'], res['NThom'])
      
          if res['Nline'] > 0 :
            for iline in range(res['Nline']):
              bytes = f.read(8*res['Nstk'])   
              arr = np.array(struct.unpack('d'*res['Nstk'],bytes))
              res['Lines'][iline]['data'][:] = arr[0:res['Nstk']]
      
          if res['NThom'] > 0 :    
            for iThom in range(res['NThom']):
              bytes = f.read(8*2) 
              arr = np.array(struct.unpack('d'*2,bytes))
              res['Thoms'][iline]['data'][:] = arr[0:2]



          try:
            bytes = f.read(4)            
            res['nspec'] = int(struct.unpack('i',bytes)[0])
            print('profile number', res['nspec'])
          except:
            res['nspec'] = 0
            print('no profile ')

          print('nspec',res['nspec'])
          if res['nspec'] > 0:
            try:

              res['prof'] = []          
              for i in range(res['nspec']):
                res['prof'].append({})
                        
              for ispec in range(res['nspec']):
                bytes = f.read(4)            
                res['prof'][ispec]['nl'] = struct.unpack('i',bytes)[0]
                res['prof'][ispec]['L'] = np.zeros(res['prof'][ispec]['nl'],np.float64)
                res['prof'][ispec]['prof'] = np.zeros((res['prof'][ispec]['nl'], \
                    res['Nstk']),np.float64)
                
                bytes = f.read(8*res['prof'][ispec]['nl']) 
                arr = np.array(struct.unpack('d'*res['prof'][ispec]['nl'],bytes))
                res['prof'][ispec]['L'][:] = arr[:]
                        
              for ispec in range(res['nspec']):
                for istk in range(res['Nstk']):
                  bytes = f.read(8*res['prof'][ispec]['nl'])   
                  arr = np.array(struct.unpack('d'*res['prof'][ispec]['nl'],bytes))    
                  res['prof'][ispec]['prof'][:,istk] = arr[:]
                
            except:
              print('error in reading profile ')



      
                            
        case b'sigp':
      
          res['type'] = 'single pixel'
      
          bytes = f.read(8*2)            
          arr = np.array(struct.unpack('d'*2,bytes))
          res['Y'] = arr[0]
          res['Z'] = arr[1]
      
          print(res['Y'], res['Z'])
      
          bytes = f.read(4)            
          res['Nstk'] = struct.unpack('i',bytes)[0]
          
          bytes = f.read(4)                
          res['Nline'] = struct.unpack('i',bytes)[0]
          print(' line number', res['Nline'])
          if res['Nline'] > 0:
            res['Lines'] = []
            for i in range(res['Nline']):
              res['Lines'].append({})

              bytes = f.read(8)   
              res['Lines'][i]['L'] = struct.unpack('d',bytes)[0]
              res['Lines'][i]['data'] = np.zeros(res['Nstk'],dtype=np.float64)
      
          bytes = f.read(4)
          res['NThom'] = struct.unpack('i',bytes)[0]
          print(' Thomson number', res['NThom'])

          if res['NThom'] > 0:
            res['Thoms'] = []
            for i in range(res['NThom']):
              res['Thoms'].append({})
              bytes = f.read(8)                
              res['Thoms'][i]['L'] = struct.unpack('d',bytes)[0]
              res['Thoms'][i]['data'] = np.zeros(2,dtype=np.float64)
              
      
      
          print(res['Nstk'], res['Nline'], res['NThom'])
      
          if res['Nline'] > 0 :
            for iline in range(res['Nline']):
              bytes = f.read(8*res['Nstk'])   
              arr = np.array(struct.unpack('d'*res['Nstk'],bytes))
              res['Lines'][iline]['data'][:] = arr[0:res['Nstk']]
      
          if res['NThom'] > 0 :    
            for iThom in range(res['NThom']):
              bytes = f.read(8*2) 
              arr = np.array(struct.unpack('d'*2,bytes))
              res['Thoms'][iline]['data'][:] = arr[0:2]



          try:
            bytes = f.read(4)            
            res['nspec'] = int(struct.unpack('i',bytes)[0])
            print('profile number', res['nspec'])
          except:
            res['nspec'] = 0
            print('no profile ')

          print('nspec',res['nspec'])
          if res['nspec'] > 0:
            try:

              res['prof'] = []          
              for i in range(res['nspec']):
                res['prof'].append({})
                        
              for ispec in range(res['nspec']):
                bytes = f.read(4)            
                res['prof'][ispec]['nl'] = struct.unpack('i',bytes)[0]
                res['prof'][ispec]['L'] = np.zeros(res['prof'][ispec]['nl'],np.float64)
                res['prof'][ispec]['prof'] = np.zeros((res['prof'][ispec]['nl'], \
                    res['Nstk']),np.float64)
                
                bytes = f.read(8*res['prof'][ispec]['nl']) 
                arr = np.array(struct.unpack('d'*res['prof'][ispec]['nl'],bytes))
                res['prof'][ispec]['L'][:] = arr[:]
                        
              for ispec in range(res['nspec']):
                for istk in range(res['Nstk']):
                  bytes = f.read(8*res['prof'][ispec]['nl'])   
                  arr = np.array(struct.unpack('d'*res['prof'][ispec]['nl'],bytes))    
                  res['prof'][ispec]['prof'][:,istk] = arr[:]
                
            except:
              print('error in reading profile ')


        case _ :
          print(2)

      
      f.close() 
      return res
    except:
        
      print("error in reading " + filename)
      f.close() 
      return {'check': False}
  
################################################################################ 

def interp3d(Data, i, j, k, ir, jr, kr):
    
    Ratio_XC = 1-ir
    Ratio_YC = 1-jr

    tmp1 = Ratio_XC*Data[i][j][k]+ir*Data[i+1][j][k]
    tmp2 = Ratio_XC*Data[i][j+1][k]+ir*Data[i+1][j+1][k]
    tmp3 = Ratio_YC*tmp1+jr*tmp2
    tmp1 = Ratio_XC*Data[i][j][k+1]+ir*Data[i+1][j][k+1]
    tmp2 = Ratio_XC*Data[i][j+1][k+1]+ir*Data[i+1][j+1][k+1]
    tmp4 = Ratio_YC*tmp1+jr*tmp2
    
    return (1-kr)*tmp3+kr*tmp4

################################################################################ 