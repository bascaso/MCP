# Begona Ascaso

import numpy as N
import pylab as PL
import sys

# Compute completeness from a catalogue searching for the counterparts of haloes

################## Need to check that where is doing the correct thing.

def completeness(cat, nbinsz, nbinsm):

# Read cat1 #rich1 = Mhalo
    id1,ra1,dec1,z1,ngal1,rich1,id2,ra2,dec2,z2,ngal2,rich2,disti,rad1=N.loadtxt(cat,dtype='|S20,float,float,float,float,float,|S20,float,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)
   
    binsz=(max(z1)-min(z1))*1./(nbinsz-1.)
    binsm=(max(rich1)-min(rich1))*1./(nbinsm-1.)

    mm = []
    zz = []
    cc = []
    
    for i in range(0,nbinsm):

	meanm=min(rich1)+(i+1./2)*binsm
	    
	minrichbin=min(rich1)+i*binsm
	maxrichbin=min(rich1)+(i+1.)*binsm  
    	
	for j in range(0,nbinsz):

		meanz=min(z1)+(j+1./2)*binsz

		minzbin=min(z1)+j*binsz
		maxzbin=min(z1)+(j+1.)*binsz  
  

    		indtotal=N.where((rich1 >= minrichbin) & (rich1 <= maxrichbin) & (z1 >= minzbin) & (z1 <= maxzbin))
		ntotal=len(z1[indtotal])

		radm = N.median(rad1[indtotal])
				
    		indmatched=N.where((disti >=0) & (disti <=2.0*radm) & (rich1 >= minrichbin) & (rich1 <= maxrichbin) & (z1 >= minzbin) & (z1 <= maxzbin))
		nmatched=len(z1[indmatched])


		if ntotal==0:
			comp=-99 
		else:
			comp = nmatched*1.0/ntotal
		
		mm.append(meanm)
		zz.append(meanz)
		cc.append(comp)
     
    return mm,zz,cc
     

# Compute purity from a catalogue searching for the counterparts of detections

def purity(cat, nbinsz, nbinsm):

# Read cat1 
    id1,ra1,dec1,z1,ngal1,rich1,id2,ra2,dec2,z2,ngal2,rich2,disti,rad1=N.loadtxt(cat,dtype='|S20,float,float,float,float,float,|S20,float,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)

    binsz=(max(z1)-min(z1))/(nbinsz-1)
    binsm=(max(rich1)-min(rich1))/(nbinsm-1)
 
    
    mm = []
    zz = []
    pp = []
    
    
    for i in range(0,nbinsm):

	meanm=min(rich1)+(i+1./2)*binsm
		    
	minrichbin=min(rich1)+i*binsm
	maxrichbin=min(rich1)+(i+1.)*binsm  
    	
	for j in range(0,nbinsz):

		meanz=min(z1)+(j+1./2)*binsz

		minzbin=min(z1)+j*binsz
		maxzbin=min(z1)+(j+1.)*binsz  

    		indtotal=N.where((rich1 >= minrichbin) & (rich1 <= maxrichbin) & (z1 >= minzbin) & (z1 <= maxzbin))
		ntotal=len(z1[indtotal])

		radm = N.median(rad1[indtotal])
				
    		indmatched=N.where((disti >=0) & (disti <=2.0*radm) & (rich1 >= minrichbin) & (rich1 <= maxrichbin) & (z1 >= minzbin) & (z1 <= maxzbin))
		nmatched=len(z1[indmatched])

		if ntotal==0:
			puri=-99. 
		else:
			puri = nmatched*1.0/ntotal
		
		mm.append(meanm)
		zz.append(meanz)
		pp.append(puri)
		
    return mm,zz,pp
     

# Compute Purity vs Completeness
          
def purityvscomp(catp, catc, nbinsz,aa0,aa1):

# Read cat1 
    id1p,ra1p,dec1p,z1p,ngal1p,rich1p,id2p,ra2p,dec2p,z2p,ngal2p,rich2p,distip,rad1p=N.loadtxt(catp,dtype='|S20,float,float,float,float,float,|S20,float,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)

# Only fit if we do not pass them the fit (in order to fit the best data set)
    if aa0<-50 and aa1<-50:
        good=N.where(distip >=0)
#        PL.plot(rich1p[good],rich2p[good],'ro')    
#       This can be of course do make better.
        aa = N.polyfit(rich1p[good],rich2p[good], 1)
        aa0=aa[1]; aa1=aa[0]
#        PL.plot(rich1p[good],aa0+aa1*rich1p[good],'bo-')
#        PL.show()

    mass1p = aa0+aa1*rich1p
    
    id1,ra1,dec1c,z1c,ngal1c,rich1c,id2c,ra2c,dec2c,z2c,ngal2c,rich2c,distic,rad1c=N.loadtxt(catc,dtype='|S20,float,float,float,float,float,|S20,float,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6,7,8,9,10,11,12,13),unpack=True)

    mass1c = rich1c

#    binsz=(max(z1p)-min(z1p))/(nbinsz-1)
    binsz=(2.0)/(nbinsz-1)
     
    mm = [N.log10(1.0e13), N.log10(3.0e13), N.log10(5.0e13), N.log10(7.0e13), N.log10(8.0e13), N.log10(1.0e14), N.log10(2.0e14), N.log10(1.0e17)]

    pur = N.zeros((nbinsz,len(mm)-1))
    com = N.zeros((nbinsz,len(mm)-1))
    
    for j in range(0,nbinsz):


	minzbin=j*binsz
	maxzbin=(j+1.)*binsz  
     
     	pp = []
     	cc = []
     
     	for i in range(len(mm)-1):

# Compute purity
        	indtotal=N.where((mass1p >= mm[i]) & (mass1p <= mm[i+1]) & (z1p >= minzbin) & (z1p <= maxzbin))
           	ntotal=len(z1p[indtotal])

		radm = N.median(rad1p[indtotal])
				
    		indmatched=N.where((distip >=0) & (distip <=2.0*radm) & (mass1p >= mm[i]) & (mass1p <= mm[i+1]) & (z1p >= minzbin) & (z1p <= maxzbin))
		nmatched=len(z1p[indmatched])
           
		if ntotal==0:
			puri=-99
		else:
			puri = nmatched*1.0/ntotal
           
		pp.append(puri)
		
# Compute completeness
        	indtotal=N.where((mass1c >= mm[i]) & (mass1c <= mm[i+1]) & (z1c >= minzbin) & (z1c <= maxzbin))
           	ntotal=len(z1c[indtotal])

		radm = N.median(rad1c[indtotal])
				
    		indmatched=N.where((distic >=0) & (distic <=2.0*radm) & (mass1c >= mm[i]) & (mass1c <= mm[i+1]) & (z1c >= minzbin) & (z1c <= maxzbin))
		nmatched=len(z1c[indmatched])


		if ntotal==0:
			comp=-99
		else:
			comp = nmatched*1.0/ntotal

		cc.append(comp)
  
    	pur[j,:] = pp  
    	com[j,:] = cc  
  
		
    return pur,com,aa0,aa1
     
          
