# Begona Ascaso

import os,sys
import numpy as N
import string as S
import astropy as A
from astropy import coordinates as C
from astropy import units as U
from astropy.cosmology import WMAP9 as cosmo
import functionsB as B


# Match using cilindrical matching

def matchC(cat1, cat2,sigmaz,printprocess):#, dz, opts):

# Read cat1 #rich1 = Mhalo
    id1,ra1,dec1,z1,ngal1,rich1,r2001=N.loadtxt(cat1,dtype='|S20,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6),unpack=True)

# Sort the catalogue putting the most massive first
    be=(-rich1).ravel().argsort()
    id1=id1[be] ; ra1=ra1[be] ; dec1=dec1[be] ; z1=z1[be]; ngal1=ngal1[be]; rich1=rich1[be]; r2001=r2001[be]

# Read cat2 # rich2 = SNR
    id2,ra2,dec2,z2,ngal2,rich2,r2002=N.loadtxt(cat2,dtype='|S20,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6),unpack=True)
  
# Sort the catalogue putting the firts ranked first
    be=(-rich2).ravel().argsort()
    id2=id2[be] ; ra2=ra2[be] ; dec2=dec2[be] ; z2=z2[be]; ngal2=ngal2[be]; rich2=rich2[be] ; r2002=r2002[be]
    
    c2 = C.SkyCoord(ra=ra2*U.degree,dec=dec2*U.degree)

    origID = []
    origRA = []
    origDEC = []
    origZ = []
    origNGAL = []
    origRICH = []

    matchID = []
    matchRA = []
    matchDEC = []
    matchZ = []
    matchNGAL = []
    matchRICH = []
    
    matchDISTI = []
    origRAD = []

    kpcmin1=cosmo.kpc_proper_per_arcmin(z1)
    r2001Mpc=kpcmin1.value*r2001/1000.
    
    null = float(-99.)

    for i in range(len(id1)):
#    for i in range(84100,84200):
#    for i in range(0,299):
    	
	if printprocess=='yes': 
         if i/100. == int(i/100.): 
             print i 
   
	c1 = C.SkyCoord(ra=ra1[i]*U.degree,dec=dec1[i]*U.degree)
	sep = c1.separation(c2) #arcmins
	sepV=sep.arcminute

#	ja=N.less_equal((N.array(sepV)-r2001[i]),0.0)*N.less_equal(abs(z2-z1[i])-2.*sigmaz*(1.+z1[i]),0.)
	ja=N.less_equal((N.array(sepV)-2.*r2001[i]),0.0)*N.less_equal(abs(z2-z1[i])-2.*sigmaz*(1.+z1[i]),0.) #2 veces radio 2R
#	ja=N.less_equal((N.array(sepV)-2.*r2001[i]),0.0) #2 veces radio 2R
	
    	id2S,ra2S,dec2S,z2S,ngal2S,rich2S,sepS=B.multicompress(ja,(id2,ra2,dec2,z2,ngal2,rich2,sepV))

	origID.append(id1[i])
	origRA.append(ra1[i])
	origDEC.append(dec1[i])
	origZ.append(z1[i])
	origNGAL.append(ngal1[i])
	origRICH.append(rich1[i])
	origRAD.append(r2001Mpc[i])

# picks up the most massive (the files are originally sort)
	if len(ra2S) >0 :

		matchID.append(id2S[0])
		matchRA.append(ra2S[0])
		matchDEC.append(dec2S[0])
		matchZ.append(z2S[0])
		matchNGAL.append(ngal2S[0])
		matchRICH.append(rich2S[0])

		kpcmin2=cosmo.kpc_proper_per_arcmin(z2S[0])
		distiMpc=kpcmin2.value*sepS[0]/1000.
		
		matchDISTI.append(distiMpc)

# I save all even if not matched (necessary for the Completeness & Purity computations)

	if len(ra2S) ==0 :

		matchID.append(null)
		matchRA.append(null)
		matchDEC.append(null)
		matchZ.append(null)
		matchNGAL.append(null)
		matchRICH.append(null)		
		matchDISTI.append(null)
    

# Check how many are matched
    ja=N.greater_equal((N.array(matchDISTI)),0.0)
    la=N.compress(ja,N.array(matchDISTI))
       
    print len(la), 'matched out of ', len(ra1)
    
    return origID,origRA,origDEC,origZ,origNGAL,origRICH,matchID,matchRA,matchDEC,matchZ,matchNGAL,matchRICH,matchDISTI,origRAD



# Match using FoF+cilindrical matching

def matchFC(cat1, cat2,sigmaz,printprocess):#, dz, opts):

# Read cat1 #rich1 = Mhalo
    id1,ra1,dec1,z1,ngal1,rich1,r2001=N.loadtxt(cat1,dtype='|S20,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6),unpack=True)


# Sort the catalogue putting the most massive first
    be=(-rich1).ravel().argsort()
    id1=id1[be] ; ra1=ra1[be] ; dec1=dec1[be] ; z1=z1[be]; ngal1=ngal1[be]; rich1=rich1[be]; r2001=r2001[be]

# Read cat2 # rich2 = SNR
    id2,ra2,dec2,z2,ngal2,rich2,r2002=N.loadtxt(cat2,dtype='|S20,float,float,float,float,float,float',usecols=(0,1,2,3,4,5,6),unpack=True)
      

    kpcmin=cosmo.kpc_proper_per_arcmin(z2)
    r2002=(z2*0.+1.0)*1000./kpcmin.value #1 Mpc in minutes
    
# Sort the catalogue putting the firts ranked first
    be=(-rich2).ravel().argsort()
    id2=id2[be] ; ra2=ra2[be] ; dec2=dec2[be] ; z2=z2[be]; ngal2=ngal2[be]; rich2=rich2[be] ; r2002=r2002[be]
    
    c2 = C.SkyCoord(ra=ra2*U.degree,dec=dec2*U.degree)


    origID = []
    origRA = []
    origDEC = []
    origZ = []
    origNGAL = []
    origRICH = []
    origRAD = []

    matchID = []
    matchRA = []
    matchDEC = []
    matchZ = []
    matchNGAL = []
    matchRICH = []
    
    matchDISTI = []

    kpcmin1=cosmo.kpc_proper_per_arcmin(z1)
    r2001Mpc=kpcmin1.value*r2001/1000.

    null = float(-99.)

    for i in range(len(id1)):
#    for i in range(50000,50100):
    	
	if printprocess=='yes': 
         if i/100. == int(i/100.): 
             print i 
    
	c1 = C.SkyCoord(ra=ra1[i]*U.degree,dec=dec1[i]*U.degree)
	sep = c1.separation(c2) #arcmins
	sepV=sep.arcminute

#	ja=N.less_equal((N.array(sepV)-r2001[i]),0.0)*N.less_equal(abs(z2-z1[i])-2.*sigmaz*(1.+z1[i]),0.)
	ja=N.less_equal((N.array(sepV)-2.*r2001[i]),0.0)*N.less_equal(abs(z2-z1[i])-2.*sigmaz*(1.+z1[i]),0.) #2 veces radio 2R
	
    	id2S,ra2S,dec2S,z2S,ngal2S,rich2S,r2002S,sepS=B.multicompress(ja,(id2,ra2,dec2,z2,ngal2,rich2,r2002,sepV))

	origID.append(id1[i])
	origRA.append(ra1[i])
	origDEC.append(dec1[i])
	origZ.append(z1[i])
	origNGAL.append(ngal1[i])
	origRICH.append(rich1[i])
	origRAD.append(r2001Mpc[i])

	if len(ra2S) >0 :

    		c2F = C.SkyCoord(ra=ra2S*U.degree,dec=dec2S*U.degree)

		# I search for Friends of Friends
		for j in range(len(id2S)):

    			matchFID = []
    			matchFRA = []
    			matchFDEC = []
    			matchFZ = []
    			matchFNGAL = []
    			matchFRICH = []
    
    			matchFDISTI = []
	
			c1F = C.SkyCoord(ra=ra2S[j]*U.degree,dec=dec2S[j]*U.degree)
			sepF = c1F.separation(c2F) #arcmins
			sepVF=sepF.arcminute
            	ja=N.less_equal((N.array(sepVF)-2.*r2002S[j]),0.0)*N.less_equal(abs(z2S-z2S[j])-2.*sigmaz*(1.+z2S[j]),0.) #2 veces radio 2R
                id2F,ra2F,dec2F,z2F,ngal2F,rich2F,r2002F,sepF=B.multicompress(ja,(id2S,ra2S,dec2S,z2S,ngal2S,rich2S,r2002S,sepS))
                
                if len(id2F)>0:
                    matchFID.append(id2F)	
                    matchFRA.append(ra2F)	
                    matchFDEC.append(dec2F)	
                    matchFZ.append(z2F)	
                    matchFNGAL.append(ngal2F)	
                    matchFRICH.append(rich2F)	

			# Append the distances of the FoF
                    kpcmin2F=cosmo.kpc_proper_per_arcmin(z2F)
                    distiMpcF=kpcmin2F.value*sepF/1000.
                    matchFDISTI.append(distiMpcF)	
                
#                print j,len(id2F),len(matchFID),r2002S[j]

    		# Sort the catalogue putting the most massive first and select the most massive
		if len(matchFID)>0: 
#                    	print 'ea',matchFID
#                     	print 'yo',ra1[i],dec1[i],z1[i]
#                     	print 'ea',matchFDISTI,matchFRA,matchFDEC
            		matchFID=N.concatenate(matchFID,axis=0); matchFRA=N.concatenate(matchFRA,axis=0) 
            		matchFDEC=N.concatenate(matchFDEC,axis=0); matchFZ=N.concatenate(matchFZ,axis=0) 
            		matchFNGAL=N.concatenate(matchFNGAL,axis=0); matchFRICH=N.concatenate(matchFRICH,axis=0) 
            		matchFDISTI=N.concatenate(matchFDISTI,axis=0)
#            		print 'ea2',matchFDISTI,matchFRA,matchFDEC
            
            		be=(-matchFRICH).ravel().argsort(); 
#            		be=(matchFDISTI).ravel().argsort();  
                      	matchFID=matchFID[be[0]] ; matchFRA=matchFRA[be[0]] ; matchFDEC=matchFDEC[be[0]] ; matchFZ=matchFZ[be[0]]
            		matchFNGAL=matchFNGAL[be[0]]; matchFRICH=matchFRICH[be[0]]; matchFDISTI=matchFDISTI[be[0]]
#                      	print matchFDISTI

            
            		matchID.append(matchFID)	
            		matchRA.append(matchFRA)	
            		matchDEC.append(matchFDEC)	
            		matchZ.append(matchFZ)	
            		matchNGAL.append(matchFNGAL)	
            		matchRICH.append(matchFRICH)	
            		matchDISTI.append(matchFDISTI)	

        	if len(matchFID) ==0:
            		matchID.append(null)
            		matchRA.append(null)
            		matchDEC.append(null)
            		matchZ.append(null)
            		matchNGAL.append(null)
            		matchRICH.append(null)		
            		matchDISTI.append(null)
		
# I save all even if not matched (necessary for the Completeness & Purity computations)

	if len(ra2S) ==0 :

		matchID.append(null)
		matchRA.append(null)
		matchDEC.append(null)
		matchZ.append(null)
		matchNGAL.append(null)
		matchRICH.append(null)		
		matchDISTI.append(null)

# Check how many are matched
    ja=N.greater_equal((N.array(matchDISTI)),0.0)
    la=N.compress(ja,N.array(matchDISTI))
       
    print len(la), 'matched out of ', len(ra1)
    
    return origID,origRA,origDEC,origZ,origNGAL,origRICH,matchID,matchRA,matchDEC,matchZ,matchNGAL,matchRICH,matchDISTI,origRAD

