#  Begona Ascaso. 07/02/17
#  Program that matches one structure to others using Cylindrical Matching and FoF+Cylindrical Matching. 
#  Also, if requested the program will compute completeness and purity rates.
#  Note that before producing the completeness and purity plots, we need to match two ways.

#  Usage:
#  python mainMCP.py MCP.pars
#  Modify input using the MCP.pars file

import numpy as N
import pylab as PL
import matchB as M
import CompPurB as CP
import functionsB as B
import sys

######### Input (read fom the file)


try:
    file_pars=sys.argv[1]#Read the input pars file
except:
#    print usage
    sys.exit()

expl,valor=N.loadtxt(file_pars,dtype='|S20,|S10',usecols=(0,1),unpack=True)

# Define the default values of the parameters

# Mean sigma of the survey
sigmaz=float(valor[0])

#Which limit in Halo mass your catalogues have?
minhalomass=valor[1]

# Are you mathing?
matching=valor[2]

# If so, with which method?
matchMethode=valor[3]

# If so, which direction?
matchDirection=valor[4]
#name=str(int(name))

# If so, do you want to print on the screen the tracking iterations?
printprocess=valor[5]
#name=str(int(name))

# Do you want to compute completeness and purity?
comppur=valor[6]
nbinsz=int(valor[7])
nbinsm=int(valor[8])


print 'MINHALOMASS= ',minhalomass
print 'MATCHING= ',matching
print 'METHODE= ',matchMethode
print 'DIRECTION= ',matchDirection

#sys.exit()

###############


if matching == 'yes':
    nameMatch = 'Match_'+matchMethode+'_'+matchDirection+'_'+minhalomass+'.txt'

if comppur == 'yes':
    nameComp = 'Completeness_'+matchMethode+'_'+minhalomass
    namePur = 'Purity_'+matchMethode+'_'+minhalomass

####################


#Match using cilindrical matching

if matching == 'yes' and matchMethode == 'CM':

    if matchDirection == 'H2D':

        idi,rai,deci,zi,ngali,richi,idf,raf,decf,zf,ngalf,richf,disti,rad=M.matchC('haloes_mock_'+minhalomass+'.txt','detection_mock.txt',sigmaz=sigmaz,printprocess=printprocess)

    if matchDirection == 'D2H':
        
        idi,rai,deci,zi,ngali,richi,idf,raf,decf,zf,ngalf,richf,disti,rad=M.matchC('detection_mock.txt','haloes_mock_'+minhalomass+'.txt',sigmaz=sigmaz,printprocess=printprocess)

       
    # Print the results into a file
    f = open(nameMatch,'w')
    for i in range(len(rai)):
	
    	columnnew="%.30s %.4f %.4f %.4f %.3f %.4f %.26s %.4f %.4f  %.4f %.3f  %.4f %.3f %.3f"%(idi[i],rai[i],deci[i],zi[i],ngali[i],richi[i],idf[i],raf[i],decf[i],zf[i],ngalf[i],richf[i],disti[i],rad[i])
    	f.write(columnnew+"\n")

    f.close()

#Match using FoF+cilindrical matching 

if matching == 'yes' and matchMethode == 'FCM':

    if matchDirection == 'H2D':
        
        idi,rai,deci,zi,ngali,richi,idf,raf,decf,zf,ngalf,richf,disti,rad=M.matchFC('haloes_mock_'+minhalomass+'.txt','detection_mock.txt',sigmaz=sigmaz,printprocess=printprocess)

    if matchDirection == 'D2H':
        
        idi,rai,deci,zi,ngali,richi,idf,raf,decf,zf,ngalf,richf,disti,rad=M.matchFC('detection_mock.txt','haloes_mock_'+minhalomass+'.txt',sigmaz=sigmaz,printprocess=printprocess)

    # Print the results into a file
    f = open(nameMatch,'w')
    for i in range(len(rai)):
	
    	columnnew="%.30s %.4f %.4f %.4f %.3f %.4f %.26s %.4f %.4f  %.4f %.3f  %.4f %.3f %.3f"%(idi[i],rai[i],deci[i],zi[i],ngali[i],richi[i],idf[i],raf[i],decf[i],zf[i],ngalf[i],richf[i],disti[i],rad[i])
    	f.write(columnnew+"\n")

    f.close()


# Compute completeness and purity from the previous matchings

if comppur == 'yes':

	##### Compute completeness
	
    	mc,zc,comp=CP.completeness('Match_'+matchMethode+'_H2D_'+minhalomass+'.txt',nbinsz=nbinsz,nbinsm=nbinsm)

	mmc=N.array(mc).reshape(nbinsm,nbinsz)
	zzc=N.array(zc).reshape(nbinsm,nbinsz)
	ccomp=N.array(comp).reshape(nbinsm,nbinsz)

	##### Plot completitud as a function of z for different mass bins
		
	zsc=zzc[1,:]	
	massc=mmc[:,1]

    	good=N.greater_equal(ccomp[0,:],0)
    	zsc2,ccomp2=B.multicompress(good,(zsc,ccomp[0,:]))
	PL.plot(zsc2,ccomp2,'ro-',label='$Rich$= %.2f'%massc[0])

    	good=N.greater_equal(ccomp[1,:],0)
    	zsc2,ccomp2=B.multicompress(good,(zsc,ccomp[1,:]))
	PL.plot(zsc2,ccomp2,'bo-',label='$Rich$= %.2f'%massc[1])

    	good=N.greater_equal(ccomp[2,:],0)
    	zsc2,ccomp2=B.multicompress(good,(zsc,ccomp[2,:]))
	PL.plot(zsc2,ccomp2,'go-',label='$Rich$= %.2f'%massc[2])

    	good=N.greater_equal(ccomp[3,:],0)
    	zsc2,ccomp2=B.multicompress(good,(zsc,ccomp[3,:]))
	PL.plot(zsc2,ccomp2,'co-',label='$Rich$= %.2f'%massc[3])

    	good=N.greater_equal(ccomp[4,:],0)
    	zsc2,ccomp2=B.multicompress(good,(zsc,ccomp[4,:]))
	PL.plot(zsc2,ccomp2,'mo-',label='$Rich$= %.2f'%massc[4])

    	good=N.greater_equal(ccomp[5,:],0)
    	zsc2,ccomp2=B.multicompress(good,(zsc,ccomp[5,:]))
	PL.plot(zsc2,ccomp2,'yo-',label='$Rich$= %.2f'%massc[5])
	
	PL.xlabel('z')
	PL.ylabel('Comp')
	PL.ylim(-0.1,1.2)
	PL.legend(loc='upper right',fontsize='small') #x-small
#	PL.show()
	PL.savefig(nameComp+'.png')
	PL.close()
	
	##### Print the results into a file
	
	f = open(nameComp+'.txt','w')	
	f.write("# z  and Mass bins: ")

	# Transform the mass into a row array

	massc=N.array(massc).reshape(1,nbinsm) 
	N.savetxt(f,massc,fmt='%.4f')
	f.write("\n")

	# Transform the zc into a column array and append it to the completitud
	
	zsc=N.array(zsc).reshape(nbinsz,1) 
 	cczz=N.append(zsc,ccomp.transpose(),axis=1)
	
	N.savetxt(f,cczz,fmt='%.4f')

	f.close()
        
	

	##### Compute purity
	 
    	rp,zp,pur=CP.purity('Match_'+matchMethode+'_D2H_'+minhalomass+'.txt',nbinsz=nbinsz,nbinsm=nbinsm)

	rrp=N.array(rp).reshape(nbinsm,nbinsz)
	zzp=N.array(zp).reshape(nbinsm,nbinsz)
	ppur=N.array(pur).reshape(nbinsm,nbinsz)

	##### Plot purity as a function of z for different richness bins
		
	zsp=zzp[1,:]	
	richp=rrp[:,1]

    	good=N.greater_equal(ppur[0,:],0)
    	zsp2,ppur2=B.multicompress(good,(zsp,ppur[0,:]))
	PL.plot(zsp2,ppur2,'ro-',label='$Rich$= %.2f'%richp[0])

    	good=N.greater_equal(ppur[1,:],0)
    	zsp2,ppur2=B.multicompress(good,(zsp,ppur[1,:]))
	PL.plot(zsp2,ppur2,'bo-',label='$Rich$= %.2f'%richp[1])

    	good=N.greater_equal(ppur[2,:],0)
    	zsp2,ppur2=B.multicompress(good,(zsp,ppur[2,:]))
	PL.plot(zsp2,ppur2,'go-',label='$Rich$= %.2f'%richp[2])

    	good=N.greater_equal(ppur[3,:],0)
    	zsp2,ppur2=B.multicompress(good,(zsp,ppur[3,:]))
	PL.plot(zsp2,ppur2,'co-',label='$Rich$= %.2f'%richp[3])

    	good=N.greater_equal(ppur[4,:],0)
    	zsp2,ppur2=B.multicompress(good,(zsp,ppur[4,:]))
	PL.plot(zsp2,ppur2,'mo-',label='$Rich$= %.2f'%richp[4])

    	good=N.greater_equal(ppur[5,:],0)
    	zsp2,ppur2=B.multicompress(good,(zsp,ppur[5,:]))
	PL.plot(zsp2,ppur2,'yo-',label='$Rich$= %.2f'%richp[5])

	PL.xlabel('z')
	PL.ylabel('Purity')
	PL.ylim(-0.1,1.2)
	PL.legend(loc='upper right',fontsize='small') #x-small
#	PL.show()
	PL.savefig(namePur+'.png')
	PL.close()
	
	##### Print the results into a file
	
	f = open(namePur+'.txt','w')	
	f.write("# z  and Richness bins: ")

	# Transform the mass into a row array

	rrp=N.array(richp).reshape(1,nbinsm) 
	N.savetxt(f,rrp,fmt='%.4f')
	f.write("\n")

	# Transform the zc into a column array and append it to the completitud
	
	zsp=N.array(zsp).reshape(nbinsz,1) 
 	ppzz=N.append(zsp,ppur.transpose(),axis=1)
	
	N.savetxt(f,ppzz,fmt='%.4f')

	f.close()
