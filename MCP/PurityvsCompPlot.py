#  Begona Ascaso. 07/02/17
#  Program that provided the Matching both ways performed with either CM or FCM, do:
# 1. computes completeness and purity,
# 2, roughly fits a log-log linear relation between mass and richness to
# 3, plot purity vs completeness


#  Usage:
#  python PurityvsCompPlot.py
#  Modify matchMethode to produce this plot for different matching methodes

import numpy as N
import pylab as PL
import CompPurB as CP
import functionsB as B


#matchMethode = 'CM'
matchMethode = 'FCM'

cutHalo = 'no'
cutDetect = 'yes'


nbinsz = 4 #9
binz=2.0/nbinsz


# Purity vs Completeness varying the cut in the Halo Catalogues

if cutHalo == 'yes':
    
    pur113,comp113,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_113_1snr.txt','Match_'+matchMethode+'_H2D_113_1snr.txt',nbinsz=nbinsz,aa0=-99,aa1=-99)
#Keep the fit of the Mass-Observable relation the one obtained from the whole sample without restrictions.
    a0Good=a0
    a1Good=a1

    pur313,comp313,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_313_1snr.txt','Match_'+matchMethode+'_H2D_313_1snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)
    pur513,comp513,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_513_1snr.txt','Match_'+matchMethode+'_H2D_513_1snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)
    pur713,comp713,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_713_1snr.txt','Match_'+matchMethode+'_H2D_713_1snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)
    pur214,comp214,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_214_1snr.txt','Match_'+matchMethode+'_H2D_214_1snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)
    pur114,comp114,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_114_1snr.txt','Match_'+matchMethode+'_H2D_114_1snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)


    subpl='22' #nbinsz x nbinsz
    zname ='<z<'
    for i in range(nbinsz): 
        PL.subplot(subpl+`i+1`)
        for j in range(6):
            if j==0: namepur=pur113; namecomp=comp113
            if j==1: namepur=pur313; namecomp=comp313
            if j==2: namepur=pur513; namecomp=comp513
            if j==3: namepur=pur713; namecomp=comp713
            if j==4: namepur=pur114; namecomp=comp114
            if j==5: namepur=pur214; namecomp=comp214
    
            good=N.greater_equal(namepur[i,:],0)*N.greater_equal(namecomp[i,:],0)
            ppp,ccc=B.multicompress(good,(namepur[i,:],namecomp[i,:]))

            if j==0:
                PL.plot(ppp,ccc,'ro-',label='$M_{lim}$= 1.0e13')
                PL.ylim((-0.15,1.15))
                PL.xlim((-0.15,1.15))
                PL.xlabel('Purity', fontsize='small')
                PL.ylabel('Completeness', fontsize='small')
            if j==1: PL.plot(ppp,ccc,'bo-',label='$M_{lim}$= 3.0e13')
            if j==2: PL.plot(ppp,ccc,'go-',label='$M_{lim}$= 5.0e13')
            if j==3: PL.plot(ppp,ccc,'co-',label='$M_{lim}$= 7.0e13')
            if j==4: PL.plot(ppp,ccc,'mo-',label='$M_{lim}$= 1.0e14')
            if j==5: PL.plot(ppp,ccc,'yo-',label='$M_{lim}$= 2.0e14')

        PL.title(`binz*i`+zname+`binz*(i+1)`, fontsize='small')        

    PL.legend(loc='lower right',fontsize='x-small') #small
    #PL.show()
    PL.savefig('PurvsComp_'+matchMethode+'_snrlim1.png')


# Purity vs Completeness varying the cut in the Halo Catalogues

if cutDetect == 'yes':
    
    pur113,comp113,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_113_1snr.txt','Match_'+matchMethode+'_H2D_113_1snr.txt',nbinsz=nbinsz,aa0=-99,aa1=-99)
#Keep the fit of the Mass-Observable relation the one obtained from the whole sample without restrictions.
    a0Good=a0
    a1Good=a1

    pur313,comp313,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_113_3snr.txt','Match_'+matchMethode+'_H2D_113_3snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)
    pur513,comp513,a0,a1=CP.purityvscomp('Match_'+matchMethode+'_D2H_113_5snr.txt','Match_'+matchMethode+'_H2D_113_5snr.txt',nbinsz=nbinsz,aa0=a0Good,aa1=a1Good)


    subpl='22' #nbinsz x nbinsz
    zname ='<z<'
    for i in range(nbinsz): 
        PL.subplot(subpl+`i+1`)
        for j in range(3):
            if j==0: namepur=pur113; namecomp=comp113
            if j==1: namepur=pur313; namecomp=comp313
            if j==2: namepur=pur513; namecomp=comp513
    
            good=N.greater_equal(namepur[i,:],0)*N.greater_equal(namecomp[i,:],0)
            ppp,ccc=B.multicompress(good,(namepur[i,:],namecomp[i,:]))

            if j==0:
                PL.plot(ppp,ccc,'ro-',label='$SNR_{lim}$= 1')
                PL.ylim((-0.15,1.15))
                PL.xlim((-0.15,1.15))
                PL.xlabel('Purity', fontsize='small')
                PL.ylabel('Completeness', fontsize='small')
            if j==1: PL.plot(ppp,ccc,'bo-',label='$SNR_{lim}$= 3')
            if j==2: PL.plot(ppp,ccc,'go-',label='$SNR_{lim}$= 5')

        PL.title(`binz*i`+zname+`binz*(i+1)`, fontsize='small')        

    PL.legend(loc='lower right',fontsize='x-small') #small
    #PL.show()
    PL.savefig('PurvsComp_'+matchMethode+'_mlim113.png')

