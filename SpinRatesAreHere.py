#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 10:20:14 2020

@author: jroquette

This is a package I built to organize and load databases of spin rates for 
low mass stars.


The package currently includes the following regions:
    - hPer
    - NGC2264
    - NGC2362
    - USco
    - NGC 6530

Last Update: 17 August 2020
"""
import numpy as np
import astropy.coordinates as ac
from astropy.table import Table
import MacOSFile

print('Version 0 - LAST UPDATE 16 June 2020')
print('Package currently contains: hPer, NGC2264, NGC2362, USco, NGC 6530')
print('Disk is >0 if star is known to have a disk, <0 if the star is is known diskless and =0 if unkown')
class hPer:
    """
    hPer database 
    
    Reference table is Moraux+2013
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/560/A13    
    Reference for Spectral types is Currie et al. 2010
    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJS/186/191
    Reference for Mass transformations: ask Julia Roquette
    Usage:
        hPer=SpinRatesAreHere.hPer()
        
        mass_type = 0 source paper mass
                    1 MESA with median extinction
                    2 Baraffe+98 with median extinction
    LAST UPDATE: 14 June 2020
    """
    def __init__(self,filename='hPer_P_m_updated_14_Jun_2020.fit',datadir='tables/',mass_type=0):
        print('Loading data,RA,Dec,Prot,Mass,Amp,SpT,EBV')
        print('This table had binary stars removed following the flags in Moraux et al. 2013')
        self.data=Table.read(datadir+filename, format='fits') 
        self.RA=self.data['RAJ2000']
        self.Dec=self.data['DEJ2000']
        self.Prot=self.data['Per']
        if mass_type==0:
            print('Using Masses from Moraux et al. 2013')            
            self.Mass=self.data['Mass'] 
        elif mass_type==1:
            print('Using MESA masses with median extinction')
            self.Mass=self.data['Mass_MESA_median'] 
        elif mass_type==2:
            print('Using B98 masses with median extinction')
            self.Mass=self.data['Mass_B98_median'] 
        self.Mass[self.Mass<=0.]=np.nan    
        self.Amp=self.data['Amp'] 
        self.SpT=self.data['SpectralType']
        self.EBV=self.data['E_B-V_']
        self.cat= ac.SkyCoord(self.RA,self.Dec, unit="deg")  
        print('PS1 for PanSTARRs DR2 stacked data')
        print('TwoMASS for 2MASS data')     
        print('OriginalPhotometry for VIc, and i_cfht data')
        print('ClusterInfo forgeneral info about the cluster')
    def OriginalPhotometry(self):
        self.V=self.data['Vmag']
        self.Ic=self.data['Icmag']
        self.i_mag=self.data['i_mag']
    def TwoMASS(self):        
        self.Jmag=self.data['Jmag']
        self.Hmag=self.data['Hmag']   
        self.Ksmag=self.data['Kmag']
    def PS1(self):
        self.g=self.data['gPSFMag']
        self.g_e=self.data['gPSFMagErr']
        self.r=self.data['rPSFMag']
        self.r_e=self.data['rPSFMagErr']
        self.i=self.data['iPSFMag']
        self.i_e=self.data['iPSFMagErr']
        self.z=self.data['zPSFMag']
        self.z_e=self.data['zPSFMagErr']
        self.y=self.data['yPSFMag']
        self.y_e=self.data['yPSFMagErr']  
    def ClusterInfo(self,datadir='tables/'):
        cluster_info=MacOSFile.pickle_load(datadir+'hPer_ClusterInfo.npy')
        print('Distance from Gaia DR2 -',cluster_info['Gaia_distance_ref'])
        self.dist=cluster_info['Gaia_distance']
        self.DM=5*np.log10(self.dist/10)
        print(self.dist,' pc, DM= ',self.DM)        
        print('Median E(B-V) from Currie et al. 2010')
        print(cluster_info['comment_AV'])
        self.EBV=cluster_info['EBV']
        print(self.EBV)
        self.Rv=3.1
        self.Av=self.Rv*self.EBV
        print('Rv=',self.Rv,' Av=',self.Av)
        print('Age from CMoraux et al. 2013')
        self.Age=cluster_info['Age']
        print(cluster_info['FeH_ref'])
        self.FeH=cluster_info['FeH']
        print(self.FeH)
        
class NGC2264:
    """
    Reference table is a combination of:
        1. Venuti et al. 2017: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/599/A23
        2. Affer et al. 2013: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/430/1433
        3. Cieza et al. 2007: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJ/671/605
        4. Lamm et al. 2005: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/430/1005
        5. Makidon et al. 2004: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/127/2228
    Reference for disks: Alana Souza (Private comunication)
    References for Mass: Ask Julia Roquette, but it is complicated. 
    
    
    mass_type = 0 source paper mass
                1 MESA with mixed - comes either from Teff-Mass or from J or iband
                2 Baraffe+98 wmixed - comes either from Teff-Mass or from J or iband
    LAST UPDATE: 14 June 2020
                
    """
    def __init__(self,filename='NGC2264_P_M_updated_14_Jun_2020.fit',datadir='tables/',mass_type=0):
        print('Reference table is 1. Venuti et al. 2017, or 2. Affer et al. 2013, or 3. Cieza et al. 2007, or 4. Lamm et al., or 5. Makidon et al. 2004')
        self.data=Table.read(datadir+filename, format='fits') 
        print('Loading data,RA,Dec,Prot,Mass,Amp,SpT,Av,Disk')
        self.RA=self.data['RA_1']
        self.Dec=self.data['DE']
        self.Prot=self.data['Prot']
        if mass_type==0:
            print('Masses from Venuti et al. 2017')
            self.Mass=self.data['M_V14']
        elif mass_type==1:
            self.Mass=self.data['Mass_MESA_mixed']
        elif mass_type==2:
            self.Mass=self.data['Mass_B98_mixed']         
        self.Mass[self.Mass<=0.]=np.nan 
        self.AV=self.data['Av_V14']
        self.SpT=self.data['spt_V14']    
        self.Disk=self.data['Disked']
        print('PS1 for PanSTARRs DR2 stacked data')
        print('TwoMASS for 2MASS data')       
        print('ClusterInfo forgeneral info about the cluster')
        
    def PS1(self):
        self.g=self.data['gPSFMag']
        self.g_e=self.data['gPSFMagErr']
        self.r=self.data['rPSFMag']
        self.r_e=self.data['rPSFMagErr']
        self.i=self.data['iPSFMag']
        self.i_e=self.data['iPSFMagErr']
        self.z=self.data['zPSFMag']
        self.z_e=self.data['zPSFMagErr']
        self.y=self.data['yPSFMag']
        self.y_e=self.data['yPSFMagErr']  
    def TwoMASS(self):            
        self.Jmag=self.data['j_m']
        self.Jmag_e=self.data['j_cmsig']
        self.Hmag=self.data['h_m']
        self.Hmag_e=self.data['h_cmsig']
        self.Ksmag=self.data['k_m']
        self.Ksmag_e=self.data['k_cmsig']
    def ClusterInfo(self,datadir='tables/'):
        cluster_info=MacOSFile.pickle_load(datadir+'NGC2264_ClusterInfo.npy')
        print('Distance from Gaia DR2 -',cluster_info['Gaia_distance_ref'])
        self.dist=cluster_info['Gaia_distance']
        self.DM=5.*np.log10(self.dist/10.)
        print(cluster_info['comment_AV'])
        self.Av=cluster_info['AV']
        print(self.Av)
        self.EBV=self.Av/3.1
        print('Age adopted by Venuti+2017')
        self.Age=cluster_info['Age']
        print(cluster_info['FeH_ref'])
        self.FeH=cluster_info['FeH']
class USco:
    """
    Reference table:
        Rebull et al. 2018  http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/AJ/155/196

    Reference for disks: Rebull et al. 2018
    References for Mass: Ask Julia Roquette, but it is complicated. 
    
    
    mass_type = 0 source paper mass
                1 MESA with individual extinction: from J, iPS1 or H
                2 Baraffe+98 wmixed with individual extinction: from J, iPS1 or H
    LAST UPDATE: 
        16 June 2020 - Updated RA,Dec from hms/dms to deg
                
    """
    def __init__(self,filename='USco_P_M_updated_16_Jun_2020.fit',datadir='tables/',mass_type=1):
        print('Reference table is 1. Rebull et al. 2018')
        self.data=Table.read(datadir+filename, format='fits') 
        print('Loading data,RA,Dec,Prot,Mass,EBV,Disk')
        self.RA=self.data['RAJ2000']
        self.Dec=self.data['DEJ2000']
        self.Prot=self.data['Per1']
        if mass_type==0:
            print('There were no masses estimated in Rebull et al. 2018!')
            self.Mass=np.nan
        elif mass_type==1:
            self.Mass=self.data['Mass_MESA']
        elif mass_type==2:
            self.Mass=self.data['Mass_B98']      
        self.Mass[self.Mass<=0.]=np.nan 
        self.EBV=self.data['E_B-V']
        self.Disk=self.data['Disked']
        self.SpT=np.nan        
        print('PS1 for PanSTARRs DR2 stacked data')
        print('TwoMASS for 2MASS data')       
        print('ClusterInfo forgeneral info about the cluster')
    def PS1(self):
        self.g=self.data['gPSFMag']
        self.g_e=self.data['gPSFMagErr']
        self.r=self.data['rPSFMag']
        self.r_e=self.data['rPSFMagErr']
        self.i=self.data['iPSFMag']
        self.i_e=self.data['iPSFMagErr']
        self.z=self.data['zPSFMag']
        self.z_e=self.data['zPSFMagErr']
        self.y=self.data['yPSFMag']
        self.y_e=self.data['yPSFMagErr']  
    def TwoMASS(self):            
        self.Jmag=self.data['j_m']
        self.Jmag_e=self.data['j_cmsig']
        self.Hmag=self.data['h_m']
        self.Hmag_e=self.data['h_cmsig']
        self.Ksmag=self.data['k_m']
        self.Ksmag_e=self.data['k_cmsig']
    def ClusterInfo(self,datadir='tables/'):
        cluster_info=MacOSFile.pickle_load(datadir+'UpperSco_ClusterInfo.npy')
        print('Distance from Gaia DR2 -',cluster_info['Gaia_distance_ref'])
        print(cluster_info['Gaia_distance_comment'])
        self.dist=cluster_info['Gaia_distance']
        self.DM=5.*np.log10(self.dist/10)
        print(cluster_info['comment_EBV'])
        self.EBV=self.luster_info['ref_EBV']
        self.Rv=3.1
        self.Av=self.Rv*self.EBV
        print('Av=',self.Av)
        print('Median Age adopted')
        self.Age=cluster_info['Age']
        print(cluster_info['FeH_ref'])
        self.FeH=cluster_info['FeH']
        
class NGC2362:
    """
    Reference table:
        Irwin et al. 2008 http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/384/675


    Reference for disks: Julia updated it using:
                            1- If MYStIX IR-excess
                            2- If [3.6]-[8.0]>0.7 (same criteria as in Irwin+2008)
    References for Mass: Median Extinction Masses  
    
    mass_type = 0 source paper mass
                1 MESA with individual extinction: from J, iPS1 or H
                2 Baraffe+98 wmixed with individual extinction: from J, iPS1 or H
    LAST UPDATE: 16 June 2020
                
    """
    def __init__(self,filename='NGC2362_P_M_updated_16_Jun_2020.fit',datadir='tables/',mass_type=0):
        print('Reference table is 1. Irwin et al. 2008')
        self.data=Table.read(datadir+filename, format='fits') 
        print('Loading data,RA,Dec,Prot,Mass,Amp,Disk')
        self.RA=self.data['RAJ2000_1']
        self.Dec=self.data['DEJ2000_2']
        self.Prot=self.data['Per']
        self.Amp=self.data['iamp']
        if mass_type==0:
            self.Mass=self.data['Mass']
        elif mass_type==1:
            self.Mass=self.data['Mass_median_MESA']
        elif mass_type==2:
            self.Mass=self.data['Mass_median_B98']      
        self.Mass[self.Mass<=0.]=np.nan 
        self.Disk=self.data['Disked']       
        print('PS1 for PanSTARRs DR2 stacked data')
        print('TwoMASS for 2MASS data')       
        print('ClusterInfo forgeneral info about the cluster')
    def PS1(self):
        self.g=self.data['gPSFMag']
        self.g_e=self.data['gPSFMagErr']
        self.r=self.data['rPSFMag']
        self.r_e=self.data['rPSFMagErr']
        self.i=self.data['iPSFMag']
        self.i_e=self.data['iPSFMagErr']
        self.z=self.data['zPSFMag']
        self.z_e=self.data['zPSFMagErr']
        self.y=self.data['yPSFMag']
        self.y_e=self.data['yPSFMagErr']  
    def TwoMASS(self):            
        self.Jmag=self.data['_Jmag_']
        self.Jmag_e=self.data['e_Jmag_']
        self.Hmag=self.data['_Hmag_']
        self.Hmag_e=self.data['e_Hmag_']
        self.Ksmag=self.data['_Ksmag_']
        self.Ksmag_e=self.data['e_Ksmag_']
    def Spitzer(self):
        self._3_6_=self.data['_3_6_']
        self._4_5_=self.data['_4_5_']
        self._5_8_=self.data['_5_8_']
        self._8_0_=self.data['_8_0_'] 
        self._24_=self.data[ '_24_']
        self.e_3_6_=self.data[ 'e_3_6_'] 
        self.e_4_5_=self.data[ 'e_4_5_']
        self.e_5_8_=self.data['e_5_8_'] 
        self.e_8_0_=self.data['e_8_0_'] 
        self.e_24_=self.data['e_24_']
    def ClusterInfo(self,datadir='tables/'):
        cluster_info=MacOSFile.pickle_load(datadir+'NGC2362_ClusterInfo.npy')
        print('Distance from Gaia DR2 -',cluster_info['Gaia_distance_ref'])
        self.dist=cluster_info['Gaia_distance']
        self.DM=5.*np.log10(self.dist/10)
        print(cluster_info['comment_AV'])
        self.EBV=cluster_info['ref_EBV']
        print('E(B-V)=',self.EBV)
        self.Rv=3.1
        self.Av=self.Rv*self.EBV
        print('Av=',self.Av)
        print('Old Literature Age adopted')
        self.Age=cluster_info['Age']
        print(self.Age)        
        print(cluster_info['FeH_ref'])
        self.FeH=cluster_info['FeH']
        
 #I just coppied pasted this NGC 6530 code here!       
class NGC6530:
    """
    Reference table: 
        Henderson& Stassun 2012
        http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/ApJ/747/51


    Reference for disks: 
        for the M-type stars - Damiani+2019 http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A+A/623/A25
    Reference for individual AVs
        Prisinzano+ 2019 https://vizier.u-strasbg.fr/viz-bin/VizieR-3?-source=J/A%2bA/623/A159/tablea2

    References for Mass: 
    
    mass_type = For 125 stars I used the Teff-Mass scale to get masses using Teffs
    from Prisinzano+2019
    For the remaining stars I used Vmag or Imag form the source paper along with median extinction
    
    LAST UPDATE: 16 13 August 2020
                
    """
       
     
    """
    Reference table is 
    """
    def __init__(self,filename='NGC6530_P_M_updated_13_Aug_2020.fit',datadir='tables/',mass_type=0):
        print('Reference table is Henderson& Stassun 2012')
        print('Loading data,RA,Dec,Prot,Mass,EBV,Disk')        
        self.data=Table.read(datadir+filename, format='fits') 
        self.RA=self.data['RAJ2000_H12']
        self.Dec=self.data['DEJ2000_H12']
        self.Prot=self.data['Per_H12']
        self.Disk=self.data['Disked']    
        self.E_B_V=self.data['E_B_V']           
        if mass_type==0:
            print('Using Masses from Henderson et al. 2012')
            self.Mass=self.data['Mass_H12']
        elif mass_type==1:
            print('Using Masses obtained with MESA')
            self.Mass=self.data['Mass_MESA_mixed']
        elif mass_type==2:
                print("I dont have Baraffe+ masses")
        self.Mass[self.Mass<=0.]=np.nan 
        print('PS1 for PanSTARRs DR2 stacked data')
        print('TwoMASS for 2MASS data')       
        print('ClusterInfo forgeneral info about the cluster')
            
    class PS1:
        def __init__(self,filename='NGC6530_Per_PS1DR2.fits',datadir='/Users/jroquette/work/data/NGC6530/'):
            data=Table.read(datadir+filename, format='fits') 
            self.g=data['gPSFMag']
            self.g_e=data['gPSFMagErr']
            self.r=data['rPSFMag']
            self.r_e=data['rPSFMagErr']
            self.i=data['iPSFMag']
            self.i_e=data['iPSFMagErr']
            self.z=data['zPSFMag']
            self.z_e=data['zPSFMagErr']
            self.y=data['yPSFMag']
            self.y_e=data['yPSFMagErr']  
            self.cat= ac.SkyCoord(data['raStack'],data['decStack'], unit="deg")          

        print('No [Fe/H] info')
    def ClusterInfo(self,datadir='tables/'):
        cluster_info=MacOSFile.pickle_load(datadir+'NGC6530_ClusterInfo.npy')
        print('Distance from Gaia DR2 -',cluster_info['Gaia_distance_ref'])
        self.dist=cluster_info['Gaia_distance']
        print(self.dist)
        self.DM=cluster_info['DM']
        print(cluster_info['comment_AV'])
        self.Rv=cluster_info['Rv']
        print('cluster Rv=',self.Rv)
        self.Av=cluster_info['AV']
        self.EBV=cluster_info['E_B_V']
        print('E(B-V)=',self.EBV)
        print('Av=',self.Av)
        print('Literature Age adopted')
        self.Age=cluster_info['Age']
        print(self.Age)        
        print(cluster_info['FeH_ref'])
        self.FeH=cluster_info['FeH']         
   

#class CepOB3b:
#    """
#    Reference table is Littlefair+2010
#    http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/MNRAS/403/545

    # Reference for disks: Julia updated it using:
    #                         1- If MYStIX IR-excess
    #                         2- If [3.6]-[8.0]>0.7 (same criteria as in Irwin+2008)
    # References for Mass: Median Extinction Masses  
    
    # mass_type = 0 source paper mass
    #             1 MESA with individual extinction: from J, iPS1 or H
    #             2 Baraffe+98 wmixed with individual extinction: from J, iPS1 or H
    # LAST UPDATE: 18 June 2020
                
    # """
    # def __init__(self,filename='CepOB3b_Per_PS1DR2.fits',datadir='/Users/jroquette/work/data/CepOB3b/'):
    #     print('Reference table is Littlefair+2010')
    #     print('Note that their table does not provide membership, masses or disk diagnosis')
    #     self.data=Table.read(datadir+filename, format='fits') 
    #     self.RA=self.data['RAJ2000']
    #     self.Dec=self.data['DEJ2000']
    #     self.Prot=self.data['Per']
    #     self.cat= ac.SkyCoord(self.RA,self.Dec, unit="deg")  
    # def OriginalPhotometry(self):
    #     self.V=self.data['Vmag']
    #     self.eV=self.data['e_Vmag']
    #     self.Ha=self.data['Ha']
    #     self.eHa=self.data['e_Ha']   
    #     self.VI=self.data['V-I']
    #     self.eVI=self.data['e_V-I']
    #     self.RHa=self.data['R-Ha']
    #     self.eRHa=self.data['e_R-Ha']
    # class PS1:
    #     def __init__(self,filename='CepOB3b_Per_PS1DR2.fits',datadir='/Users/jroquette/work/data/CepOB3b/'):
    #         data=Table.read(datadir+filename, format='fits') 
    #         self.g=data['gPSFMag']
    #         self.g_e=data['gPSFMagErr']
    #         self.r=data['rPSFMag']
    #         self.r_e=data['rPSFMagErr']
    #         self.i=data['iPSFMag']
    #         self.i_e=data['iPSFMagErr']
    #         self.z=data['zPSFMag']
    #         self.z_e=data['zPSFMagErr']
    #         self.y=data['yPSFMag']
    #         self.y_e=data['yPSFMagErr']  
    #         self.cat= ac.SkyCoord(data['raStack'],data['decStack'], unit="deg")    
    # def ClusterInfo(self):
    #     print('Distance from Gaia DR2 -     ')
    #     self.dist=819.
    #     self.DM=rt.DisttoDM(self.dist)
    #     print('Median E(B-V) from Littlefair+2010')
    #     self.EBV=0.79
    #     self.Rv=3.26 # this is the value that Litlefair+2010 uses
    #     self.Av=self.Rv*self.EBV
    #     print('E(B-V)=',self.EBV,'Av=',self.Av)
    #     print('Age from Bell et al. 2012')
    #     self.Age=6. 
    #     print(self.Age,' Myrs')
    #     print('No [Fe/H] info')
    #     self.FeH=np.nan   
    # def readJeromePPVI(self,filename='CepOB3b.permass',datadir='/Users/jroquette/work/data/JEROME_Data_from_PPVI/'):
    #     ppvi=Table.read(datadir+filename, format='ascii') 
    #     self.M_ppvi=ppvi['M']
    #     self.P_ppvi=ppvi['P']           
                                              
                                      