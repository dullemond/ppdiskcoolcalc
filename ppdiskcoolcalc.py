#
# This is a (hopefully user-friendly) Python tool to compute the cooling
# time and the relaxation time in a protoplanetary disk. Based on the
# paper "Razor-thin dust layers in protoplanetary disks:
# Limits on the Vertical Shear Instability" by Dullemond, Ziampras
# & Ostertag (A&A submitted). 
#
# (c) 2022 Cornelis Dullemond
#
# 2022-08-27
#

import numpy as np
import matplotlib.pyplot as plt
from natconst import *


class coolcalc(object):
    def __init__(self,
                 mstar=MS,          # Stellar mass [g]
                 lstar=LS,          # Stellar luminosity [erg/s]
                 rin=1*au,          # Inner radius of the grid [cm]
                 rout=1e2*au,       # Outer radius of the grid [cm]
                 nr=100,            # Number of radial grid points in r-grid
                 r0=1*au,           # Some fiducial radius r0 [cm] (to scale sigmag0 etc)
                 sigmag0=100,       # Sigma_gas at r=r0 [g/cm^2]
                 plsig=-1,          # Powerlaw of Sigma_gas(r) ~ r^plsig
                 mdisk=None,        # Mass of the disk (optional), overrides sigmag0
                 flang=0.05,        # Flaring angle to compute midplane temperature tmid(r)
                 tmid0=None,        # T_mid(r=r0) [K], overrides flang
                 pltmid=None,       # Powerlaw of T_mid(r) ~ r^pltmid, overrides flang
                 gamma=1.6666,      # Adiabatic coefficient
                 khat=6.2832,       # The dimensionless k_x for the opt thick relaxation time
                 agrain=1e-5,       # The grain radius in cm
                 dustmdens=1.97,    # The grain material density in g/cm^3
                 dtg=1e-2,          # Small grain dust-to-gas ratio
                 opac='default',    # Opacity model. 'default'=standard, 'belllin'=Bell & Lin (1997), 'optool' uses optool
                 optool='optool -c pyr-mg70 0.696 -c c-z 0.104 -m h2o-w 0.2 -a 0.1'
                 ):
        """Tool to compute the cooling time and the relaxation time for
        a protoplanetary disk. Particularly useful for computing where
        the disk can become unstable to the VSI (Vertical Shear Instability).

        Optional parameters (all have default values):

          mstar    : default = MS = Solar mass
                     Stellar mass [g]
             
          lstar    : default = LS = Solar luminosity
                     Stellar luminosity [erg/s]
             
          rin      : default = 1 au
                     Inner radius of the grid [cm]
     
          rout     : default = 1e3 au
                     Outer radius of the grid [cm]
             
          nr       : default = 100
                     Number of radial grid points in the r-grid
             
          r0       : default = 1 au
                     Some fiducial radius r0 [cm] (to scale sigmag0 etc)
          
          sigmag0  : default = 100 g/cm^2
                     Gas surface density Sigma_gas at r=r0 [g/cm^2]
          
          plsig    : default = -1
                     Powerlaw of Sigma_gas(r) ~ r^plsig
          
          mdisk    : default = None
                     Mass of the disk (optional). If this is set, it overrides
                     the value of sigmag0
          
          flang    : defaut = 0.05
                     Flaring angle used to compute midplane temperature T_mid(r).
  
          tmid0    : default = None
                     T_mid(r=r0) [K], overrides flang
          
          pltmid   : default = None
                     Powerlaw of T_mid(r) ~ r^pltmid, overrides flang

          gamma    : defaul = 5/3
                     Adiabatic coefficient. Default is 5/3 because in the
                     cold outer disk regions the temperature is not hot
                     enough to excite rotation of H2. If it is, then 7/5
                     is more appropriate.
        
          dtg      : default = 1e-2
                     Small grain dust-to-gas ratio

          khat     : default = 2*pi
                     The dimensionless wavenumber for the perturbation analysis
                     for the optically thick relaxation time such that
                     k_x = khat / h_p where h_p is the pressure scale height.

          agrain   : default = 1e-1 micron = 1e-5 cm
                     The grain radius [cm]

          dustmdens: default = 1.97 g/cm^3
                     The dust material density [g/cm^3]
        
          opac     : default = None
                     Opacity model for the small-grained dust. Possibilities:
                      'default': Standard opacity model of the Dullemond, Ziampras,
                                 Ostertag (2022) paper.
                      'belllin': Bell & Lin (1997) opacity model
                      'optool':  Use Carsten Dominik's optool code (see below)
                                 to compute the opacity model
  
          optool   : default = 'optool -c pyr-mg70 0.696 -c c-z 0.104 -m h2o-w 0.2 -a 0.1'
                     (only when opac=='optool', and make sure to install optool,
                     see https://github.com/cdominik/optool). This is the command
                     to create the small-grain dust opacity with optool. You can
                     change this command line as you see fit. All options are
                     described in the optool manual or with optool --help. Also the
                     grain size is specified in this command line.

        This object computes:

          self.tcool              : The radiative cooling time [s]
          self.trelax             : The radiative relaxation time [s] including everything
          self.trelax_vsi_limit   : The critical relaxation time for VSI [s]
          self.beta               : = Omega_K t_relax
          self.beta_vsi_limit     : = Omega_K t_relax_vsi_limit

        Internal calculations also compute:

          self.trelax_thin        : The optically thin relaxation time [s]
          self.trelax_thick       : The optically thin relaxation time [s]
          self.t_dustgas          : The dust-gas coupling time [s]
          self.kapplaw            : The value of dlog(kappa_P)/dlog(T)

        You can also modify self.sigmag, self.sigmad, and/or self.tmid after the
        initial call. You then must call self.recompute() to recompute everything.
        Example:

          model = coolcalc()
          model.sigmad *= (model.r/(30*au))**(-0.2)
          model.recompute()

        """
        self.mstar     = mstar
        self.lstar     = lstar
        self.gamma     = gamma
        self.khat      = khat
        self.agrain    = agrain
        self.dustmdens = dustmdens

        # Setup the grid

        ri         = rin * (rout/rin)**np.linspace(0,1,nr+1)
        r          = np.sqrt(ri[:-1]*ri[1:])
        dsurf      = np.pi*(ri[1:]**2-ri[:-1]**2)
        self.r     = r
        self.ri    = ri
        self.dsurf = dsurf

        # Compute the midplane temperature

        if flang is not None:
            tmid   = (0.5*flang*lstar/(4*pi*r**2*ss))**0.25
        if tmid0 is not None and pltmid is not None:
            tmid   = tmid0 * (r/r0)**pltmid
        self.tmid  = tmid

        # Set up the sigma_gas (a powerlaw model)

        if mdisk is not None:
            sigmag = (r/r0)**plsig
            md     = (sigmag*dsurf).sum()
            sigmag *= mdisk/md
        elif sigmag0 is not None:
            sigmag = sigmag0 * (r/r0)**plsig
        else:
            raise ValueError('Could not compute Sigma_gas(r): Too little information')
        self.sigmag = sigmag
        self.mdisk  = (sigmag*dsurf).sum()

        # Set up the dust

        self.sigmad = self.sigmag * dtg
        
        # Set up the Rosseland and Planck mean opacities

        if opac=='default':
            # Use opacity model of Dullemond, Ziampras & Ostertag (2022)
            assert self.agrain<1e-3, 'The opacity model of Dullemond, Ziampras & Ostertag (2022) is only valid for agrain<10 micron.'
            kapplaw= 1.6
            kap_p  = (tmid/1.7)**kapplaw
            kap_r  = (tmid/2.25)**kapplaw
        elif opac=='belllin':
            # Use opacity model of Bell & Lin
            assert self.agrain<1e-3, 'The Bell & Lin opacity model is only valid for agrain<10 micron.'
            assert self.tmid<150, 'The Bell & Lin opacity model used here is only valid for T<150 K.'
            kapplaw= 2
            kap_p  = 2e-2 * tmid**kapplaw
            kap_r  = 2e-2 * tmid**kapplaw
            assert tmid<150, 'This Bell & Lin opacity powerlaw is only valid for T<150K'
        elif opac=='optool':
            raise ValueError('For now the optool option is not implemented yet.')
        self.kapplaw = kapplaw
        self.kap_p   = kap_p
        self.kap_r   = kap_r

        # Now do the calculations

        self.compute_kepler_cs_hp_cv()
        self.compute_vsi_trelax_limit()
        self.compute_cooling_radiative()

    def compute_kepler_cs_hp_cv(self):
        self.omk    = np.sqrt(GG*self.mstar/self.r**3)
        self.vk     = self.omk*self.r
        self.cs     = np.sqrt(kk*self.tmid/(2.3*mp))
        self.hp     = self.cs/self.omk
        self.cvg    = kk/((self.gamma-1)*2.3*mp) # c_V of the gas
        self.cvd    = 1e7  # Rough guess of c_V of the dust

    def compute_q_tmid(self):
        q_tmid      = (np.log(self.tmid[2:])-np.log(self.tmid[:-2]))/ \
                      (np.log(self.r[2:])-np.log(self.r[:-2]))
        q_tmid      = np.hstack([q_tmid[0],q_tmid,q_tmid[-1]])
        self.q_tmid = q_tmid

    def compute_rhomid(self):
        self.rhomidg = self.sigmag / ( np.sqrt(2*pi) * self.hp )
        self.rhomidd = self.sigmad / ( np.sqrt(2*pi) * self.hp )

    def compute_vsi_trelax_limit(self):
        self.compute_q_tmid()
        self.trelax_vsi_limit = np.abs(self.q_tmid)*(self.hp/self.r)/((self.gamma-1)*self.omk)
        self.beta_vsi_limit   = self.omk*self.trelax_vsi_limit

    def compute_dust_gas_coupling(self):
        self.compute_rhomid()
        mgr         = (4*pi/3)*self.dustmdens*self.agrain**3
        sgr         = 4*pi*self.agrain**2
        eth_g       = self.cvg*self.tmid
        dtg         = self.sigmad/self.sigmag
        gamma       = self.gamma
        tmid        = self.tmid
        Cbar_H      = ((gamma+1)/(gamma-1)) * np.sqrt(kk*tmid/(2*pi*2.3*mp)) \
                      * kk / (2*2.3*mp)
        srfpermass  = sgr/mgr
        qdg         = Cbar_H * srfpermass * self.rhomidg * tmid # Exchange of heat per gram of dust
        self.t_dustgas = eth_g /qdg / dtg
        self.Cbar_H = Cbar_H
        
    def compute_cooling_radiative(self):
        self.compute_rhomid()
        self.compute_dust_gas_coupling()
        qcool_thin  = 4*self.kap_p*ss*self.tmid**4
        mugas       = 2.3
        eth_g       = self.cvg*self.tmid
        dtg         = self.sigmad/self.sigmag
        tcool_thin  = (1/dtg)*eth_g/qcool_thin
        b           = self.kapplaw
        trelax_thin = tcool_thin/(4+b)
        k           = self.khat/self.hp
        eta         = 16*ss*self.tmid**3 / (3*self.kap_r*self.rhomidd*self.rhomidg*self.cvg)
        trelax_thick= 1/(eta*k**2)
        self.qcool_thin   = qcool_thin
        self.trelax_thin  = trelax_thin
        self.trelax_thick = trelax_thick
        self.trelax       = trelax_thin + trelax_thick + self.t_dustgas
        self.beta         = self.omk*self.trelax

    def recompute(self):
        self.mdisk  = (self.sigmag*self.dsurf).sum()
        self.compute_kepler_cs_hp_cv()
        self.compute_vsi_trelax_limit()
        self.compute_cooling_radiative()
