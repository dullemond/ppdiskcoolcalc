# Protoplanetary Disk Cooling Time Scale Calculator

## Goal

For colleagues in the astrophysical community who study protoplanetary disks,
this little tool allows to compute the cooling and thermal relaxation time
scales of the gas, caused by the (small grain) dust. The thermal relaxation time
scale (in dimensionless form called beta=Omega_K*t_relax) is an important time
scale for the dynamics of the disk. For beta>>1 the disk behaves adiabatically,
while for beta<<1 the disk behaves locally isothermally.

Since the proper computation of t_relax is not entirely trivial, and there
has been some confusion about it, the goal of this tool is to make it easier
to make a quick computation of this value.

This tool is based on the paper

   Razor-thin dust layers in protoplanetary disks:
   Limits on the Vertical Shear Instability
   C.P. Dullemond, A. Ziampras and D. Ostertag
   A&A submitted 2022

## Usage

ppdiskcoolcalc is a python tool that can be used directly from a python
prompt or from within a python notebook. In Python type:

    from ppdiskcoolcalc import *

You can make the standard model by typing:

    model=coolcalc()

You can make a plot of beta:

    plt.figure()
    plt.loglog(model.r/au,model.beta)
    plt.xlabel('r [au]')
    plt.ylabel(r'$\beta$')

or of the surface density

    plt.figure()
    plt.loglog(model.r/au,model.sigmag)
    plt.xlabel('r [au]')
    plt.ylabel(r'$\Sigma\;[g/cm^2]$')

or of the midplane temperature

    plt.figure()
    plt.loglog(model.r/au,model.tmid)
    plt.xlabel('r [au]')
    plt.ylabel('T [K]')

To find out where the disk is unstable to the vertical shear instability (VSI):

    plt.figure()
    plt.loglog(model.r/au,model.beta,label=r'$\beta$-coefficient')
    plt.loglog(model.r/au,model.beta_vsi_limit,':',label='VSI limit')
    plt.xlabel('r [au]')
    plt.ylabel(r'$\beta$')
    plt.legend()

This model is, admittedly, extremely simplified. It can be modified by setting
any of the numerous parameters to a value other than the default. Just type

    coolcalc?

to get an overview. Example for a Herbig Ae star with a disk of 0.1 solar mass,
but a dust-to-gas ratio of 1e-3:

    model=coolcalc(mstar=2.4*MS,lstar=50*LS,mdisk=0.1*MS,dtg=1e-4)

