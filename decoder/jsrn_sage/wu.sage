r"""Wu list decoder for Reed-Solomon and Binary Goppa codes.
Parameter choices. See also ratinter for the core rational
interpolation functionality."""
#Copyright Johan S. R. Nielsen
#This file is part of the codinglib for Sage.
#This library is free software released under the GPL v3 or any later version as
#published by the Free Software Foundation.


from codinglib.ratinter import *


#####################################################
### CHOOSING PARAMETERS
### specialised for decoding various codes
#####################################################

# GRS codes
def wu_w(n,k,tau):
    """For GRS codes, return the w factor for these parameters"""
    return 2*tau-(n-k+1)

def wu_params_rs(n,k,tau):
    """For GRS codes, return my closed expressions for satisfiable s and l,
    from the Goppa paper"""
    return rat_params(n, wu_w(n,k,tau), tau)

def wu_params_old_rs(n,k,tau):
    """For GRS codes, calculate the s and l parameters given n, k and tau,
    using the closed expressions of Wu's original paper (and my thesis)"""
    return rat_params_old(n, wu_w(n,k,tau), tau)

def wu_params_trifonov_rs(n,k,tau):
    """For GRS codes, return Trivonov's values for the multiplicity and
    list size."""
    return rat_params(n, wu_w(n,k,tau), tau)

def wu_params_min_rs(n,k,tau):
    """For GRS codes, return minimal values for the multiplicity and
    list size."""
    return rat_params_min_list(n, wu_w(n,k,tau), tau)

def wu_system_size_rs(n,k,tau,(s,l)):
    """For GRS codes, return the size of the linear system to solve for the
    given stats."""
    return rat_system_size(n, wu_w(n,k,tau), tau, (s,l))

def wu_stats_rs(n,k,tau):
    """For GRS codes, return some stats on the rational interpolation problem
    to solve"""
    return rat_stats(n, wu_w(n,k,tau), tau)

def wu_decoding_radius_rs(n,k, l=None, s=None):
    """Return the maximal decoding radius of the Wu decoder and minimal
    parameters attaining this. If ell is set, the maximal radius using this
    ell; same for s."""
    def get_tau(s,l):
        return gs_decoding_radius(n,k,l=l,s=l-s)
    if l==None and s==None:
        tau = codinglib.listdecoding.list_decoding_range(n,n-k+1)[1]
        return (tau, wu_params_min_rs(n,k,tau))
    if l!=None and s!=None:
        return (get_tau(s,l), (s,l))
    if s!= None:
        raise NotImplemented
    if l!= None:
        raise NotImplemented


# Goppa codes

def wu_params_goppa(n,gdeg,tau):
    """Return my closed expressions for satisfiable s and l, from the Goppa
    paper"""
    return rat_params(n,tau-gdeg-1/2,tau)

def wu_params_min_goppa(n,gdeg,tau):
    """For Goppa codes, return minimal values for the multiplicity and
    list size."""
    return rat_params_min_list(n, tau-gdeg-1/2, tau)

def wu_system_size_goppa(n,gdeg,tau,(s,l)):
    """For Goppa codes, return the size of the linear system to solve for the
    given stats."""
    return rat_system_size(n, tau-gdeg-1/2, tau, (s,l))

def wu_stats_goppa(n,gdeg,tau):
    """For Goppa codes, return some stats on the rational interpolation problem
    to solve"""
    return rat_stats(n,tau-gdeg-1/2, tau)
