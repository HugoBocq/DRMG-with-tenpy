# -*- coding: utf-8 -*-

"""Definition of a model: the chiral hamiltonian."""

from tenpy.networks.site import SpinHalfSite
from tenpy.models.lattice import Triangular
from tenpy.models.model import CouplingModel, MPOModel
from tenpy.tools.params import asConfig
from copy import copy
from tenpy.networks.mpo import MPOEnvironment
from tenpy.simulations.measurement import m_energy_MPO

class CharalityOnTheTriangular(CouplingModel, MPOModel):
#my model inheritis from CouplingModel, NearestNeighbourModel, and MPOMolde
    def __init__(self, model_params):
        # the lattice defines the geometry
        self.model_params=model_params
        model_parameters=asConfig(self.model_params, "CharalityOnTheTriangular")
        
        nb_rings=model_parameters.get('nb_rings',2)
        nb_sites_per_ring=model_parameters.get('nb_sites_per_ring',2)
        JOnJx=model_parameters.get('JOnJx',0)
        
        lat = Triangular(nb_rings, nb_sites_per_ring, site=SpinHalfSite(conserve='Sz',
        sort_charge=True), bc=['periodic', 'periodic'], bc_MPS='infinite')
        #not conserve s^z in site i #infinite bc_MPS makes first lattice direction infinite
        CouplingModel.__init__(self, lat)
        # add terms of the Hamiltonian
        #exchange
        self.add_coupling(-JOnJx*0.5, 0, 'Sp', 0, 'Sm', [-1,1])  # strength, site1 w.r.t. unit cell, operator1, site2 w.r.t. unit cell, operator 2, dx 
        self.add_coupling(-JOnJx*0.5, 0, 'Sm', 0, 'Sp', [-1,1])  # Sp_i Sm_{i-1}
        self.add_coupling(-JOnJx, 0, "Sz", 0, "Sz", [-1,1])
        
       
        self.add_coupling(-JOnJx*0.5, 0, "Sp", 0, "Sm", [0,1])  
        self.add_coupling(-JOnJx*0.5, 0, "Sm", 0, "Sp", [0,1])  
        self.add_coupling(-JOnJx, 0, "Sz", 0, "Sz", [0,1])

        self.add_coupling(-JOnJx*0.5, 0, "Sp", 0, "Sm", [1,0])  
        self.add_coupling(-JOnJx*0.5, 0, "Sm", 0, "Sp", [1,0])  
        self.add_coupling(-JOnJx, 0, "Sz", 0, "Sz", [1,0])
        
        # add chirality scalar
        #top left
        j=[0,0]
        k=[-1,1]
        l=[0,1]
        self.add_multi_coupling(1j/2, [('Sz', l, 0), ('Sp',j, 0),('Sm',k,0)])
        self.add_multi_coupling(-1j/2, [('Sz', l, 0), ('Sm',j,0),('Sp',k,0)])
        self.add_multi_coupling(1j/2, [('Sz', k, 0), ('Sp',l, 0),('Sm',j,0)])
        self.add_multi_coupling(-1j/2, [('Sz', k, 0), ('Sm',l,0),('Sp',j,0)])
        self.add_multi_coupling(1j/2, [('Sz', j, 0), ('Sp',k, 0),('Sm',l,0)])
        self.add_multi_coupling(-1j/2, [('Sz', j, 0), ('Sm',k,0),('Sp',l,0)])
        
        #top right
        j=[0,0]
        k=[0,1]
        l=[1,0]
        self.add_multi_coupling(1j/2, [('Sz', l, 0), ('Sp',j, 0),('Sm',k,0)])
        self.add_multi_coupling(-1j/2, [('Sz', l, 0), ('Sm',j,0),('Sp',k,0)])
        self.add_multi_coupling(1j/2, [('Sz', k, 0), ('Sp',l, 0),('Sm',j,0)])
        self.add_multi_coupling(-1j/2, [('Sz', k, 0), ('Sm',l,0),('Sp',j,0)])
        self.add_multi_coupling(1j/2, [('Sz', j, 0), ('Sp',k, 0),('Sm',l,0)])
        self.add_multi_coupling(-1j/2, [('Sz', j, 0), ('Sm',k,0),('Sp',l,0)])
        

        # finish initialization
        # generate MPO for DMRG
        MPOModel.__init__(self, lat, self.calc_H_MPO())




def m_scalar_chirality(results, psi, model, simulation):
    """
        Measuring the scalar chirality per ring
    """
    model_params_chirality_only=copy(model.model_params)
    model_params_chirality_only['JOnJx']=0
    model_chirality_only=CharalityOnTheTriangular(model_params_exchange_only)
    results['scalar_chirality'] =model_chirality_only.H_MPO.expectation_value(psi)

