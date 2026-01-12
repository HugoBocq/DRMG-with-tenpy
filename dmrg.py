# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 17:24:24 2024

@author: Hugo2
"""
from tenpy.algorithms import dmrg
from model import CharalityOnTheTriangular
import tenpy.networks.mps as MPS

#create a model: lattice (model.lat), mpo
model_params={}
model_params['nb_rings']=2
model_params['nb_sites_per_ring']=4
model_params['JxOnJ']=0

model=CharalityOnTheTriangular(model_params)

#create an initial state: a product state of local Neel state
options = {'method': 'lat_product_state',

           'product_state' : [[['up'], ['down']],

                              [['down'], ['up']]],

           }
psi = MPS.InitialStateBuilder(model.lat, options).run()
print(model.lat.possible_couplings(0,0,[0,1]))

#run dmrg
dmrg_params = {"trunc_params": {"chi_max": 20, "svd_min": 1.e-10}, "mixer": True}
#a=dmrg.run(psi,model, dmrg_params)
print(a['E'])
