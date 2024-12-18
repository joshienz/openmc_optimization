#! coding:utf-8

import numpy as np
import matplotlib.pyplot as plt
import os
library_path = "/home/nichojo/codes/openmc/endfb-viii.0-hdf5/cross_sections.xml"
os.environ["OPENMC_CROSS_SECTIONS"] = library_path
import openmc
import pandas as pd
import openmc.mgxs as mgxs

model = openmc.Model()

inf_medium = openmc.Material(name='moderator')
inf_medium.set_density('g/cc', 5.)
inf_medium.add_nuclide('H1',0.028999667)
inf_medium.add_nuclide('O16', 0.01450188)
inf_medium.add_nuclide('U235', 0.000114142)
inf_medium.add_nuclide('U238', 0.006886019)

model.materials = openmc.Materials([inf_medium])

left = openmc.XPlane(boundary_type='reflective',x0=-0.63)
right = openmc.XPlane(boundary_type='reflective', x0=0.63)
up = openmc.YPlane(boundary_type='reflective', y0=-0.63)
down = openmc.YPlane(boundary_type='reflective', y0=0.63)

cell = openmc.Cell(cell_id=1,name='cell')
cell.region = +left & -right & +up & -down
cell.fill = inf_medium

root_universe = openmc.Universe(name='root universe', cells=[cell])

model.geometry = openmc.Geometry(root_universe)

batches = 50
inactive = 15
particles = 10000

settings = openmc.Settings()
settings.batches = batches
settings.inactive = inactive
settings.particles = particles
settings.generations_per_batch = 5
settings.output = {'tallies':True}

bounds = [-.63, -.63, -.63, 0.63, 0.63, 0.63]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable = True)
settings.source = openmc.IndependentSource(space= uniform_dist)

model.settings = settings

twoKGroup = np.logspace(-5,7,num=2000)
groups = mgxs.EnergyGroups(group_edges=twoKGroup)

total = mgxs.TotalXS(domain=cell,energy_groups=groups)
absorption = mgxs.AbsorptionXS(domain=cell,energy_groups=groups)
scatteringMultFactor = mgxs.ScatterXS(domain=cell, energy_groups=groups, nu=True)
scattering = mgxs.ScatterXS(domain=cell, energy_groups=groups)

scatteringMatrix = mgxs.ScatterMatrixXS(domain=cell,energy_groups=groups)
fission = mgxs.FissionXS(domain=cell,energy_groups=groups)
fissionMultFactor = mgxs.FissionXS(domain=cell,energy_groups=groups, nu=True)

print("break")
absorption.tallies

tallies = openmc.Tallies()

tallies += total.tallies.values()
tallies += absorption.tallies.values()
tallies+= scattering.tallies.values()
tallies += scatteringMultFactor.tallies.values()
tallies += scatteringMatrix.tallies.values()
tallies += fission.tallies.values()
tallies += fissionMultFactor.tallies.values()

model.tallies = tallies

statepoint_xs = model.run()

sp = openmc.StatePoint(statepoint_xs)


total.load_from_statepoint(sp)
absorption.load_from_statepoint(sp)
scattering.load_from_statepoint(sp)
scatteringMatrix.load_from_statepoint(sp)
scatteringMultFactor.load_from_statepoint(sp)
fission.load_from_statepoint(sp)
fissionMultFactor.load_from_statepoint(sp)

sp.close()

total.print_xs()
df = scattering.get_pandas_dataframe()
df.head(10)