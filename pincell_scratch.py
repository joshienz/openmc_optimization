#matplotlib inline
import os 
library_path = "/home/nichojo/codes/openmc/endfb-viii.0-hdf5/cross_sections.xml"
os.environ["OPENMC_CROSS_SECTIONS"] = library_path
# Python imports
import openmc  # for the python API of OpenMC 

import numpy as np  # for mathematical tools
import matplotlib.pyplot as plt  # to plot results


# OpenMC Model instantiation
model = openmc.Model()



# Path to the OpenMC executable
openmc_exec = "/home/nichojo/anaconda3/envs/openmc/bin/openmc"

u235 = openmc.Nuclide('U235')
u238 = openmc.Nuclide('U238')
o16 = openmc.Nuclide('O16')
h1 = openmc.Nuclide('H1')

uo2 = openmc.Material(1,"uo2")
print(uo2)

mat = openmc.Material()
print(mat)

uo2.add_nuclide(u235,0.03)
uo2.add_nuclide(u238,.97)
uo2.add_nuclide(o16,2.0)
uo2.set_density('g/cm3',10.8)

water = openmc.model.borated_water(boron_ppm=508,density=0.72)
water.add_nuclide(h1,2.0)
water.add_nuclide(o16,1.0)
water.set_density('g/cm3', 0.72)
water.add_s_alpha_beta('c_H_in_H2O')

mats = openmc.Materials([uo2, water])
mats = openmc.Materials()
mats.append(uo2)
mats += [ water]
isinstance(mats, list)

mats.export_to_xml()

fuel_or = openmc.ZCylinder(r=0.4)
fuel_region = -fuel_or

fuel = openmc.Cell(1,'fuel')
fuel.fill = uo2
fuel.region = fuel_region

pitch = 1.25984
left = openmc.XPlane(x0=-pitch/2, boundary_type='reflective')
right = openmc.XPlane(x0=pitch/2, boundary_type='reflective')
bottom = openmc.YPlane(y0=-pitch/2, boundary_type='reflective')
top = openmc.YPlane(y0=pitch/2, boundary_type='reflective')
#up = openmc.ZPlane(z0=pitch/2,boundary_type='reflective')
#down = openmc.ZPlane(z0=-pitch/2,boundary_type='reflective')

water_region = +left & -right & +bottom & -top & +fuel_or #+down & -top & 

moderator = openmc.Cell(4, 'moderator')
moderator.fill = water
moderator.region = water_region

root = openmc.Universe(cells=(fuel, moderator))

geom = openmc.Geometry()
geom.root_universe = root
geom.export_to_xml()
point = openmc.stats.Point((0, 0, 0))
src = openmc.Source(space=point)
settings = openmc.Settings()
settings.source = src
settings.batches = 100
settings.inactive = 10
settings.particles = 1000
settings.export_to_xml()


cell_filter = openmc.CellFilter(fuel)

t = openmc.Tally(1)
t.filters = [cell_filter]

t.nuclides = ['U235']
t.scores = ['total', 'fission', 'absorption', '(n,gamma)']

tallies = openmc.Tallies([t])
tallies.export_to_xml()

openmc.run()

