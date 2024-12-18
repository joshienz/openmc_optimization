#! coding:utf-8

"""PWR pincell example with OpenMC.

The geometry is based on the VERA benchmark problem 2b without the gap in
the rod.
"""
import os  # to create folder and environment variable
# Set the path to the cross-sections data
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
# ============================================================================
# Material definition
# ============================================================================

# Uranium dioxide (UO2) enriched at 3ao% for the fuel
uo2 = openmc.Material(1, "uo2")
uo2.add_nuclide('U235', 0.03)
uo2.add_nuclide('U238', 0.97)
uo2.add_nuclide('O16', 2.0)
uo2.set_density('g/cm3', 10.8)

# Zirconium for the rod cladding
#zirconium = openmc.Material(name="zirconium")
#zirconium.add_element('Zr', 1.0)
#zirconium.set_density('g/cm3', 6.6)

# Water for the coolant
water = openmc.model.borated_water(boron_ppm=508)
water.add_nuclide("O16", 1.0)
water.set_density("g/cm3", 0.72)
water.add_s_alpha_beta("c_H_in_H2O")  # Water thermal scattering cross-section

# Register the material list in the model
model.materials = openmc.Materials([uo2, water])

# ============================================================================
# Geometry definition
# ============================================================================

# Define geometric shapes for the fuel and cladding
fuel_radius = openmc.ZCylinder(r=0.4)
#clad_radius = openmc.ZCylinder(r=0.46)

# Define regions for the fuel and cladding
fuel_region = -fuel_radius
#clad_region = +fuel_radius & -clad_radius

# Create a fuel cell
fuel = openmc.Cell(name="fuel")
fuel.fill = uo2
fuel.region = fuel_region

# Create a cladding cell
#clad = openmc.Cell(name="clad")
#clad.fill = zirconium
#clad.region = clad_region

# Define the box containing the pincell
pitch = 1.25984
box = openmc.model.RectangularPrism(
    width=pitch, height=pitch, boundary_type="reflective")

# Define the water region
water_region = -box

# Create a moderator cell
moderator = openmc.Cell(name="moderator")
moderator.fill = water
moderator.region = water_region

# Define the root universe from the cells
root_universe = openmc.Universe(cells=(fuel, moderator))

# Register the geometry in the model
model.geometry = openmc.Geometry(root_universe)

# ============================================================================
# Settings
# ============================================================================

# Settings instantiation
settings = openmc.Settings()

# Run mode
settings.run_mode = "eigenvalue"

# Run strategy
settings.batches = 150
settings.inactive = 15
settings.particles = 10000

# Temperature
settings.temperature = {'default' : 600.0}

# Photon transport
settings.photon_transport = True
settings.delayed_photon_scaling = True

# Initial source distribution
bounds = [-0.63, -0.63, -1, 0.63, 0.63, 1]
uniform_dist = openmc.stats.Box(bounds[:3], bounds[3:], only_fissionable=True)
settings.source = openmc.source.Source(space=uniform_dist)

# Register the settings in the model
model.settings = settings

# ============================================================================
# Tallies
# ============================================================================

# Energy filter
min_energy = 1.0E-5
max_energy = 3.0E+7
energy_number = 11
energy_filter = openmc.EnergyFilter(
    np.logspace(np.log10(min_energy), np.log10(max_energy), energy_number))

# Nuclide filter
nuclides = [
    "U238", "H1", "U235", "O16"]

# Particle filter
particle_filter = openmc.ParticleFilter(
    ["neutron", "photon", "electron", "positron"])

# Mesh filter
mesh = openmc.RegularMesh()
mesh.dimension = [100, 100]
mesh.lower_left = [-0.63, -0.63]
mesh.upper_right = [0.63, 0.63]
mesh_filter = openmc.MeshFilter(mesh)

# Flux tally
tally_flux = openmc.Tally(name='flux')
tally_flux.scores = ['flux']

# Flux tally with a mesh
tally_flux_mesh = openmc.Tally(name='flux-mesh')
tally_flux_mesh.filters = [mesh_filter]
tally_flux_mesh.scores = ['flux']

# Heating tally per nuclides
tally_heating = openmc.Tally(name='heating')
tally_heating.nuclides = nuclides
tally_heating.scores = ['heating']

# Heating tally per nuclides for the energy grid
tally_heating_energy = openmc.Tally(name='heating-energy')
tally_heating_energy.filters = [energy_filter]
tally_heating_energy.nuclides = nuclides
tally_heating_energy.scores = ['heating']

# Create tally list
tallies = openmc.Tallies([
    tally_heating, tally_heating_energy, tally_flux, tally_flux_mesh])

# Register the tally list in the model
model.tallies = tallies

# ============================================================================
# Plots
# ============================================================================

# Simple 2D plot
plot = openmc.Plot()
plot.basis = 'xy'
plot.origin = (0., 0., 0.)
plot.width = (1.26, 1.26)
plot.pixels = (800, 800)
plot.color_by = "material"
plot.colors = {
    water: "lightskyblue",
    fuel: "goldenrod"
}

# Register the plot list in the model
model.plots = openmc.Plots([plot])

# ============================================================================
# Run OpenMC
# ============================================================================

# Create run folder
output_path = "pwr_pincell_openmc"
os.makedirs(output_path)
model.export_to_xml(output_path)



# Plot the geometry
openmc.plot_geometry(cwd=output_path, openmc_exec=openmc_exec)

# Run
openmc.run(cwd=output_path, openmc_exec=openmc_exec)

# ============================================================================
# Plot Results
# ============================================================================

# Retrieve OpenMC results from a statepoint file
statepoint_path = "pwr_pincell_openmc/statepoint.50.h5"
sp = openmc.StatePoint(statepoint_path)

# Retrieve tally
mesh_size = 100
tally_flux = sp.get_tally(id=2)
flux = tally_flux.get_slice(scores=["flux"])
flux.mean.shape = (mesh_size,mesh_size)
flux.std_dev.shape = (mesh_size,mesh_size)

# Heating values
fig = plt.figure(figsize=(12, 5))
fig = plt.subplot(121)
plt.imshow(flux.mean)
plt.tick_params(left=False, labelleft=False, labelbottom=False, bottom=False)
clbar = plt.colorbar()
clbar.ax.set_ylabel("Flux")

# Standard deviations
fig = plt.subplot(122)
plt.imshow(flux.std_dev, cmap='coolwarm')
plt.tick_params(left=False, labelleft=False, labelbottom=False, bottom=False)
clbar = plt.colorbar()
clbar.ax.set_ylabel("Standard deviation")
plt.tight_layout()

plt.savefig("pwr_pincell_openmc/mesh_flux.png")
