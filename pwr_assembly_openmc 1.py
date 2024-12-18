#! coding:utf-8

"""PWR assembly example with OpenMC.

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
uo2.set_density('g/cm3', 10.0)

# Zirconium for the rod cladding
zirconium = openmc.Material(name="zirconium")
zirconium.add_element('Zr', 1.0)
zirconium.set_density('g/cm3', 6.6)

# Water for the coolant
water = openmc.Material(name="h2o")
water.add_nuclide("H1", 2.0)
water.add_nuclide("O16", 1.0)
water.set_density("g/cm3", 1.0)
water.add_s_alpha_beta("c_H_in_H2O")  # Water thermal scattering cross-section

# Register the material list in the model
model.materials = openmc.Materials([uo2, zirconium, water])

# ============================================================================
# Geometry definition
# ============================================================================

# ----------------------------------------------------------------------------
# Fuel universe
# ----------------------------------------------------------------------------

# Define geometric shapes for the fuel and cladding
fuel_radius = openmc.ZCylinder(r=0.39)
clad_radius = openmc.ZCylinder(r=0.46)

# Define regions for the fuel and cladding
fuel_region = -fuel_radius
clad_region = +fuel_radius & -clad_radius

# Create a fuel cell
fuel = openmc.Cell(name="fuel")
fuel.fill = uo2
fuel.region = fuel_region

# Create a cladding cell
clad = openmc.Cell(name="clad")
clad.fill = zirconium
clad.region = clad_region

# Define the water region
water_region = +clad_radius

# Create a moderator cell
moderator = openmc.Cell(name="moderator")
moderator.fill = water
moderator.region = water_region

# Fuel universe
u_fuel = openmc.Universe(cells=(fuel, clad, moderator))

# ----------------------------------------------------------------------------
# Guide tube universe
# ----------------------------------------------------------------------------

# Define geometric shapes for the fuel and cladding
gt_inner_radius = openmc.ZCylinder(r=0.56)
gt_outer_radius = openmc.ZCylinder(r=0.60)

# Guide tube universe
u_gt = openmc.model.pin(
    [gt_inner_radius, gt_outer_radius],
    [water, zirconium, water]
)

# ----------------------------------------------------------------------------
# Assembly lattice
# ----------------------------------------------------------------------------

assembly_left   = openmc.XPlane(x0=-10.71, boundary_type='reflective')
assembly_right  = openmc.XPlane(x0= 10.71, boundary_type='reflective')
assembly_back   = openmc.YPlane(y0=-10.71, boundary_type='reflective')
assembly_front  = openmc.YPlane(y0= 10.71, boundary_type='reflective')

 # Instantiate a Lattice
lattice = openmc.RectLattice()
lattice.lower_left = [-10.71, -10.71]
lattice.pitch = [1.26, 1.26]
                        
lattice.universes  = [
    #  1        2       3       4       5        6       7       8       9      10     11       12      13     14      15      16      17
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 1
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 2
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 3
    [u_fuel, u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_fuel], # 4
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 5
    [u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel], # 6
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 7
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 8
    [u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel], # 9
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 10
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 11
    [u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel], # 12
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 13
    [u_fuel, u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_fuel], # 14
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_gt  , u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 15
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel], # 16
    [u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel, u_fuel]  # 17
]

lat_inner = openmc.Cell(
    name='assembly inside',
    fill=lattice,
    region=+assembly_left & -assembly_right & +assembly_back & -assembly_front
)

# fill assembly with cells
u_assembly = openmc.Universe(name='infinity assembly universe 2a')
u_assembly.add_cells([lat_inner])

# Register the geometry in the model
model.geometry = openmc.Geometry(u_assembly)

# ============================================================================
# Settings
# ============================================================================

# Settings instantiation
settings = openmc.Settings()

# Run mode
settings.run_mode = "eigenvalue"

# Run strategy
settings.batches = 50
settings.inactive = 15
settings.particles = 10000

# Temperature
settings.temperature = {'default' : 600.0}

# Photon transport
settings.photon_transport = True
settings.delayed_photon_scaling = True

# Initial source distribution
bounds = [-10.71, -10.71, -1, 10.71, 10.71, 1]
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
    "U238", "H1", "U235", "O16", "Zr90", "Zr91", "Zr92", "Zr94", "Zr96"]

# Particle filter
particle_filter = openmc.ParticleFilter(
    ["neutron", "photon", "electron", "positron"])

# Mesh filter
mesh = openmc.RegularMesh()
mesh.dimension = [500, 500]
mesh.lower_left = [-10.71, -10.71]
mesh.upper_right = [10.71, 10.71]
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

# Create tallies list
tallies = openmc.Tallies([
    tally_heating, tally_heating_energy, tally_flux, tally_flux_mesh])

# Register the tallies list in the model
model.tallies = tallies

# ============================================================================
# Plots
# ============================================================================

# Entire assembly
plot1 = openmc.Plot()
plot1.basis = 'xy'
plot1.origin = (0., 0., 0.)
plot1.width = (21.42, 21.42)
plot1.pixels = (1600, 1600)
plot1.color_by = "material"
plot1.colors = {
    water: "lightskyblue",
    zirconium: "gray",
    fuel: "goldenrod"
}

# Guide tube universe
plot2 = openmc.Plot()
plot2.basis = 'xy'
plot2.origin = (0., 0., 0.)
plot2.width = (1.26, 1.26)
plot2.pixels = (400, 400)
plot2.color_by = "material"
plot2.colors = {
    water: "lightskyblue",
    zirconium: "gray",
    fuel: "goldenrod"
}

# Fuel rod universe
plot3 = openmc.Plot()
plot3.basis = 'xy'
plot3.origin = (1.26, 0., 0.)
plot3.width = (1.26, 1.26)
plot3.pixels = (400, 400)
plot3.color_by = "material"
plot3.colors = {
    water: "lightskyblue",
    zirconium: "gray",
    fuel: "goldenrod"
}

# Register the plot list in the model
model.plots = openmc.Plots([plot1, plot2, plot3])

# ============================================================================
# Run OpenMC
# ============================================================================

# Create run folder
output_path = "pwr_assembly_openmc"
os.makedirs(output_path)
model.export_to_xml(output_path)
print('__path__')
# Plot the geometry
openmc.plot_geometry(cwd=output_path, openmc_exec=openmc_exec)

# Run
openmc.run(cwd=output_path, openmc_exec=openmc_exec)

# ============================================================================
# Plot Results
# ============================================================================

# Retrieve OpenMC results from a statepoint file
statepoint_path = "pwr_assembly_openmc/statepoint.50.h5"
sp = openmc.StatePoint(statepoint_path)

# Retrieve tally
mesh_size = 500
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

plt.savefig("pwr_assembly_openmc/mesh_flux.png")
