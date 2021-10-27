--[[

Starting setup for FCC Aluminium

This is the starter for bigger simulations

]]--

systemid = "al"
system_title = "Aluminium test for LSMS 3"
pot_in_type = 1
pot_out_type = 1

relativity = "s"
core_relativity = "f"

num_atoms = 1
nspin = 1
mtasa = 0

xcFunctional = { 0, 1 }

iprint = 0
default_iprint = -1
print_node = 0
istop = "main"

nscf = 40
rmsTolerance = 1.0e-10

max_iterations = 5
iso_relax = true
box_relax = true

energyContour = { npts = 31,
                  grid = 2,
                  ebot = -0.6,
                  etop = 0.0,
                  eitop = 0.825,
                  eibot = 0.0025 }

a = 4.03 / 0.529177

bravais = {}
bravais[1] = { 0.5 * a, 0.5 * a, 0.0 * a }
bravais[2] = { 0.0 * a, 0.5 * a, 0.5 * a }
bravais[3] = { 0.5 * a, 0.0 * a, 0.5 * a }

site_default = { lmax = 3,
                 rLIZ = 12.5,
                 rsteps = { 89.5, 91.5, 93.2, 99.9 },
                 atom = "Al",
                 Z = 13,
                 Zc = 2,
                 Zs = 8,
                 Zv = 3,
                 rad = 2 }

mixing = { { quantity = "potential", algorithm = "broyden", mixing_parameter = 0.05 } }

numberOfMixQuantities = 0
for k, v in pairs(mixing) do
    numberOfMixQuantities = numberOfMixQuantities + 1
end

site = {}
for i = 1, num_atoms do
    site[i] = {}
end

site[1].pos = { 0, 0, 0 }
site[1].evec = { 0, 0, 1 }

for i = 1, num_atoms do
    for k, v in pairs(site_default) do
        if (site[i][k] == nil) then
            site[i][k] = v
        end
    end
end


