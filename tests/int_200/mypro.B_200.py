from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

print('Loading...')
pdb = PDBFile('4qvf_chainB_cleaned.pdb')
f=open('results.B_200.csv','w')
forcefield = ForceField('amber99sb.xml', 'tip3p.xml')
modeller = Modeller(pdb.topology, pdb.positions)
print('Adding hydrogens...')
modeller.addHydrogens(forcefield)
print('Adding solvent...')
modeller.addSolvent(forcefield, model='tip3p', padding=1*nanometer)
print('Minimizing...')
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
platform = Platform.getPlatformByName('OpenCL')


simulation = Simulation(modeller.topology, system, integrator,platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy()
simulation.reporters.append(PDBReporter('output.B_200.pdb', 1000))

print('Saving...')
simulation.reporters.append(StateDataReporter(f,200,step=True,kineticEnergy=True,totalEnergy=True,
						potentialEnergy=True,density=True,temperature=True))
simulation.step(30000)
print('Done')
