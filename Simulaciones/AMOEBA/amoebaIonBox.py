from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout
from openmmtools import *

pdb = PDBFile(sys.argv[1])
forcefield = ForceField("amoeba2013.xml")
temperatura = 300.0
timeStep = 0.0005
stepImp = 2000
pasosTotales = 100

system = forcefield.createSystem(
    pdb.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=0.7 * nanometer,
    vdwCutoff=0.9 * nanometer,
    constraints=None,
    rigidWater=False,
    polarization="mutual",
    mutualInducedTargetEpsilon=0.00001,
)

platform = Platform.getPlatformByName("CUDA")
properties = {"CudaPrecision": "mixed"}

integrator = LangevinIntegrator(
    temperatura * kelvin, 1 / picosecond, timeStep * picoseconds
)
simulation = Simulation(pdb.topology, system, integrator, platform, properties)
simulation.context.setPositions(pdb.positions)
simulation.context.setVelocitiesToTemperature(temperatura * kelvin)

print("Minimizacion....")
simulation.minimizeEnergy()
print("Termina Minimizacion")

forces = [system.getForce(force_index) for force_index in range(system.getNumForces())]
forces = [force for force in forces if isinstance(force, openmm.AmoebaMultipoleForce)]
force = forces[0]
# fileInducedDipole = open('DipoloInducido.txt','w')
# fileTotDipole = open('DipolosTotales.txt','w')
fileSystMomMultipolar = open("MomentoMultipolarSistema.txt", "w")
# filePermanentDipoles = open('DipolosPermanentes.txt','w')


class ForceReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, "w")
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, False, True, False, None)

    def report(self, simulation, state):
        forces = state.getForces().value_in_unit(kilojoules / mole / nanometer)
        self._out.write(
            "Time: %s    Step: %s\n"
            % (str(state.getTime()), str(simulation.currentStep))
        )
        self._out.write("\n")
        for f in forces:
            self._out.write("%+20.8f\t%+20.8f\t%+20.8f\n" % (f[0], f[1], f[2]))


class VelReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, "w")
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, False, True, False, False, None)

    def report(self, simulation, state):
        pos = state.getVelocities().value_in_unit(nanometers / picoseconds)
        self._out.write(
            "Time: %s    Step: %s\n"
            % (str(state.getTime()), str(simulation.currentStep))
        )
        self._out.write("\n")
        for v in pos:
            self._out.write("%+20.8f\t%+20.8f\t%+20.8f\n" % (v[0], v[1], v[2]))


class PosReporter(object):
    def __init__(self, file, reportInterval):
        self._out = open(file, "w")
        self._reportInterval = reportInterval

    def __del__(self):
        self._out.close()

    def describeNextReport(self, simulation):
        steps = self._reportInterval - simulation.currentStep % self._reportInterval
        return (steps, True, False, False, False, False)

    def report(self, simulation, state):
        forces = state.getPositions().value_in_unit(nanometer)
        self._out.write(
            "Time: %s    Step: %s\n"
            % (str(state.getTime()), str(simulation.currentStep))
        )
        self._out.write("\n")
        for r in forces:
            self._out.write(
                "%+20.8f\t%+20.8f\t%+20.8f\n" % (r[0] * 10.0, r[1] * 10.0, r[2] * 10.0)
            )


def imprimeMultipolos(filename, arrangesize, currentStep, timeStep):
    filename.write(
        "Time:  %f ps  \t Step: %f\n\n" % (currentStep * timeStep, currentStep)
    )
    for i in arrangesize:
        for j in i:
            filename.write(
                "{:+16.12f}\t".format(j),
            )
        filename.write("\n ")


# Momento multipolar del sistema
def imprimeMomentoMultipolar(filename, arrangesize, currentStep, timeStep):
    filename.write(
        "Time:  %f ps  \t Step: %f\n\n" % (currentStep * timeStep, currentStep)
    )
    contador = 1
    for i in arrangesize:
        if contador == 1:
            filename.write(
                "{:+16.12f}\t".format(i),
            )
            filename.write("\n ")
            contador = 2
        else:
            filename.write(
                "{:+16.12f}\t".format(i),
            )
            if contador == 4 or contador == 7 or contador == 10 or contador == 13:
                filename.write("\n ")
            contador = contador + 1

#stepImp=1
# simulation.reporters.append(VelReporter('velocidades.dat', stepImp))
# simulation.reporters.append(ForceReporter('fuerzas.dat',stepImp))
simulation.reporters.append(PDBReporter("trayectoria.pdb", stepImp))
simulation.reporters.append(
    StateDataReporter(
        "dinamica.dat",
        stepImp,
        step=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        speed=True,
        time=True,
        remainingTime=True,
        totalSteps=pasosTotales,
        separator="\t",
    )
)
print("Simulacion inicia")
#simulation.step(10)

#exit(0)
for i in range(pasosTotales):
    simulation.saveState("output.xml")
    with open("system.xml", "w") as outfile:
        outfile.write(XmlSerializer.serialize(system))
    simulation.step(stepImp)
    print(stepImp*(i+1)*timeStep)

    # Dipolos inducidos
    #    inducedDipoleMoments = force.getInducedDipoles(simulation.context)
    #    imprimeMultipolos(fileInducedDipole,inducedDipoleMoments,simulation.currentStep,timestep)

    # Dipolos totales
    #    totaldipoles = force.getTotalDipoles(simulation.context)
    #    imprimeMultipolos(fileTotDipole,totaldipoles,simulation.currentStep,timestep)

    # Momento multipolar del sistema
    systemMultipolarMoment = force.getSystemMultipoleMoments(simulation.context)
    imprimeMomentoMultipolar(
        fileSystMomMultipolar,
        systemMultipolarMoment,
        simulation.currentStep,
        timeStep,
    )

    # Dipolos Permanentes
    #    permanentDipole = force.getLabFramePermanentDipoles(simulation.context)
    #    imprimeMultipolos(filePermanentDipoles,permanentDipole,simulation.currentStep,timestep)

print("Simulacion Termina")