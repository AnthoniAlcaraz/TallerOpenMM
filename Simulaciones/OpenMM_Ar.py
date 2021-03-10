#!/usr/bin/env python
# -*- coding: utf-8 -*-
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools import integrators,testsystems
import os
import sys
from sys import stdout

argon = testsystems.IdealGas()
########################################################################
# Datos a modificar para las simulaciones
cutoff = 206.227*0.5  # en nm
ensamble = "NVE"  # tipos: NVE, NVT, NPT
dt = 0.002  # picosegundos
T = 298.0  # kelvin
P = 1.0  # atmosferas
intervaloBarost = 25  # frecuencia de activacion del barostato
stepImp = 50  # frecuencia de impresion (adimensional)
pasosTotales = 1000  # este valor por stepImp por dt da la longitud de la trayectoria
VelocityVerlet = False  # Solo para NVE ponerlo en True
########################################################################

system = argon.system
if ensamble == "NVE":
    if VelocityVerlet == True:
        # utilizando openmmtools
        integrator = integrators.VelocityVerletIntegrator(dt * picoseconds)
        integrator.addConstrainPositions()
        integrator.addConstrainVelocities()
    else:
        integrator = VerletIntegrator(dt * picoseconds)

elif ensamble == "NVT":
    # temperatura, frecuencia de colision, paso de integracion
    integrator = LangevinIntegrator(T * kelvin, 1 / picosecond, dt * picoseconds)
elif ensamble == "NPT":
    # temperatura, frecuencia de colision, paso de integracion
    integrator = LangevinIntegrator(T * kelvin, 1 / picosecond, dt * picoseconds)
    # Presion, temperatura y frecuencia del barostato
    system.addForce(MonteCarloBarostat(P * atmospheres, T * kelvin, intervaloBarost))
else:
    print("Ensamble no reconocido {:}.....".format(ensamble))
    sys.exit(0)

# Util cuando los enlaces estan restringidos
integrator.setConstraintTolerance(1e-8)

platform = Platform.getPlatformByName("OpenCL")
properties = {"OpenCLPrecision": "mixed"}

simulation = Simulation(argon.topology, system, integrator)
# simulation = Simulation(argon.topology, system, integrator, platform, properties)
simulation.context.setPositions(argon.positions)
simulation.context.setVelocitiesToTemperature(T * kelvin)

state = simulation.context.getState(getEnergy=True,getPositions=True)
# energiaPotencial=state.getPotentialEnergy().value_in_unit(kilojoule / mole)
posiciones = state.getPositions().value_in_unit(angstrom)
app.PDBFile.writeFile(argon.topology,posiciones, open("configuracionOptimizada.pdb", 'w'))

# Minimizacion
#simulation.minimizeEnergy(tolerance=0.0001 * kilojoule / mole, maxIterations=1000)

# Equilibrio
# simulation.step(1000)

# Informacion a guardar de la trayectoria
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


# simulation.reporters.append(VelReporter("velocidades.dat", stepImp))
# simulation.reporters.append(ForceReporter('fuerzas.dat',stepImp))
# simulation.reporters.append(PosReporter('posiciones.xyz',stepImp))
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

# Corre la trayectoria
for i in range(pasosTotales):
    # Guarda posiciones y velocidades
    simulation.saveState("output.xml")
    # Guarda el contexto de la simulacion
    with open("system.xml", "w") as outfile:
        outfile.write(XmlSerializer.serialize(system))
    simulation.step(stepImp)
