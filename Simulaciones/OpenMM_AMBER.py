#!/usr/bin/env python
# -*- coding: utf-8 -*-
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from openmmtools import integrators
import os
import sys
from sys import stdout

# Archivo leido por consola
prmtop = AmberPrmtopFile(sys.argv[1])
inpcrd = AmberInpcrdFile(sys.argv[2])

# Datos a modificar para las simulaciones
cutoff = 0.9
ensamble = "NVE"
dt = 0.002  # picosegundos
T = 300.0  # kelvin
P = 1.0  # atmosferas
intervaloBarost = 25  # frecuencia de activacion del barostato
stepImp = 50  # frecuencia de impresion (adimensional)
pasosTotales = 1000  # este valor por stepImp por dt da la longitud de la trayectoria
VelocityVerlet = False  # Solo para NVE ponerlo en True
####
system = prmtop.createSystem(
    nonbondedMethod=PME, nonbondedCutoff=cutoff * nanometer, constraints=None
)

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

# simulation = Simulation(prmtop.topology, system, integrator)
simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocitiesToTemperature(T * kelvin)

if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

# Minimizacion
simulation.minimizeEnergy()

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
