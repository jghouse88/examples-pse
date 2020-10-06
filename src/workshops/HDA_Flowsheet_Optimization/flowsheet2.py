#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:44:20 2020

@author: andrew
"""

from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import (PressureChanger,
                                        Mixer,
                                        Separator as Splitter,
                                        Heater,
                                        StoichiometricReactor,
                                        Flash)

from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom

# Import idaes logger to set output levels
import idaes.logger as idaeslog

from idaes.generic_models.properties.core.generic.generic_property import GenericParameterBlock
from BTHM_ideal import configuration
import hda_reaction as reaction_props

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.thermo_params = GenericParameterBlock(default=configuration)

m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                           "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})

m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                            "has_pressure_change": False,
                            "has_phase_equilibrium": True})

m.fs.R101 = StoichiometricReactor(
            default={"property_package": m.fs.thermo_params,
                     "reaction_package": m.fs.reaction_params,
                     "has_heat_of_reaction": False,
                     "has_heat_transfer": True,
                     "has_pressure_change": False})

m.fs.F101 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

m.fs.S101 = Splitter(default={"property_package": m.fs.thermo_params,
                               "ideal_separation": False,
                               "outlet_list": ["purge", "recycle"]})


m.fs.C101 = PressureChanger(default={
            "property_package": m.fs.thermo_params,
            "compressor": True,
            "thermodynamic_assumption": ThermodynamicAssumption.isothermal})

# m.fs.F102 = Flash(default={"property_package": m.fs.thermo_params,
#                            "has_heat_transfer": True,
#                            "has_pressure_change": True})

m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
m.fs.s09 = Arc(source=m.fs.C101.outlet,
               destination=m.fs.M101.vapor_recycle)
# m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)

print(degrees_of_freedom(m))

m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
m.fs.M101.toluene_feed.temperature.fix(303.2)
m.fs.M101.toluene_feed.pressure.fix(350000)

m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
m.fs.M101.hydrogen_feed.temperature.fix(303.2)
m.fs.M101.hydrogen_feed.pressure.fix(350000)

m.fs.H101.outlet.temperature.fix(600)

m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

m.fs.R101.conv_constraint = Constraint(
    expr=m.fs.R101.conversion*m.fs.R101.inlet.
    flow_mol_phase_comp[0, "Vap", "toluene"] ==
    (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
     m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

m.fs.R101.conversion.fix(0.75)
m.fs.R101.heat_duty.fix(0)
# m.fs.R101.outlet.temperature.fix(600)

m.fs.F101.vap_outlet.temperature.fix(325.0)
m.fs.F101.deltaP.fix(0)

# m.fs.F102.vap_outlet.temperature.fix(370)  # reduced T here
# m.fs.F102.deltaP.fix(-200000)

m.fs.S101.split_fraction[0, "purge"].fix(0.2)
m.fs.C101.outlet.pressure.fix(350000)

print(degrees_of_freedom(m))

seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 2

# Using the SD tool
G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

for o in heuristic_tear_set:
    print(o.name)

for o in order:
    print(o[0].name)

tear_guesses = {
        "flow_mol_phase_comp": {
                (0, "Vap", "benzene"): 1e-5,
                (0, "Vap", "toluene"): 1e-5,
                (0, "Vap", "hydrogen"): 0.30,
                (0, "Vap", "methane"): 0.02,
                (0, "Liq", "benzene"): 1e-5,
                (0, "Liq", "toluene"): 0.30},
        "temperature": {0: 303},
        "pressure": {0: 350000}}

# Pass the tear_guess to the SD tool
seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

def function(unit):
    unit.initialize(outlvl=idaeslog.INFO)

seq.run(m, function)

# m.fs.H101.initialize(
#     state_args={
#         "flow_mol_phase_comp": {
#                 ("Vap", "benzene"): 1e-5,
#                 ("Vap", "toluene"): 1e-5,
#                 ("Vap", "hydrogen"): 0.30,
#                 ("Vap", "methane"): 0.02,
#                 ("Liq", "benzene"): 1e-5,
#                 ("Liq", "toluene"): 0.30},
#         "temperature": 303,
#         "pressure": 350000},
#     outlvl=idaeslog.DEBUG)

solver = SolverFactory('ipopt')
solver.solve(m, tee=True)

# m.fs.F101.report()
# m.fs.F102.report()

# m.fs.thermo_params.benzene.display()
# m.fs.thermo_params.temperature_ref.display()
# # m.fs.R101.control_volume.properties_in[0].enth_mol.display()
# m.fs.R101.control_volume.properties_in[0].temperature.display()
# # m.fs.R101.control_volume.properties_out[0].temperature.display()
# m.fs.R101.control_volume.properties_in[0].pressure.display()
# # m.fs.R101.control_volume.properties_out[0].pressure.display()
# m.fs.R101.control_volume.properties_in[0].enth_mol_phase_comp.display()
# m.fs.R101.control_volume.properties_out[0].enth_mol_phase_comp.display()
# for k in m.fs.R101.control_volume.properties_in[0].enth_mol_phase_comp:
#     # print(m.fs.R101.control_volume.properties_in[0].enth_mol_phase_comp[k]._expr)
#     print(k, value(
#         m.fs.R101.control_volume.properties_in[0].enth_mol_phase_comp[k] -
#         m.fs.R101.control_volume.properties_out[0].enth_mol_phase_comp[k]))
m.fs.R101.report()
m.fs.F101.report()
