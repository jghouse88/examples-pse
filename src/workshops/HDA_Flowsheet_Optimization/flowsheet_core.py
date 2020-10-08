#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:44:20 2020

@author: andrew
"""

import pyomo.environ as pyo

from idaes.core import FlowsheetBlock

m = pyo.ConcreteModel()

m.fs = FlowsheetBlock(default={"dynamic": False})

from idaes.generic_models.properties.activity_coeff_models.\
    BTX_activity_coeff_VLE import BTXParameterBlock
import hda_reaction as reaction_props
from hda_ideal_VLE import HDAParameterBlock

m.fs.bthm_params = HDAParameterBlock()
m.fs.bt_params = BTXParameterBlock(default={
    "valid_phase": ('Liq', 'Vap'),
    "activity_coeff_model": "Ideal"})

m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.bthm_params})

# -----------------------------------------------------------------------------
from idaes.generic_models.unit_models import (Feed,
                                              Flash,
                                              Heater,
                                              Mixer,
                                              Product,
                                              PressureChanger,
                                              Separator as Splitter,
                                              StoichiometricReactor,
                                              Translator)
from idaes.generic_models.unit_models.distillation import TrayColumn

from idaes.generic_models.unit_models.distillation.condenser \
    import CondenserType, TemperatureSpec
from idaes.generic_models.unit_models.pressure_changer \
    import ThermodynamicAssumption

m.fs.hydrogen_feed = Feed(default={
    "property_package": m.fs.bthm_params})

m.fs.toluene_feed = Feed(default={
    "property_package": m.fs.bthm_params})

m.fs.M101 = Mixer(default={
    "property_package": m.fs.bthm_params,
    "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})

m.fs.H101 = Heater(default={
    "property_package": m.fs.bthm_params,
    "has_pressure_change": False,
    "has_phase_equilibrium": True})

m.fs.R101 = StoichiometricReactor(default={
    "property_package": m.fs.bthm_params,
    "reaction_package": m.fs.reaction_params,
    "has_heat_of_reaction": True,
    "has_heat_transfer": True,
    "has_pressure_change": False})

# -----------------------------------------------------------------------------
m.fs.F101 = Flash(default={"property_package": m.fs.bthm_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

m.fs.S101 = Splitter(default={"property_package": m.fs.bthm_params,
                               "ideal_separation": False,
                               "outlet_list": ["purge", "recycle"]})


m.fs.P101 = PressureChanger(default={
            "property_package": m.fs.bthm_params,
            "compressor": True,
            "thermodynamic_assumption": ThermodynamicAssumption.isothermal})

m.fs.purge = Product(default={
    "property_package": m.fs.bthm_params})

# -----------------------------------------------------------------------------
m.fs.translator = Translator(default={
    "inlet_property_package": m.fs.bthm_params,
    "outlet_property_package": m.fs.bt_params})

# Translator constraints linking outlet state variables to inlet
# state variables
m.fs.translator.eq_total_flow = pyo.Constraint(
    expr=m.fs.translator.outlet.flow_mol[0] ==
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"])

m.fs.translator.eq_temperature = pyo.Constraint(
    expr=m.fs.translator.outlet.temperature[0] ==
    m.fs.translator.inlet.temperature[0])

m.fs.translator.eq_pressure = pyo.Constraint(
    expr=m.fs.translator.outlet.pressure[0] ==
    m.fs.translator.inlet.pressure[0])

m.fs.translator.eq_mole_frac_benzene = pyo.Constraint(
    expr=m.fs.translator.outlet.mole_frac_comp[0, "benzene"] ==
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] /
    (m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
     m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))

m.fs.translator.eq_mole_frac_toluene = pyo.Constraint(
    expr=m.fs.translator.outlet.mole_frac_comp[0, "toluene"] ==
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"] /
    (m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
     m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))

# -----------------------------------------------------------------------------
m.fs.H102 = Heater(default={
    "property_package": m.fs.bt_params,
    "has_pressure_change": True,
    "has_phase_equilibrium": True})

# m.fs.C101 = TrayColumn(default={
#     "property_package": m.fs.bt_params,
#     "number_of_trays": 10,
#     "feed_tray_location": 5,
#     "condenser_type": CondenserType.totalCondenser,
#     "condenser_temperature_spec": TemperatureSpec.atBubblePoint,
#     "has_heat_transfer": False,
#     "has_pressure_change": False})

# m.fs.benzene_product = Product(default={
#     "property_package": m.fs.bt_params})
# m.fs.toluene_product = Product(default={
#     "property_package": m.fs.bt_params})

# -----------------------------------------------------------------------------
from pyomo.network import Arc, SequentialDecomposition

m.fs.s01 = Arc(source=m.fs.hydrogen_feed.outlet,
               destination=m.fs.M101.hydrogen_feed)
m.fs.s02 = Arc(source=m.fs.toluene_feed.outlet,
               destination=m.fs.M101.toluene_feed)
m.fs.s03 = Arc(source=m.fs.M101.outlet,
               destination=m.fs.H101.inlet)
m.fs.s04 = Arc(source=m.fs.H101.outlet,
               destination=m.fs.R101.inlet)
m.fs.s05 = Arc(source=m.fs.R101.outlet,
               destination=m.fs.F101.inlet)
m.fs.s06 = Arc(source=m.fs.F101.vap_outlet,
               destination=m.fs.S101.inlet)
m.fs.s07 = Arc(source=m.fs.S101.purge,
               destination=m.fs.purge.inlet)
m.fs.s08 = Arc(source=m.fs.S101.recycle,
               destination=m.fs.P101.inlet)
m.fs.s09 = Arc(source=m.fs.P101.outlet,
               destination=m.fs.M101.vapor_recycle)
m.fs.s10a = Arc(source=m.fs.F101.liq_outlet,
                destination=m.fs.translator.inlet)
m.fs.s10b = Arc(source=m.fs.translator.outlet,
                destination=m.fs.H102.inlet)
# m.fs.s11 = Arc(source=m.fs.H102.outlet,
#                 destination=m.fs.C101.feed)
# m.fs.s12 = Arc(source=m.fs.C101.condenser.distillate,
#                 destination=m.fs.benzene_product.inlet)
# m.fs.s13 = Arc(source=m.fs.C101.reboiler.bottoms,
#                 destination=m.fs.toluene_product.inlet)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

# -----------------------------------------------------------------------------
m.fs.R101.conversion = pyo.Var(initialize=0.75, bounds=(0, 1))

m.fs.R101.conv_constraint = pyo.Constraint(
    expr=m.fs.R101.conversion*m.fs.R101.inlet.
    flow_mol_phase_comp[0, "Vap", "toluene"] ==
    (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
     m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

# # Add operating cost
# m.fs.cooling_cost = pyo.Expression(
#     expr=0.212e-7*(-m.fs.F101.heat_duty[0]) +
#     0.212e-7*(-m.fs.R101.heat_duty[0]) +
#     0.212e-7*(-m.fs.C101.condenser.heat_duty[0]))
# m.fs.heating_cost = pyo.Expression(
#     expr=2.2e-7*m.fs.H101.heat_duty[0] +
#     1.2e-7*m.fs.H102.heat_duty[0] +
#     1.9e-7*m.fs.C101.reboiler.heat_duty[0])
# m.fs.operating_cost = pyo.Expression(
#     expr=(3600*24*365*(m.fs.heating_cost + m.fs.cooling_cost)))

# -----------------------------------------------------------------------------
from idaes.core.util.model_statistics import degrees_of_freedom

print(degrees_of_freedom(m))

m.fs.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
m.fs.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
m.fs.toluene_feed.temperature.fix(303.2)
m.fs.toluene_feed.pressure.fix(350000)

m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
m.fs.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
m.fs.hydrogen_feed.temperature.fix(303.2)
m.fs.hydrogen_feed.pressure.fix(350000)

# -----------------------------------------------------------------------------
m.fs.H101.outlet.temperature.fix(600)

m.fs.R101.conversion.fix(0.75)
m.fs.R101.heat_duty.fix(0)

m.fs.F101.vap_outlet.temperature.fix(325.0)
m.fs.F101.deltaP.fix(0)

m.fs.S101.split_fraction[0, "purge"].fix(0.2)
m.fs.P101.outlet.pressure.fix(350000)

# Set heater outlet conditions (set outlet 2 phase)
m.fs.H102.outlet.temperature.fix(375)
m.fs.H102.deltaP.fix(-200000)

# distillation level inputs
# m.fs.C101.condenser.reflux_ratio.fix(0.5)
# m.fs.C101.condenser.condenser_pressure.fix(150000)
# m.fs.C101.reboiler.boilup_ratio.fix(0.5)

print(degrees_of_freedom(m))

# -----------------------------------------------------------------------------
seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"

seq.options.iterLim = 5

# Using the SD tool
G = seq.create_graph(m)
heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
order = seq.calculation_order(G)

for o in heuristic_tear_set:
    print(o.name)
for o in order:
    print(o[0].name)

# -----------------------------------------------------------------------------
tear_guesses = {
        "flow_mol_phase_comp": {
                (0, "Vap", "benzene"): 1e-5,
                (0, "Vap", "toluene"): 1e-5,
                (0, "Vap", "hydrogen"): 0.30,
                (0, "Vap", "methane"): 0.02,
                (0, "Liq", "benzene"): 1e-5,
                (0, "Liq", "toluene"): 0.30,
                (0, "Liq", "hydrogen"): 1e-5,
                (0, "Liq", "methane"): 1e-5},
        "temperature": {0: 303},
        "pressure": {0: 350000}}

# Pass the tear_guess to the SD tool
seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

def function(unit):
    unit.initialize()

seq.run(m, function)

solver = pyo.SolverFactory('ipopt')
solver.solve(m, tee=True)

# -----------------------------------------------------------------------------
m.fs.C101 = TrayColumn(default={
    "property_package": m.fs.bt_params,
    "number_of_trays": 10,
    "feed_tray_location": 5,
    "condenser_type": CondenserType.totalCondenser,
    "condenser_temperature_spec": TemperatureSpec.atBubblePoint,
    "has_heat_transfer": False,
    "has_pressure_change": False})

# m.fs.benzene_product = Product(default={
#     "property_package": m.fs.bt_params})
m.fs.toluene_product = Product(default={
    "property_package": m.fs.bt_params})

m.fs.s11 = Arc(source=m.fs.H102.outlet,
                destination=m.fs.C101.tray[5].feed)
# m.fs.s12 = Arc(source=m.fs.C101.condenser.distillate,
#                 destination=m.fs.benzene_product.inlet)
# m.fs.s13 = Arc(source=m.fs.C101.reboiler.bottoms,
#                 destination=m.fs.toluene_product.inlet)

pyo.TransformationFactory("network.expand_arcs").apply_to(m)

# Add operating cost
m.fs.cooling_cost = pyo.Expression(
    expr=0.212e-7*(-m.fs.F101.heat_duty[0]) +
    0.212e-7*(-m.fs.R101.heat_duty[0]) +
    0.212e-7*(-m.fs.C101.condenser.heat_duty[0]))
m.fs.heating_cost = pyo.Expression(
    expr=2.2e-7*m.fs.H101.heat_duty[0] +
    1.2e-7*m.fs.H102.heat_duty[0] +
    1.9e-7*m.fs.C101.reboiler.heat_duty[0])
m.fs.operating_cost = pyo.Expression(
    expr=(3600*24*365*(m.fs.heating_cost + m.fs.cooling_cost)))

# distillation level inputs
m.fs.C101.condenser.reflux_ratio.fix(0.5)
m.fs.C101.condenser.condenser_pressure.fix(150000)
m.fs.C101.reboiler.boilup_ratio.fix(0.5)

from idaes.core.util.initialization import propagate_state
propagate_state(m.fs.s11)
m.fs.C101.initialize()
# propagate_state(m.fs.s12)
# propagate_state(m.fs.s13)
# m.fs.benzene_product.initialize()
# m.fs.toluene_product.initialize()

# -----------------------------------------------------------------------------
m.fs.objective = pyo.Objective(expr=m.fs.operating_cost)

m.fs.overhead_loss = pyo.Constraint(
    expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] <=
    0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"])
m.fs.product_flow = pyo.Constraint(
    expr=m.fs.C101.condenser.distillate.flow_mol[0] >=
    0.15)
m.fs.product_purity = pyo.Constraint(
    expr=m.fs.C101.condenser.distillate.mole_frac_comp[0, "benzene"] >=
    0.95)

# -----------------------------------------------------------------------------
# Unfix degrees of freedom
m.fs.H101.outlet.temperature.unfix()
m.fs.R101.heat_duty.unfix()
m.fs.F101.vap_outlet.temperature.unfix()
m.fs.H102.outlet.temperature.unfix()
m.fs.C101.condenser.condenser_pressure.unfix()
m.fs.C101.condenser.reflux_ratio.unfix()
m.fs.C101.reboiler.boilup_ratio.unfix()

# Set bounds
m.fs.H101.outlet.temperature[0].setlb(500)
m.fs.H101.outlet.temperature[0].setub(600)

m.fs.R101.outlet.temperature[0].setlb(600)
m.fs.R101.outlet.temperature[0].setub(800)

m.fs.F101.vap_outlet.temperature[0].setlb(298.0)
m.fs.F101.vap_outlet.temperature[0].setub(450.0)

m.fs.H102.outlet.temperature[0].setlb(350)
m.fs.H102.outlet.temperature[0].setub(400)

m.fs.C101.condenser.condenser_pressure.setlb(101325)
m.fs.C101.condenser.condenser_pressure.setub(150000)

m.fs.C101.condenser.reflux_ratio.setlb(0.1)
m.fs.C101.condenser.reflux_ratio.setub(1.5)

m.fs.C101.reboiler.boilup_ratio.setlb(0.5)
m.fs.C101.reboiler.boilup_ratio.setub(1.5)

solver.solve(m, tee=True)
# -----------------------------------------------------------------------------

m.fs.M101.report()
m.fs.H101.report()
m.fs.R101.report()
m.fs.F101.report()
m.fs.S101.report()
m.fs.P101.report()
m.fs.H102.report()
m.fs.C101.condenser.report()
m.fs.C101.reboiler.report()

# -----------------------------------------------------------------------------
