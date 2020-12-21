#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:44:20 2020

@author: andrew
"""

from pyomo.environ import (ConcreteModel,
                           Constraint,
                           Expression,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           Var,
                           TerminationCondition)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.infeasible import log_infeasible_bounds, \
    log_infeasible_constraints
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import (Feed,
                                              PressureChanger,
                                              Mixer,
                                              Separator as Splitter,
                                              Heater,
                                              CSTR,
                                              Flash,
                                              Translator)
from idaes.generic_models.unit_models.distillation import TrayColumn
from idaes.generic_models.unit_models.distillation.condenser \
    import CondenserType, TemperatureSpec

from idaes.generic_models.unit_models.pressure_changer \
    import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state

from idaes.generic_models.properties.activity_coeff_models.\
    BTX_activity_coeff_VLE import BTXParameterBlock
from idaes.generic_models.properties.core.generic.generic_property \
    import GenericParameterBlock
from BTHM_ideal import configuration
import hda_reaction_kinetic as reaction_props

m = ConcreteModel()
m.fs = FlowsheetBlock(default={"dynamic": False})

m.fs.bthm_params = GenericParameterBlock(default=configuration)
m.fs.bt_params = BTXParameterBlock(default={
    "valid_phase": ('Liq', 'Vap'),
    "activity_coeff_model": "Ideal"})
m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.bthm_params})

m.fs.toluene_feed = Feed(default={"property_package": m.fs.bthm_params})
m.fs.hydrogen_feed = Feed(default={"property_package": m.fs.bthm_params})

m.fs.M101 = Mixer(default={
    "property_package": m.fs.bthm_params,
    "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})

m.fs.H101 = Heater(default={"property_package": m.fs.bthm_params,
                            "has_pressure_change": False,
                            "has_phase_equilibrium": True})

m.fs.R101 = CSTR(
            default={"property_package": m.fs.bthm_params,
                     "reaction_package": m.fs.reaction_params,
                     "has_heat_of_reaction": False,
                     "has_heat_transfer": True,
                     "has_pressure_change": False})

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

# -----------------------------------------------------------------------------
# Translator and H102
m.fs.translator = Translator(default={
                             "inlet_property_package": m.fs.bthm_params,
                             "outlet_property_package":
                                 m.fs.bt_params})

m.fs.H102 = Heater(default={"property_package": m.fs.bt_params,
                            "has_pressure_change": True,
                            "has_phase_equilibrium": True})

# Translator constraints linking outlet state variables to inlet
# state variables
m.fs.translator.eq_total_flow = Constraint(
    expr=m.fs.translator.outlet.flow_mol[0] ==
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"])

m.fs.translator.eq_temperature = Constraint(
    expr=m.fs.translator.outlet.temperature[0] ==
    m.fs.translator.inlet.temperature[0])

m.fs.translator.eq_pressure = Constraint(
    expr=m.fs.translator.outlet.pressure[0] ==
    m.fs.translator.inlet.pressure[0])

m.fs.translator.eq_mole_frac_benzene = Constraint(
    expr=m.fs.translator.outlet.mole_frac_comp[0, "benzene"] ==
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] /
    (m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
     m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))

m.fs.translator.eq_mole_frac_toluene = Constraint(
    expr=m.fs.translator.outlet.mole_frac_comp[0, "toluene"] ==
    m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"] /
    (m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "benzene"] +
     m.fs.translator.inlet.flow_mol_phase_comp[0, "Liq", "toluene"]))
# -----------------------------------------------------------------------------

# Arcs
m.fs.s01 = Arc(source=m.fs.hydrogen_feed.outlet,
               destination=m.fs.M101.hydrogen_feed)
m.fs.s02 = Arc(source=m.fs.toluene_feed.outlet,
               destination=m.fs.M101.toluene_feed)
m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.P101.inlet)
m.fs.s09 = Arc(source=m.fs.P101.outlet,
               destination=m.fs.M101.vapor_recycle)
m.fs.s10a = Arc(source=m.fs.F101.liq_outlet,
                destination=m.fs.translator.inlet)
m.fs.s10b = Arc(source=m.fs.translator.outlet,
                destination=m.fs.H102.inlet)

TransformationFactory("network.expand_arcs").apply_to(m)

print(degrees_of_freedom(m))

t_feed = m.fs.toluene_feed
t_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
t_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
t_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
t_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
t_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
t_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
t_feed.temperature.fix(303.2)
t_feed.pressure.fix(350000)

h_feed = m.fs.hydrogen_feed
h_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
h_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
h_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
h_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
h_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
h_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
h_feed.temperature.fix(303.2)
h_feed.pressure.fix(350000)

m.fs.H101.outlet.temperature.fix(600)

m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

m.fs.R101.conv_constraint = Constraint(
    expr=m.fs.R101.conversion*m.fs.R101.inlet.
    flow_mol_phase_comp[0, "Vap", "toluene"] ==
    (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
     m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

m.fs.R101.conversion.fix(0.75)
m.fs.R101.heat_duty.fix(0)

m.fs.F101.vap_outlet.temperature.fix(325.0)
m.fs.F101.deltaP.fix(0)

m.fs.S101.split_fraction[0, "purge"].fix(0.2)
m.fs.P101.outlet.pressure.fix(350000)

m.fs.H102.outlet.temperature.fix(375)
m.fs.H102.deltaP.fix(-200000)

print(degrees_of_freedom(m))

seq = SequentialDecomposition()
seq.options.select_tear_method = "heuristic"
seq.options.tear_method = "Wegstein"
seq.options.iterLim = 5

# Using the SD tool
G = seq.create_graph(m)
edges = G.copy().edges()
for e in edges:
    if all(_.parent_block() is not m.fs for _ in e):
        G.remove_edge(*e)
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
    unit.initialize()


seq.run(m, function)

solver = SolverFactory('ipopt')
res = solver.solve(m, tee=True)

# -----------------------------------------------------------------------------
# Distillation
m.fs.C101 = TrayColumn(default={
    "number_of_trays": 10,
    "feed_tray_location": 5,
    "condenser_type": CondenserType.totalCondenser,
    "condenser_temperature_spec": TemperatureSpec.atBubblePoint,
    "property_package": m.fs.bt_params,
    "has_heat_transfer": False,
    "has_pressure_change": False})

m.fs.s11 = Arc(
    source=m.fs.H102.outlet, destination=m.fs.C101.feed)

# distillation level inputs
m.fs.C101.condenser.reflux_ratio.fix(0.5)
m.fs.C101.condenser.condenser_pressure.fix(150000)
m.fs.C101.reboiler.boilup_ratio.fix(0.5)

TransformationFactory("network.expand_arcs").apply_to(m)

propagate_state(m.fs.s11)
m.fs.C101.initialize()

res = solver.solve(m, tee=True)

assert res.solver.termination_condition == TerminationCondition.optimal

m.fs.H101.report()
m.fs.R101.report()
m.fs.F101.report()
m.fs.H102.report()
m.fs.C101.condenser.report()
m.fs.C101.reboiler.report()

# -----------------------------------------------------------------------------
# optimize
# Add operating cost
m.fs.cooling_cost = Expression(
    expr=0.212e-7 * (-m.fs.F101.heat_duty[0]) +
    0.212e-7 * (-m.fs.R101.heat_duty[0]) +
    0.212e-7 * (-m.fs.C101.condenser.heat_duty[0]))

m.fs.heating_cost = Expression(
    expr=2.2e-7 * m.fs.H101.heat_duty[0] +
    1.2e-7 * m.fs.H102.heat_duty[0] +
    1.9e-7 * m.fs.C101. reboiler.heat_duty[0])

m.fs.operating_cost = Expression(
    expr=(3600 * 24 * 365 * (m.fs.heating_cost + m.fs.cooling_cost)))

m.fs.capital_cost = Expression(
    expr=4e5*m.fs.R101.volume[0])

m.fs.total_cost = Expression(
    expr=m.fs.operating_cost+m.fs.capital_cost)

m.fs.objective = Objective(expr=m.fs.total_cost)

# Change reactor configuration
m.fs.R101.conversion.unfix()
# m.fs.R101.volume.fix()

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
# m.fs.R101.heat_duty.setlb(-10000)
m.fs.R101.heat_duty.setub(0)
m.fs.R101.volume.setlb(0.05)
m.fs.R101.volume.setub(0.5)

m.fs.F101.vap_outlet.temperature[0].setlb(298.0)
m.fs.F101.vap_outlet.temperature[0].setub(450.0)

m.fs.H102.outlet.temperature[0].setlb(350)
m.fs.H102.outlet.temperature[0].setub(400)

m.fs.C101.condenser.condenser_pressure.setlb(101325)
m.fs.C101.condenser.condenser_pressure.setub(150000)

m.fs.C101.condenser.reflux_ratio.setlb(0.1)
m.fs.C101.condenser.reflux_ratio.setub(1.5)

m.fs.C101.reboiler.boilup_ratio.setlb(0.1)
m.fs.C101.reboiler.boilup_ratio.setub(1.5)

m.fs.overhead_loss = Constraint(
    expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] <=
    0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"])
m.fs.product_flow = Constraint(
    expr=m.fs.C101.condenser.distillate.flow_mol[0] >=
    0.15)
m.fs.product_purity = Constraint(
    expr=m.fs.C101.condenser.
    distillate.mole_frac_comp[0, "benzene"] >= 0.95)

# m.fs.product_flow = Constraint(
#     expr=m.fs.H102.outlet.flow_mol[0] >= 0.25)
# m.fs.product_purity = Constraint(
#     expr=m.fs.H102.outlet.mole_frac_comp[0, "benzene"] >= 0.75)

res = solver.solve(m, tee=True)

if not res.solver.termination_condition == TerminationCondition.optimal:
    log_infeasible_constraints(m)

raise Exception(res)
# -----------------------------------------------------------------------------
m.fs.F101.report()
m.fs.R101.report()
m.fs.H102.report()
m.fs.C101.condenser.report()
m.fs.C101.reboiler.report()

m.fs.operating_cost.display()
m.fs.capital_cost.display()
m.fs.total_cost.display()
m.fs.R101.conversion.display()
