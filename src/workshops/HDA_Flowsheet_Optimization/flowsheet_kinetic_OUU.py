#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  2 13:44:20 2020

@author: andrew
"""

import argparse
import pytest
from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value,
                           TerminationCondition)
from pyomo.network import Arc, SequentialDecomposition
from pyomo.util.infeasible import log_infeasible_constraints, \
    log_infeasible_bounds
from idaes.core import FlowsheetBlock
from idaes.generic_models.unit_models import (PressureChanger,
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
from idaes.core.util.testing import get_default_solver
# Import idaes logger to set output levels
import idaes.logger as idaeslog

from idaes.generic_models.properties.core.generic.generic_property \
    import GenericParameterBlock
from idaes.generic_models.properties.activity_coeff_models.\
    BTX_activity_coeff_VLE import BTXParameterBlock
from BTHM_ideal import configuration
# import hda_reaction as reaction_props
import hda_reaction_kinetic as reaction_props
from hda_ideal_VLE import HDAParameterBlock


def get_init_model(outlvl=idaeslog.WARNING):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    m.fs.thermo_params = HDAParameterBlock()
    m.fs.bt_properties = BTXParameterBlock(default={
                                           "valid_phase":
                                           ('Liq', 'Vap'),
                                           "activity_coeff_model":
                                           "Ideal"})

    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})
    m.fs.translator = Translator(default={
                                 "inlet_property_package": m.fs.thermo_params,
                                 "outlet_property_package":
                                     m.fs.bt_properties})

    m.fs.pre_heater = Heater(default={"property_package": m.fs.bt_properties,
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


    m.fs.M101 = Mixer(default={
        "property_package": m.fs.thermo_params,
        "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]})

    m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                                "has_pressure_change": False,
                                "has_phase_equilibrium": True})

    m.fs.R101 = CSTR(
                default={"property_package": m.fs.thermo_params,
                         "reaction_package": m.fs.reaction_params,
                         "has_heat_of_reaction": True,
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
                "thermodynamic_assumption":
                ThermodynamicAssumption.isothermal})

    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet,
                   destination=m.fs.M101.vapor_recycle)

    m.fs.s10 = Arc(source=m.fs.F101.liq_outlet,
                   destination=m.fs.translator.inlet)
    m.fs.s11 = Arc(source=m.fs.translator.outlet,
                   destination=m.fs.pre_heater.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.temperature.fix(303.2)
    m.fs.M101.toluene_feed.pressure.fix(350000)

    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
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

    m.fs.F101.vap_outlet.temperature.fix(325.0)
    m.fs.F101.deltaP.fix(0)

    m.fs.S101.split_fraction[0, "purge"].fix(0.2)
    m.fs.C101.outlet.pressure.fix(350000)

    m.fs.pre_heater.outlet.temperature.fix(375)
    m.fs.pre_heater.deltaP.fix(-200000)

    assert degrees_of_freedom(m) == 0

    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 3

    # Using the SD tool
    G = seq.create_graph(m)
    edges = G.copy().edges()
    for e in edges:
        if all(_.parent_block() is not m.fs for _ in e):
            G.remove_edge(*e)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)

    if not outlvl == idaeslog.WARNING:
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
                    (0, "Liq", "toluene"): 0.30,
                    (0, "Liq", "hydrogen"): 1e-5,
                    (0, "Liq", "methane"): 1e-5},
            "temperature": {0: 303},
            "pressure": {0: 350000}}

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

    def function(unit):
        unit.initialize(outlvl=outlvl)

    seq.run(m, function)

    solver = get_default_solver()
    solver.solve(m, tee=False)

    m.fs.distillation = TrayColumn(default={
                                   "number_of_trays": 10,
                                   "feed_tray_location": 5,
                                   "condenser_type":
                                       CondenserType.totalCondenser,
                                   "condenser_temperature_spec":
                                       TemperatureSpec.atBubblePoint,
                                       "property_package": m.fs.bt_properties,
                                       "has_heat_transfer": False,
                                       "has_pressure_change": False})

    m.fs.s12 = Arc(source=m.fs.pre_heater.outlet,
                   destination=m.fs.distillation.feed)

    # distillation level inputs
    m.fs.distillation.condenser.reflux_ratio.fix(0.5)
    m.fs.distillation.condenser.condenser_pressure.fix(150000)

    m.fs.distillation.reboiler.boilup_ratio.fix(0.5)

    TransformationFactory("network.expand_arcs").apply_to(m)

    propagate_state(m.fs.s12)

    m.fs.distillation.initialize(outlvl=outlvl)

    assert degrees_of_freedom(m) == 0
    solver.solve(m, tee=False)

    return m
    # assert value(m.fs.H101.outlet.temperature[0]) == \
    #     pytest.approx(600, rel=1e-3)
    # if args.prop == "oldprop":
    #     assert value(m.fs.R101.outlet.temperature[0]) == \
    #         pytest.approx(771.84, rel=1e-3)
    # else:
    #     assert value(m.fs.R101.outlet.temperature[0]) == \
    #         pytest.approx(760.58, rel=1e-3)
    # assert value(m.fs.F101.vap_outlet.temperature[0]) == \
    #     pytest.approx(325, rel=1e-3)
    # assert value(m.fs.pre_heater.outlet.temperature[0]) == \
    #     pytest.approx(375, rel=1e-3)
    # assert value(m.fs.distillation.condenser.condenser_pressure[0]) == \
    #     pytest.approx(150000, rel=1e-2)
    # assert value(m.fs.distillation.condenser.reflux_ratio) == \
    #     pytest.approx(0.5, rel=1e-3)
    # assert value(m.fs.distillation.reboiler.boilup_ratio) == \
    #     pytest.approx(0.5, rel=1e-3)


def get_opt_model(m):
    # m - an initialized model

    # Add operating cost
    # Removed reactor cooling duty from expression, as it will now remain fixed
    m.fs.cooling_cost = Expression(expr=0.25e-7 * (-m.fs.F101.heat_duty[0]) +
                                   0.2e-7 * (-m.fs.distillation.
                                             condenser.heat_duty[0]))
    m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0] +
                                   1.2e-7 * m.fs.pre_heater.heat_duty[0] +
                                   1.9e-7 * m.fs.distillation.
                                   reboiler.heat_duty[0])
    m.fs.operating_cost = Expression(expr=(3600 * 24 * 365 *
                                           (m.fs.heating_cost +
                                            m.fs.cooling_cost)))

    # Add expresion for capital cost
    # Random cost coefficient
    m.fs.capital_cost = Expression(expr=1e5*m.fs.R101.volume[0])

    # Added capital cost to objective
    m.fs.objective = Objective(expr=m.fs.operating_cost+m.fs.capital_cost)

    # Unfix degrees of freedom
    m.fs.H101.outlet.temperature.unfix()
    # Do not unfix reactor heat duty - trying to optimize this causes issues
    m.fs.R101.conversion.unfix()  # Unfix conversion constraint
    m.fs.F101.vap_outlet.temperature.unfix()
    m.fs.pre_heater.outlet.temperature.unfix()
    m.fs.distillation.condenser.condenser_pressure.unfix()
    m.fs.distillation.condenser.reflux_ratio.unfix()
    m.fs.distillation.reboiler.boilup_ratio.unfix()

    # Set bounds
    m.fs.H101.outlet.temperature[0].setlb(500)
    m.fs.H101.outlet.temperature[0].setub(600)

    m.fs.R101.outlet.temperature[0].setlb(600)
    m.fs.R101.outlet.temperature[0].setub(800)
    # New bounds
    m.fs.R101.volume[0].setlb(0)
    # Adding an upper bound on volume causes solver failures for some reason.

    m.fs.F101.vap_outlet.temperature[0].setlb(298)
    m.fs.F101.vap_outlet.temperature[0].setub(450.0)

    m.fs.pre_heater.outlet.temperature[0].setlb(350)
    m.fs.pre_heater.outlet.temperature[0].setub(400)

    m.fs.distillation.condenser.condenser_pressure.setlb(101325)
    m.fs.distillation.condenser.condenser_pressure.setub(150000)

    m.fs.distillation.condenser.reflux_ratio.setlb(0.1)
    m.fs.distillation.condenser.reflux_ratio.setub(1.5)

    m.fs.distillation.reboiler.boilup_ratio.setlb(0.1)
    m.fs.distillation.reboiler.boilup_ratio.setub(2)

    m.fs.overhead_loss = Constraint(
        expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] <=
        0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"])
    m.fs.product_flow = Constraint(
        expr=m.fs.distillation.condenser.distillate.flow_mol[0] >=
        0.18)
    m.fs.product_purity = Constraint(
        expr=m.fs.distillation.condenser.
        distillate.mole_frac_comp[0, "benzene"] >= 0.99)

    return m


def report_optimal(m):

    print("Optimal operating variables:")
    print("optimal H101 temperature is ",
          value(m.fs.H101.outlet.temperature[0]))
    print("optimal R101 temperature is ",
          value(m.fs.R101.outlet.temperature[0]))
    print("optimal F101 outlet temperature is ",
          value(m.fs.F101.vap_outlet.temperature[0]))
    print("optimal pre_heater outlet temperature is ",
          value(m.fs.pre_heater.outlet.temperature[0]))
    print("optimal condenser pressure is ",
          value(m.fs.distillation.condenser.condenser_pressure[0]))
    print("optimal reflux ratio is ",
          value(m.fs.distillation.condenser.reflux_ratio))
    print("optimal boil up ratio is ",
          value(m.fs.distillation.reboiler.boilup_ratio))

    # print("Column Report")
    m.fs.H101.report()
    m.fs.R101.report()
    m.fs.R101.conversion.display()
    m.fs.F101.report()
    m.fs.pre_heater.report()
    m.fs.distillation.condenser.report()
    m.fs.distillation.reboiler.report()

    # Check expected optimal values
    # assert value(m.fs.H101.outlet.temperature[0]) == \
    #     pytest.approx(500, rel=1e-3)
    # assert value(m.fs.R101.outlet.temperature[0]) == \
    #     pytest.approx(696.11, rel=1e-3)
    # assert value(m.fs.F101.vap_outlet.temperature[0]) == \
    #     pytest.approx(301.87, rel=1e-3)
    # assert value(m.fs.pre_heater.outlet.temperature[0]) == \
    #     pytest.approx(373.63, rel=1e-3)
    # assert value(m.fs.distillation.condenser.condenser_pressure[0]) == \
    #     pytest.approx(101325, rel=1e-2)
    # assert value(m.fs.distillation.condenser.reflux_ratio) == \
    #     pytest.approx(0.777, rel=1e-3)
    # assert value(m.fs.distillation.reboiler.boilup_ratio) == \
    #     pytest.approx(0.984, rel=1e-3)


if __name__ == "__main__":
    m = get_init_model()
    m = get_opt_model(m)
    solver = get_default_solver()
    res = solver.solve(m, tee=True)
    report_optimal(m)
