# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import os
import sys

import argparse as ap
import numpy as np

from importlib.metadata import version

from sandlerprops.properties import PropertiesDatabase
from sandlermisc.statereporter import StateReporter

from .eos import CubicEOS
from .idealgas import IdealGasEOS
from .soaveredlichkwong import SoaveRedlichKwongEOS
from .pengrobinson import PengRobinsonEOS
from .vanderwaals import VanDerWaalsEOS

banner = r"""
 __                 _ _           
/ _\ __ _ _ __   __| | | ___ _ __ 
\ \ / _` | '_ \ / _` | |/ _ \ '__|
_\ \ (_| | | | | (_| | |  __/ |   
\__/\__,_|_| |_|\__,_|_|\___|_|   
               _     _            
     ___ _   _| |__ (_) ___ ___   
    / __| | | | '_ \| |/ __/ __|  
   | (__| |_| | |_) | | (__\__ \  
    \___|\__,_|_.__/|_|\___|___/  v""" + version("sandlercubics") + """

"""

def bar_to_unit(value_in_bar: float, to_unit: str) -> float:
    """
    Convert pressure from bar to specified unit.
    
    Parameters
    ----------
    value_in_bar : float
        Pressure value in bar.
    to_unit : str
        Target unit ('mpa', 'MPa', 'kpa', 'kPa', 'pa', 'Pa', 'bar', 'atm').

    Returns
    -------
    float
        Converted pressure value.
    """
    conversions = {
        'mpa': 0.1,
        'MPa': 0.1,
        'kpa': 100.0,
        'kPa': 100.0,
        'pa': 1e5,
        'Pa': 1e5,
        'bar': 1.0,
        'atm': 0.986923,
    }
    if to_unit not in conversions:
        raise ValueError(f"Unsupported pressure unit conversion to {to_unit}.")
    return value_in_bar * conversions[to_unit]

def reporters(eos: CubicEOS) -> str:
    """
    Generate state and property reports for the given EOS.
    
    Parameters
    ----------
    eos : CubicEOS
        The equation of state object.

    Returns
    -------
    str
        Formatted state and property reports.
    """

    result = StateReporter({})
    
    result.add_property('T', eos.T, 'K', fstring="{: .2f}")
    pu = eos.R._capitalizations.get(eos.pressure_unit, eos.pressure_unit)
    result.add_property('P', eos.P, pu, fstring="{: .2f}")
    result.add_property('Z', ', '.join([f"{z: 4g}" for z in eos.Z]), '', fstring=None)
    vu = eos.R._capitalizations.get(eos.volume_unit, eos.volume_unit)
    result.add_property('v', ', '.join([f"{vol: 6g}" for vol in eos.v]), f'{vu}/mol', fstring=None)
    result.add_property('h', ', '.join([f"{h: 6g}" for h in eos.h]), 'J/mol', fstring=None)
    result.add_property('s', ', '.join([f"{s: 6g}" for s in eos.s]), 'J/mol-K', fstring=None)
    result.add_property('hdep', ', '.join([f"{hdep: 6g}" for hdep in eos.h_departure]), 'J/mol', fstring=None)
    result.add_property('sdep', ', '.join([f"{sdep: 6g}" for sdep in eos.s_departure]), 'J/mol-K', fstring=None)
    if eos.T < eos.Tc:
        pu = eos.R._capitalizations.get(eos.pressure_unit, eos.pressure_unit)
        if eos.Pvap is not np.nan:
            result.add_property(f'Pvap({eos.T:.2f} K)', eos.Pvap, pu, fstring="{:.2f}")
        if eos.Hvap is not np.nan:
            result.add_property(f'Hvap({eos.T:.2f} K)', eos.Hvap, 'J/mol', fstring="{:.2f}")
        if eos.Svap is not np.nan:
            result.add_property(f'Svap({eos.T:.2f} K)', eos.Svap, 'J/mol-K', fstring="{:.4f}")
        if eos.P < eos.Pc:
            pu = eos.R._capitalizations.get(eos.pressure_unit, eos.pressure_unit)
            if eos.Tsat is not np.nan:
                result.add_property(f'Tsat({eos.P:.2f} {pu})', eos.Tsat, 'K', fstring="{:.2f}")
        if eos.x is not None:
            result.add_property('Vapor fraction x', eos.x, '', fstring="{:4g}")
    prop = StateReporter({})
    prop.add_property('Tc', eos.Tc, 'K', fstring="{:.2f}")
    pu = eos.R._capitalizations.get(eos.pressure_unit, eos.pressure_unit)
    prop.add_property('Pc', eos.Pc, pu, fstring="{:.2f}")
    prop.add_property('omega', eos.omega, '', fstring="{:.3f}")

    prop.add_property('Tref', eos.Tref, 'K', fstring="{:.2f}")
    prop.add_property('Pref', eos.Pref_MPa, 'MPa', fstring="{:.2f}")
    # print(eos.Cp)
    prop.pack_Cp(eos.Cp, fmts=["{:.2f}", "{:.3e}", "{:.3e}", "{:.3e}"])
    return result.report(),  prop.report()

def state(args):
    """
    Calculate and report the state for a single condition using the specified EOS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    pressure_unit = args.pressure_unit
    volume_unit = args.volume_unit
    if args.n is not None:
        db = PropertiesDatabase()
        component = db.get_compound(args.n)
        if component is None:
            print(f"Component '{args.n}' not found in database.")
            return
    match args.eos_type:
        case 'ideal':
            eos = IdealGasEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
        case 'vdw':
            eos = VanDerWaalsEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
        case 'pr':
            eos = PengRobinsonEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
        case 'srk':
            eos = SoaveRedlichKwongEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
    if component is not None:
        eos.set_compound(component)
    if args.Tc is not None:
        eos.Tc = args.Tc
    if args.Pc is not None:
        eos.Pc = args.Pc
    if args.w is not None:
        eos.omega = args.w
    if args.Cp is not None:
        eos.Cp = args.Cp
    if args.eos_type == 'ideal':
        eos.solve(T=args.T, P=args.P, v=args.v)
    else:
        eos.solve(T=args.T, P=args.P, v=args.v, h=args.h, s=args.s, u=args.u)
    state_report, prop_report = reporters(eos)
    print(f'State report for {component.Name} using {eos.description}:')
    print(state_report)
    if prop_report and args.show_props:
        print("\nConstants used for calculations:")
        print(prop_report)

def delta(args):
    """
    Calculate and report property differences between two states using the specified EOS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    pressure_unit = args.pressure_unit
    volume_unit = args.volume_unit
    if args.n is not None:
        db = PropertiesDatabase()
        component = db.get_compound(args.n) # critical pressures are in bars!
        if component is None:
            print(f"Component '{args.n}' not found in database.")
            return
    match args.eos_type:
        case 'ideal':
            eos1 = IdealGasEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
            eos2 = IdealGasEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
        case 'vdw':
            eos1 = VanDerWaalsEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
            eos2 = VanDerWaalsEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
        case 'pr':
            eos1 = PengRobinsonEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
            eos2 = PengRobinsonEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
        case 'srk':
            eos1 = SoaveRedlichKwongEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
            eos2 = SoaveRedlichKwongEOS(pressure_unit=pressure_unit, volume_unit=volume_unit)
    if component is not None:
        eos1.set_compound(component)
        eos2.set_compound(component)
    if args.Tc is not None:
        eos1.Tc = args.Tc
        eos2.Tc = args.Tc
    if args.Pc is not None:
        eos1.Pc = args.Pc
        eos2.Pc = args.Pc
    if args.w is not None:
        eos1.omega = args.w
        eos2.omega = args.w
    if args.Cp is not None:
        eos1.Cp = args.Cp
        eos2.Cp = args.Cp
    if args.eos_type == 'ideal':
        eos1.solve(T=args.T1, P=args.P1, v=args.v1)
        eos2.solve(T=args.T2, P=args.P2, v=args.v2)
    else:
        eos1.solve(T=args.T1, P=args.P1, v=args.v1, h=args.h1, s=args.s1, u=args.u1)
        eos2.solve(T=args.T2, P=args.P2, v=args.v2, h=args.h2, s=args.s2, u=args.u2)

    delta_State = StateReporter({})
    delta_H = eos2.h - eos1.h
    delta_S = eos2.s - eos1.s
    delta_U = eos2.u - eos1.u
    delta_State.add_property('Δh', ', '.join(f'{x: 7g}' for x in delta_H), 'J/mol', fstring=None)
    delta_State.add_property('Δs', ', '.join(f'{x: 7g}' for x in delta_S), 'J/mol-K', fstring=None)
    delta_State.add_property('Δu', ', '.join(f'{x: 7g}' for x in delta_U), 'J/mol', fstring=None)
    state_1, _ = reporters(eos1)
    state_2, consts = reporters(eos2)
    print(f"State-change calculations for {component.Name} using {eos1.description}:")
    if args.show_states:
        print()
        two_states = ["State 1:                       State 2:"]
        for line1, line2 in zip(state_1.splitlines(), state_2.splitlines()):
            two_states.append(f"{line1:<26s}     {line2}")
        print("\n".join(two_states))
        print()
        print("Property changes:")
    print(delta_State.report())
    if args.show_props:
        print("\nConstants used for calculations:")
        print(consts)
    
def cli():
    """
    Command-line interface for sandlercubics package.
    """

    subcommands = {
        'state': dict(
            func = state,
            help = 'work with a cubic equation of state for a single state'
        ),
        'delta': dict(
            func = delta,
            help = 'work with property differences between two states (not implemented yet)'
        ),
    }
    parser = ap.ArgumentParser(
        prog='sandlercubics',
        description="Sandlercubics: A collection of computational tools using cubic equations of state based on Chemical, Biochemical, and Engineering Thermodynamics (5th edition) by Stan Sandler",
        epilog="(c) 2026, Cameron F. Abrams <cfa22@drexel.edu>"
    )
    parser.add_argument(
        '-b',
        '--banner',
        default=False,
        action=ap.BooleanOptionalAction,
        help='toggle banner message'
    )
    parser.add_argument(
        '-v',
        '--version',
        action='version',
        version=f'sandlercubics version {version("sandlercubics")}',
        help='show program version and exit'
    )
    subparsers = parser.add_subparsers(
        title="subcommands",
        dest="command",
        metavar="<command>",
        required=True,
    )
    command_parsers={}
    for k, specs in subcommands.items():
        command_parsers[k] = subparsers.add_parser(
            k,
            help=specs['help'],
            add_help=False,
            formatter_class=ap.RawDescriptionHelpFormatter
        )
        command_parsers[k].set_defaults(func=specs['func'])
        command_parsers[k].add_argument(
            '--help',
            action='help',
            help=specs['help']
        )

    options = [
        ('eos', 'eos_type', 'type of cubic equation of state to use', str, 'vdw', ['ideal', 'vdw', 'pr', 'srk']),
        ('pu', 'pressure_unit', 'pressure unit (mpa, bar, atm)', str, 'mpa', ['mpa', 'bar', 'atm']),
        ('vu', 'volume_unit', 'volume unit (m3, l, cm3)', str, 'm3', ['m3', 'l', 'cm3']),
    ]
    for short, long, desc, typ, default, choices in options:
        command_parsers['state'].add_argument(
            f'-{short}',
            f'--{long}',
            dest=long,
            type=typ,
            default=default,
            choices=choices,
            help=desc
        )
        command_parsers['delta'].add_argument(
            f'-{short}',
            f'--{long}',
            dest=long,
            type=typ,
            default=default,
            choices=choices,
            help=desc
        )

    crit_args = [
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, False),
        ('Pc', 'critical_pressure', 'critical pressure (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
        ('w', 'acentric_factor', 'acentric factor omega (if component not specified)', float, False),
    ]

    cp_args = [
        ('CpA', 'CpA', 'heat capacity coefficient A (J/mol-K) (if component not specified)', float, False),
        ('CpB', 'CpB', 'heat capacity coefficient B (J/mol-K^2) (if component not specified)', float, False),
        ('CpC', 'CpC', 'heat capacity coefficient C (J/mol-K^3) (if component not specified)', float, False),
        ('CpD', 'CpD', 'heat capacity coefficient D (J/mol-K^4) (if component not specified)', float, False),
    ]

    state_args = [
        ('P', 'pressure', 'pressure', float, False),
        ('T', 'temperature', 'temperature in K (always in K)', float, False),
        ('v', 'molar_volume', 'molar volume', float, False),
        ('h', 'enthalpy', 'molar enthalpy in J/mol', float, False),
        ('s', 'entropy', 'molar entropy in J/mol-K', float, False),
        ('u', 'internal_energy', 'molar internal energy in J/mol', float, False),
    ]
    for prop, long_arg, explanation, arg_type, required in state_args + crit_args:
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    
    command_parsers['state'].add_argument(
        '--Cp',
        nargs=4,
        type=float,
        metavar=('CpA', 'CpB', 'CpC', 'CpD'),
        help='heat capacity polynomial coefficients A, B, C, D (J/mol-K, J/mol-K^2, J/mol-K^3, J/mol-K^4) (if component not specified)',
        default=None
    )
    command_parsers['state'].add_argument(
        '--show-props',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show all critical properties and Cp coefficients used'
    )

    delta_args = [
        ('P1', 'pressure1', 'pressure of state 1', float, False),
        ('T1', 'temperature1', 'temperature of state 1 in K (always in K)', float, False),
        ('v1', 'molar_volume1', 'molar volume of state 1', float, False),
        ('h1', 'enthalpy1', 'molar enthalpy of state 1 in J/mol', float, False),
        ('s1', 'entropy1', 'molar entropy of state 1 in J/mol-K', float, False),
        ('u1', 'internal_energy1', 'molar internal energy of state 1 in J/mol', float, False),
        ('P2', 'pressure2', 'pressure of state 2', float, False),
        ('T2', 'temperature2', 'temperature of state 2 in K (always in K)', float, False),
        ('v2', 'molar_volume2', 'molar volume of state 2', float, False),
        ('h2', 'enthalpy2', 'molar enthalpy of state 2 in J/mol', float, False),
        ('s2', 'entropy2', 'molar entropy of state 2 in J/mol-K', float, False),
        ('u2', 'internal_energy2', 'molar internal energy of state 2 in J/mol', float, False),
    ]
    for prop, long_arg, explanation, arg_type, required in delta_args + crit_args:
        command_parsers['delta'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )

    command_parsers['delta'].add_argument(
        '--show-states',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show the full states for state 1 and state 2'
    )
    command_parsers['delta'].add_argument(
        '--show-props',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show all critical properties and Cp coefficients used'
    )
    command_parsers['delta'].add_argument(
        '--Cp',
        nargs=4,
        type=float,
        metavar=('CpA', 'CpB', 'CpC', 'CpD'),
        help='heat capacity polynomial coefficients A, B, C, D (J/mol-K, J/mol-K^2, J/mol-K^3, J/mol-K^4) (if component not specified)',
        default=None
    )

    args = parser.parse_args()

    if args.func == state:
        nprops = 0
        for prop, _, _, _, _ in state_args:
            if hasattr(args, prop) and getattr(args, prop) is not None:
                nprops += 1
        if nprops > 2:
            parser.error('At most two of P, T, x, v, u, h, and s may be specified for "state" subcommand')

    if args.banner:
        print(banner)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')
    if args.banner:
        print('Thanks for using sandlercubics!')