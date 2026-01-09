
import os
import sys

import argparse as ap

from importlib.metadata import version

from sandlerprops.properties import PropertiesDatabase
from sandlermisc.statereporter import StateReporter

from .eos import CubicEOS, IdealGasEOS, SoaveRedlichKwongEOS, PengRobinsonEOS, VanDerWaalsEOS

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

def apply_schema(args: ap.Namespace) -> ap.Namespace:
    """
    Apply any necessary transformations to the parsed arguments.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    argparse.Namespace
        Transformed command-line arguments.
    """
    # Currently, no transformations are needed.
    if args.command == 'state':
        # must specify component, and two of T, P, v
        assert args.n is not None, "Component name must be specified."
        specified = sum([args.T is not None, args.P is not None, args.v is not None])
        assert specified == 2, "Please specify exactly two of temperature, pressure, or molar volume."
    elif args.command == 'delta':
        # must specify component, and T1, P1, T2, P2
        assert args.n is not None, "Component name must be specified."
        specified1 = sum([args.T1 is not None, args.P1 is not None, args.v1 is not None])
        specified2 = sum([args.T2 is not None, args.P2 is not None, args.v2 is not None])
        assert specified1 == 2, "For state 1, please specify exactly two of temperature, pressure, or molar volume."
        assert specified2 == 2, "For state 2, please specify exactly two of temperature, pressure, or molar volume."

    return args

def bar_to_unit(value_in_bar: float, to_unit: str) -> float:
    """
    Convert pressure from bar to specified unit.
    
    Parameters
    ----------
    value_in_bar : float
        Pressure value in bar.
    to_unit : str
        Target unit ('mpa', 'bar', 'atm', 'pa').

    Returns
    -------
    float
        Converted pressure value.
    """
    conversions = {
        'mpa': 0.1,
        'bar': 1.0,
        'atm': 0.986923,
        'pa': 1e5,
    }
    if to_unit not in conversions:
        raise ValueError(f"Unsupported pressure unit conversion to {to_unit}.")
    return value_in_bar * conversions[to_unit]

def reporters(eos: CubicEOS, eos_type: str, Cp: float | list[float] | dict [str, float] = None) -> str:
    """
    Generate state and property reports for the given EOS.
    
    Parameters
    ----------
    eos : CubicEOS
        The equation of state object.
    eos_type : str
        Type of EOS ('ideal', 'vdw', 'pr', etc.).
    Cp : float | list[float] | dict[str, float], optional
        Heat capacity constants, by default None.

    Returns
    -------
    str
        Formatted state and property reports.
    """

    result = StateReporter({})
    result.add_property('EOS', eos_type)
    result.add_property('T', eos.T, eos.temperature_unit, fstring="{:.2f}")
    result.add_property('P', eos.P, eos.pressure_unit, fstring="{:.2f}")
    if eos_type != 'ideal':
        if hasattr(eos.Z, '__len__'):
            result.add_property('Z (roots)', ', '.join([f"{z:.4f}" for z in eos.Z]), '', fstring=None)
        else:
            result.add_property('Z', eos.Z, '', fstring="{:.4f}")
    if hasattr(eos.v, '__len__'):
        result.add_property('v (roots)', ', '.join([f"{vol:6g}" for vol in eos.v]), f'{eos.volume_unit}/mol', fstring=None)
    else:
        result.add_property('v', eos.v, f'{eos.volume_unit}/mol', fstring="{:6g}")
    if eos_type != 'ideal':
        if hasattr(eos.h, '__len__'):
            result.add_property('h (roots)', ', '.join([f"{h:.2f}" for h in eos.h]), 'J/mol', fstring=None)
        else:
            result.add_property('h', eos.h, 'J/mol', fstring="{:.2f}")
        if hasattr(eos.s, '__len__'):
            result.add_property('s (roots)', ', '.join([f"{s:.2f}" for s in eos.s]), 'J/mol-K', fstring=None)
        else:
            result.add_property('s', eos.s, 'J/mol-K', fstring="{:.2f}")
        if hasattr(eos.h_departure, '__len__'):
            result.add_property('hdep (roots)', ', '.join([f"{hdep:.2f}" for hdep in eos.h_departure]), 'J/mol', fstring=None)
        else:
            result.add_property('hdep', eos.h_departure, 'J/mol', fstring="{:.2f}")
        if hasattr(eos.s_departure, '__len__'):
            result.add_property('sdep (roots)', ', '.join([f"{sdep:.2f}" for sdep in eos.s_departure]), 'J/mol-K', fstring=None)
        else:
            result.add_property('sdep', eos.s_departure, 'J/mol-K', fstring="{:.2f}")
        if eos.T < eos.Tc:
            result.add_property(f'Pvap({eos.T:.2f} K)', eos.Pvap, eos.pressure_unit, fstring="{:.2f}")
            result.add_property(f'Hvap({eos.T:.2f} K)', eos.Hvap, 'J/mol', fstring="{:.2f}")
            result.add_property(f'Svap({eos.T:.2f} K)', eos.Svap, 'J/mol-K', fstring="{:.4f}")
        if eos.P < eos.Pc:
            result.add_property(f'Tsat({eos.P:.2f} {eos.pressure_unit})', eos.Tsat, 'K', fstring="{:.2f}")
        if eos.x is not None:
            result.add_property('Vapor fraction x', eos.x, '', fstring="{:4g}")
    prop = StateReporter({})
    if eos_type != 'ideal':
        prop.add_property('Tc', eos.Tc, 'K', fstring="{:.2f}")
        prop.add_property('Pc', eos.Pc, 'MPa', fstring="{:.2f}")
        if eos_type != 'vdw':
            prop.add_property('omega', eos.omega, '', fstring="{:.3f}")

    if Cp is not None:
        prop.add_property('Tref', eos.Tref, 'K', fstring="{:.2f}")
        prop.add_property('Pref', eos.Pref_MPa, 'MPa', fstring="{:.2f}")
        prop.pack_Cp(Cp, fmts=["{:.2f}", "{:.3e}", "{:.3e}", "{:.3e}"])
    return result.report(), prop.report()

def state(args):
    """
    Calculate and report the state for a single condition using the specified EOS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """

    db = PropertiesDatabase()
    component = db.get_compound(args.n)
    pressure_unit = args.pressure_unit
    volume_unit = args.volume_unit
    Pc_bar = component.Pc
    Tc_K = component.Tc
    Pc = bar_to_unit(Pc_bar, pressure_unit)
    Cp = [component.CpA, component.CpB, component.CpC, component.CpD]
    if component is None:
        print(f"Component '{args.n}' not found in database.")
        return
    match args.eos_type:
        case 'ideal':
            eos = IdealGasEOS(pressure_unit=pressure_unit, volume_unit=volume_unit, Cp=Cp)
        case 'vdw':
            eos = VanDerWaalsEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                Cp = Cp
            )
        case 'pr':
            eos = PengRobinsonEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                omega = component.Omega,
                Cp = Cp
            )
        case 'srk':
            eos = SoaveRedlichKwongEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                omega = component.Omega,
                Cp = Cp
            )
    eos.solve(T=args.T, P=args.P, v=args.v, h=args.h, s=args.s, u=args.u)
    state_report, prop_report = reporters(eos, args.eos_type, Cp=Cp)
    print(state_report)
    if prop_report:
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
    db = PropertiesDatabase()
    component = db.get_compound(args.n) # critical pressures are in bars!
    pressure_unit = args.pressure_unit
    volume_unit = args.volume_unit
    Pc_bar = component.Pc
    Pc = bar_to_unit(Pc_bar, pressure_unit)
    Tc_K = component.Tc
    Cp = [component.CpA, component.CpB, component.CpC, component.CpD]
    if component is None:
        print(f"Component '{args.n}' not found in database.")
        return
    match args.eos_type:
        case 'ideal':
            eos1 = IdealGasEOS(pressure_unit=pressure_unit, volume_unit=volume_unit, Cp=Cp)
            eos2 = IdealGasEOS(pressure_unit=pressure_unit, volume_unit=volume_unit, Cp=Cp)
        case 'vdw':
            eos1 = VanDerWaalsEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                Cp = Cp     
            )
            eos2 = VanDerWaalsEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                Cp = Cp
            )
        case 'pr':
            eos1 = PengRobinsonEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                omega = component.Omega,
                Cp = Cp
            )
            eos2 = PengRobinsonEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                omega = component.Omega,
                Cp = Cp
            )
        case 'srk':
            eos1 = SoaveRedlichKwongEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                omega = component.Omega,
                Cp = Cp
            )
            eos2 = SoaveRedlichKwongEOS(pressure_unit=pressure_unit, volume_unit=volume_unit,
                Tc = Tc_K,
                Pc = Pc,
                omega = component.Omega,
                Cp = Cp
            )

    eos1.solve(T=args.T1, P=args.P1, v=args.v1, h=args.h1, s=args.s1, u=args.u1)
    eos2.solve(T=args.T2, P=args.P2, v=args.v2, h=args.h2, s=args.s2, u=args.u2)
    delta_State = StateReporter({})
    Cp = [component.CpA, component.CpB, component.CpC, component.CpD]
    delta_H = eos2.DeltaH(eos1, Cp)
    delta_S = eos2.DeltaS(eos1, Cp)
    delta_U = eos2.DeltaU(eos1, Cp)
    delta_State.add_property('Delta H', delta_H, 'J/mol', fstring="{:.2f}")
    delta_State.add_property('Delta S', delta_S, 'J/mol-K', fstring="{:.2f}")
    delta_State.add_property('Delta U', delta_U, 'J/mol', fstring="{:.2f}")
    if args.show_states:
        print("State 1:")
        state_1, _ = reporters(eos1, args.eos_type, Cp=Cp)
        print(state_1)
        print("\nState 2:")
        state_2, consts = reporters(eos2, args.eos_type, Cp=Cp)
        print(state_2)
        print("\nProperty differences:")
    print(delta_State.report())
    print("\nConstants used for calculations:")
    print(consts)
    
def cli():
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
        epilog="(c) 2025, Cameron F. Abrams <cfa22@drexel.edu>"
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
        ('Pc', 'critical_pressure', 'critical pressure (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
        ('w', 'acentric_factor', 'acentric factor omega (if component not specified)', float, False),
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, True),
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