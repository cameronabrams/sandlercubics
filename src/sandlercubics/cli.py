from unittest import result
from .eos import *

from sandlerprops.properties import PropertiesDatabase
from sandlermisc.statereporter import StateReporter

import argparse as ap
banner = """
 __                 _ _           
/ _\ __ _ _ __   __| | | ___ _ __ 
\ \ / _` | '_ \ / _` | |/ _ \ '__|
_\ \ (_| | | | | (_| | |  __/ |   
\__/\__,_|_| |_|\__,_|_|\___|_|   
               _     _            
     ___ _   _| |__ (_) ___ ___   
    / __| | | | '_ \| |/ __/ __|  
   | (__| |_| | |_) | | (__\__ \  
    \___|\__,_|_.__/|_|\___|___/  

(c) 2025, Cameron F. Abrams <cfa22@drexel.edu>
"""
def reporters(eos: CubicEOS, eos_type: str, Cp: float | list[float] | dict [str, float] = None) -> str:
    result = StateReporter({})
    result.add_property('EOS', eos_type)
    result.add_property('T', eos.T, 'K', fstring="{:.2f}")
    result.add_property('P', eos.P, 'MPa', fstring="{:.2f}")
    if eos_type != 'ideal':
        result.add_property('Z', eos.Z, '', fstring="{:.2f}")
    result.add_property('v', eos.v, f'{eos.volume_unit}/mol', fstring="{:.6f}")
    if eos_type != 'ideal':
        result.add_property('Hdep', eos.h_departure, 'J/mol', fstring="{:.2f}")
        result.add_property('Sdep', eos.s_departure, 'J/mol-K', fstring="{:.2f}")
    prop = StateReporter({})
    if eos_type != 'ideal':
        prop.add_property('Tc', eos.Tc, 'K', fstring="{:.2f}")
        prop.add_property('Pc', eos.Pc, 'MPa', fstring="{:.2f}")
        if eos_type != 'vdw':
            prop.add_property('omega', eos.omega, '', fstring="{:.3f}")
    if Cp is not None:
        prop.pack_Cp(Cp, fmts=["{:.2f}", "{:.3e}", "{:.3e}", "{:.3e}"])
    return result.report(), prop.report()

def state(args):
    db = PropertiesDatabase()
    component = db.get_compound(args.n)
    if component is None:
        print(f"Component '{args.n}' not found in database.")
        return
    match args.eos_type:
        case 'ideal':
            eos = IdealGasEOS(pressure_unit='mpa', volume_unit='m3')
        case 'vdw':
            eos = GeneralizedVDWEOS(pressure_unit='mpa', volume_unit='m3',
                Tc = component.Tc,
                Pc = component.Pc/10
            )
        case 'pr':
            eos = PengRobinsonEOS(pressure_unit='mpa', volume_unit='m3',
                Tc = component.Tc,
                Pc = component.Pc/10,
                omega = component.Omega
            )
    if args.T is not None and args.P is not None:
        eos.T = args.T
        eos.P = args.P
        state_report, prop_report = reporters(eos, args.eos_type)
        print(state_report)
        if prop_report:
            print(prop_report)
    else:
        print("Please provide both temperature and pressure to calculate the state.")

def delta(args):
    db = PropertiesDatabase()
    component = db.get_compound(args.n) # pressures are in bars!
    if component is None:
        print(f"Component '{args.n}' not found in database.")
        return
    match args.eos_type:
        case 'ideal':
            eos1 = IdealGasEOS(pressure_unit='mpa', volume_unit='m3')
            eos2 = IdealGasEOS(pressure_unit='mpa', volume_unit='m3')
        case 'vdw':
            eos1 = GeneralizedVDWEOS(pressure_unit='mpa', volume_unit='m3',
                Tc = component.Tc,
                Pc = component.Pc/10
            )
            eos2 = GeneralizedVDWEOS(pressure_unit='mpa', volume_unit='m3',
                Tc = component.Tc,
                Pc = component.Pc/10
            )
        case 'pr':
            eos1 = PengRobinsonEOS(pressure_unit='mpa', volume_unit='m3',
                Tc = component.Tc,
                Pc = component.Pc/10,
                omega = component.Omega
            )
            eos2 = PengRobinsonEOS(pressure_unit='mpa', volume_unit='m3',
                Tc = component.Tc,
                Pc = component.Pc/10,
                omega = component.Omega
            )

    if args.T1 is not None and args.P1 is not None:
        eos1.T = args.T1
        eos1.P = args.P1
    if args.T2 is not None and args.P2 is not None:
        eos2.T = args.T2
        eos2.P = args.P2
    if args.T1 is None or args.P1 is None or args.T2 is None or args.P2 is None:
        print("Please provide temperature and pressure for both states to calculate property differences.")
        return
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
        state_1, _ = reporters(eos1, args.eos_type)
        print(state_1)
        print("\nState 2:")
        state_2, consts = reporters(eos2, args.eos_type, Cp)
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
        description='Interact with cubic equations of state in Sandler\'s textbook'
    )
    parser.add_argument(
        '-b',
        '--banner',
        default=False,
        action=ap.BooleanOptionalAction,
        help='toggle banner message'
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
            formatter_class=ap.RawDescriptionHelpFormatter
        )
        command_parsers[k].set_defaults(func=specs['func'])
    state_args = [
        ('P', 'pressure', 'pressure in MPa', float, True),
        ('T', 'temperature', 'temperature in K', float, True),
        ('Pc', 'critical_pressure', 'critical pressure in MPa (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
        ('w', 'acentric_factor', 'acentric factor omega (if component not specified)', float, False),
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, False)
    ]
    for prop, long_arg, explanation, arg_type, required in state_args:
        command_parsers['state'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    
    delta_args = [
        ('P1', 'pressure1', 'pressure of state 1 in MPa', float, True),
        ('T1', 'temperature1', 'temperature of state 1 in K', float, True),
        ('P2', 'pressure2', 'pressure of state 2 in MPa', float, True),
        ('T2', 'temperature2', 'temperature of state 2 in K', float, True),
        ('Pc', 'critical_pressure', 'critical pressure in MPa (if component not specified)', float, False),
        ('Tc', 'critical_temperature', 'critical temperature in K (if component not specified)', float, False),
        ('w', 'acentric_factor', 'acentric factor omega (if component not specified)', float, False),
        ('n', 'component', 'component name (e.g., methane, ethane, etc.)', str, False)
    ]
    for prop, long_arg, explanation, arg_type, required in delta_args:
        command_parsers['delta'].add_argument(
            f'-{prop}',
            f'--{long_arg}',
            dest=prop,
            type=arg_type,
            required=required,
            help=explanation
        )
    command_parsers['state'].add_argument(
        '-eos',
        '--eos-type',
        type=str,
        choices=['ideal', 'vdw', 'pr'],
        default='vdw',
        help='type of cubic equation of state to use'
    )
    command_parsers['delta'].add_argument(
        '-eos',
        '--eos-type',
        type=str,
        choices=['ideal', 'vdw', 'pr'],
        default='vdw',
        help='type of cubic equation of state to use'
    )
    command_parsers['delta'].add_argument(
        '--show-states',
        default=False,
        action=ap.BooleanOptionalAction,
        help='also show the full states for state 1 and state 2'
    )
    args = parser.parse_args()
    if args.banner:
        print(banner)
    if hasattr(args, 'func'):
        args.func(args)
    else:
        my_list = ', '.join(list(subcommands.keys()))
        print(f'No subcommand found. Expected one of {my_list}')
    if args.banner:
        print('Thanks for using sandlercubics!')