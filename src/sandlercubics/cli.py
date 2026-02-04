# Author: Cameron F. Abrams, <cfa22@drexel.edu>

import os
import sys

import argparse as ap
import numpy as np

from importlib.metadata import version

from sandlerprops.properties import PropertiesDatabase

from .eos import CubicEOS
from .idealgas import IdealGasEOS
from .soaveredlichkwong import SoaveRedlichKwongEOS
from .pengrobinson import PengRobinsonEOS
from .vanderwaals import VanDerWaalsEOS

from sandlermisc import ureg, StateReporter

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

def state_subcommand(args):
    """
    Calculate and report the state for a single condition using the specified EOS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    match args.eos_type:
        case 'ideal':
            eos = IdealGasEOS()
        case 'vdw':
            eos = VanDerWaalsEOS()
        case 'pr':
            eos = PengRobinsonEOS()
        case 'srk':
            eos = SoaveRedlichKwongEOS()
    if args.n is not None:
        eos.set_compound(args.n)
    if args.Tc is not None:
        eos.Tc = args.Tc
    if args.Pc is not None:
        eos.Pc = args.Pc
    if args.w is not None:
        eos.omega = args.w
    if args.Cp is not None:
        eos.Cp = args.Cp
    for p in 'TPhsuv':
        v = getattr(args, p, None)
        if v is not None:
            setattr(eos, p, v)
    additional_vars = ['Z']
    if eos.T < eos.Tc:
        additional_vars.extend(['Pvap', 'Hvap', 'Svap'])
    if eos.P < eos.Pc:
        additional_vars.extend(['Tsat'])
    property_notes = {
        'Pvap': f'at {eos.T.to(ureg.kelvin):g}',
        'Hvap': f'at {eos.T.to(ureg.kelvin):g}',
        'Svap': f'at {eos.T.to(ureg.kelvin):g}',
        'Tsat': f'at {eos.P.to(ureg.megapascal):g}',
    }
    print(eos.report(additional_vars=additional_vars, show_parameters=args.show_props, property_notes=property_notes))

def delta_subcommand(args):
    """
    Calculate and report property differences between two states using the specified EOS.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.
    """
    match args.eos_type:
        case 'ideal':
            eos1 = IdealGasEOS()
            eos2 = IdealGasEOS()
        case 'vdw':
            eos1 = VanDerWaalsEOS()
            eos2 = VanDerWaalsEOS()
        case 'pr':
            eos1 = PengRobinsonEOS()
            eos2 = PengRobinsonEOS()
        case 'srk':
            eos1 = SoaveRedlichKwongEOS()
            eos2 = SoaveRedlichKwongEOS()
    if args.n is not None:
        eos1.set_compound(args.n)
        eos2.set_compound(args.n)
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

    for p in 'TPhsuv':
        v = getattr(args, f'{p}1', None)
        if v is not None:
            setattr(eos1, p, v)
        v = getattr(args, f'{p}2', None)
        if v is not None:
            setattr(eos2, p, v)
    state_1 = eos1.report(additional_vars=['Z'], show_parameters=args.show_props)
    state_2 = eos2.report(additional_vars=['Z'], show_parameters=args.show_props)
    delta = eos1.delta(eos2, additional_vars=['Pv', 'Z'])
    print(f"State-change calculations for {args.n} using {eos1.description}:")
    if args.show_states:
        print()
        two_states = ["State 1:                                 State 2:"]
        for line1, line2 in zip(state_1.splitlines(), state_2.splitlines()):
            two_states.append(f"{line1:<36s}     {line2}")
        print("\n".join(two_states))
        print()
        print("Property changes:")
    for p in ['T', 'P', 'h', 's', 'u', 'v', 'Pv', 'Z']:
        if p in delta:
            val = delta[p]
            eq = ' =' if p == 'Pv' else '  ='
            print(f'Δ{p}{eq} {val: 6g}')
    
def cli():
    """
    Command-line interface for sandlercubics package.
    """

    subcommands = {
        'state': dict(
            func = state_subcommand,
            help = 'work with a cubic equation of state for a single state'
        ),
        'delta': dict(
            func = delta_subcommand,
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
        ('P', 'pressure', 'pressure in MPa', float, False),
        ('T', 'temperature', 'temperature in K (always in K)', float, False),
        ('v', 'molar_volume', 'molar volume in m3/mol', float, False),
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
        ('P1', 'pressure1', 'pressure of state 1 in MPa', float, False),
        ('T1', 'temperature1', 'temperature of state 1 in K (always in K)', float, False),
        ('v1', 'molar_volume1', 'molar volume of state 1 in m3/mol', float, False),
        ('h1', 'enthalpy1', 'molar enthalpy of state 1 in J/mol', float, False),
        ('s1', 'entropy1', 'molar entropy of state 1 in J/mol-K', float, False),
        ('u1', 'internal_energy1', 'molar internal energy of state 1 in J/mol', float, False),
        ('P2', 'pressure2', 'pressure of state 2 in MPa', float, False),
        ('T2', 'temperature2', 'temperature of state 2 in K (always in K)', float, False),
        ('v2', 'molar_volume2', 'molar volume of state 2 in m3/mol', float, False),
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

    if args.func == state_subcommand:
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