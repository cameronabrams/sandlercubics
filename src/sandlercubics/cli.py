from .eos import *

from sandlerprops.properties import PropertiesDatabase

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
def state(args):
    db = PropertiesDatabase()
    component = db.get_compound(args.component)
    if component is None:
        print(f"Component '{args.component}' not found in database.")
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
    if args.temperature is not None and args.pressure is not None:
        eos.T = args.temperature
        eos.P = args.pressure
        print(f"At T={args.temperature} K and P={args.pressure} MPa, the molar volume is {eos.v:.6f} m^3/mol")
    else:
        print("Please provide both temperature and pressure to calculate the state.")

def cli():
    subcommands = {
        'state': dict(
            func = state,
            help = 'work with a cubic equation of state for a single state'
        ),
        'delta': dict(
            func = None,
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
    command_parsers['state'].add_argument(
        '-T',
        '--temperature',
        type=float,
        help='temperature in K'
    )
    command_parsers['state'].add_argument(
        '-P',
        '--pressure',
        type=float,
        help='pressure in MPa'
    )
    command_parsers['state'].add_argument(
        '-Tc',
        '--critical-temperature',
        type=float,
        help='critical temperature in K (if component not specified)'
    )
    command_parsers['state'].add_argument(
        '-Pc',
        '--critical-pressure',
        type=float,
        help='critical pressure in MPa (if component not specified)'
    )
    command_parsers['state'].add_argument(
        '-w',
        '--acentric-factor',
        type=float,
        help='acentric factor omega (if component not specified)'
    )
    command_parsers['state'].add_argument(
        '-n',
        '--component',
        type=str,
        help='component name (e.g., methane, ethane, etc.)'
    )
    command_parsers['state'].add_argument(
        '-eos',
        '--eos-type',
        type=str,
        choices=['ideal', 'vdw', 'pr'],
        default='vdw',
        help='type of cubic equation of state to use'
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