# Sandlercubics

> Digitized cubic equations of state from Sandler's 5th ed.

Sandlercubics implements a python interface to the cubic equations of state found in  _Chemical, Biochemical, and Engineering Thermodynamics_ (5th edition) by Stan Sandler (Wiley, USA). It should be used for educational purposes only.

Currently only the generalized vdW and generalized Peng-Robison equations for pure substances are implemented.

## Installation 

sandlercubics is available via `pip`:

```sh
pip install sandlercubics
```

## Usage

### Command-line interface

A volumetric calculation example on methane using Peng-Robinson:

```sh
$ sandlercubics state -T 400 -P 0.5 -eos pr -n methane
EOS  = pr
T    = 400.00 K
P    = 0.50 MPa
Z    = 1.00
v    = 0.006628 m3/mol
Hdep = -54.75 J/mol
Sdep = -0.11 J/mol-K
Tc    = 190.40 K
Pc    = 4.60 MPa
omega = 0.011
```

State calculation example:
```sh
$ sandlercubics delta -T1 350 -P1 7.5 -T2 400 -P2 15.5 -n methane -eos pr --show-states
State 1:
EOS  = pr
T    = 350.00 K
P    = 7.50 MPa
Z    = 0.93
v    = 0.000359 m3/mol
Hdep = -989.93 J/mol
Sdep = -2.13 J/mol-K

State 2:
EOS  = pr
T    = 400.00 K
P    = 15.50 MPa
Z    = 0.95
v    = 0.000204 m3/mol
Hdep = -1412.32 J/mol
Sdep = -2.86 J/mol-K

Property differences:
Delta H = 2416.63 J/mol
Delta S = 0.02 J/mol-K
Delta U = 2883.75 J/mol

Constants used for calculations:
Tc    = 190.40 K
Pc    = 4.60 MPa
omega = 0.011
CpA   = 19.25 J/mol-K
CpB   = 5.213e-02 J/mol-K^2
CpC   = 1.197e-05 J/mol-K^3
CpD   = -1.132e-08 J/mol-K^4
```
### API

Below we create a `PengRobinsonEOS` object to reproduce the above calculation:

```python
>>> from sandlercubics.eos import PengRobsinsonEOS
>>> from sandlerprops.properties import PropertiesDatabase
>>> db = ProperitesDatabase()
>>> m = db.get_compound('methane')
>>> s1 = PengRobinsonEOS(Tc=m.Tc, Pc=m.Pc/10, omega=m.Omega)
>>> s1.T = 400
>>> s1.P = 0.5
>>> s1.v.item()  # it is a np float
0.0066279171348771915
```

## Release History

* 0.2.0
    * uses `StateReporter`
    * `delta` subcommand implemented
* 0.1.1
    * fixed erroneous thank-you message
* 0.1.0
    * Initial version, implements vdw and PengRobinson

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/sandlercubics/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
