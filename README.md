# Sandlersteam

> Digitized cubic equations of state from Sandler's 5th ed.

Sandlercubics implements a python interface to the cubic equations of state found in  _Chemical, Biochemical, and Engineering Thermodynamics_ (5th edition) by Stan Sandler (Wiley, USA). It should be used for educational purposes only.

Currently only the generalized vdW and generalized Peng-Robison equations for pure substances are implemented.

## Installation 

Sandlersteam is available via `pip`:

```sh
pip install sandlercubics
```

## Usage

### Command-line interface

A volumetric calculation on methane using Peng-Robninson:

```bash
$ sandlercubics state -T 400 -P 0.5 -eos pr -n methane
At T=400.0 K and P=0.5 MPa, the molar volume is 0.0066279 m^3/mol
```

### API

Below we create a `PengRobinsonEOS` object to reproduce the above calculation

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

* 0.1.0
    * Initial version, implements vdw and PengRobinson

## Meta

Cameron F. Abrams – cfa22@drexel.edu

Distributed under the MIT license. See ``LICENSE`` for more information.

[https://github.com/cameronabrams](https://github.com/cameronabrams/)

## Contributing

1. Fork it (<https://github.com/cameronabrams/sandlersteam/fork>)
2. Create your feature branch (`git checkout -b feature/fooBar`)
3. Commit your changes (`git commit -am 'Add some fooBar'`)
4. Push to the branch (`git push origin feature/fooBar`)
5. Create a new Pull Request
