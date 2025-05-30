# Branch Tests

## Syntax Overview

Today we're starting to build larger branch networks and analyze the data that comes out of it. This file explains the notation used in files such as `./2025-05-21/branch_networks.csv`

Syntax is

`Network = { Wire1 { Wire2[Termination], Subnetwork } }`

For example this network

`{WIRE_A{WIRE_B[T],WIRE_C{WIRE_D[S],WIRE_E{WIRE_F[O]}}}}`

maps to the corresponding network

```ASCII
-> WIRE A -|-> WIRE_C -|-> WIRE_E -> WIRE_F -> [open]
	 |	         |
  WIRE_B 	   WIRE_D
	 |	         |
[terminated]  [short]
```

Where [X] is the termination type, which can be any of the following:

* O: Open
* S: Short
* T: Terminated
* R=_: Resistor value
* L=_: Inductor value
* C=_: Capacitance value

Series elements are not supported, and every unattached wire must specify it's ending type

A rib looks like

`{A00{A01{A02{A03{A04{B05[O],B04[O]},B03[O]},B02[O]},B01[O]},B00[O]}}`

or this,

`{A00{B00[O],A01{B01[O],A02{B02[O],A03{B03[O],A04{B04[O],B05[O]}}}}}}`

or this,

```
{
	A00 {
		B00[O],
		A01 {
			B01[O],
			A02 {
				B02[O],
				A03 {
					B03[O],
					A04 {
						B04[O],
						B05[O]
					}
				}
			}
		}
	}
}
```

They're all the same network

This format should be easily parseable by code for further analysis.

# Rules

* The top level network has to be contained within brackets `{}`
* Every wire must have a connection at its end
  * Other wires (1 or 2)
  * Termination type
* Each wire must have it's own .lws recording for proper simulation
