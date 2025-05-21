# Branch Tests

Today we're starting to build larger branch networks and analyze the data that comes out of it. This file explains the notation used in `branch_networks.csv`

Syntax is

`{WIRE_A{WIRE_B}{WIRE_C{WIRE_D}{WIRE_E}}}`

For the corresponding network

```ASCII
LiveWire SSTDR -> WIRE A -|-> WIRE_C -|-> WIRE_E
			  |	      |
			WIRE_B 	    WIRE_D
```

A rib looks like

`{A00{A01{A02{A03{A04{B05}{B04}}{B03}}{B02}}{B01}}{B00}}`

or this,

`{A00{B00}{A01{B01}{A02{B02}{A03{B03}{A04{B04}{B05}}}}}}`

or this,

```
{A00
	{B00}
	{A01
		{B01}
		{A02
			{B02}
			{A03
				{B03}
				{A04
					{B04}
					{B05}
				}
			}
		}
	}
}
```

They're all the same network

This format should be easily parseable by code for further analysis.

In the event of a 4-way split, one can just list an additional wire split-off, like `{A00{B00}{C00}{D00}}`, but this usage is unlikely and not guarenteed to be supported.
This syntax does not support loops, and likely neither will the code, however it should be extensible to add features on wires, such as open/close signals.

`{A00{B00[S]}}{A01{B01[O]}{B02[T]}}}`

With codes 'S', 'O', and 'T' for short, open, and terminated.
