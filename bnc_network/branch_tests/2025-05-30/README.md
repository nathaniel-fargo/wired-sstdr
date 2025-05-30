# Branching Network Data
In this data set our goal is to measure how deep into a network we can 'see' into by changing the termination type. We'll need to create a sufficiently large network, and then terminate different 'legs' of the network with either open, terminated, or shorted ends. By comparing data across termination types we can accurately quantify our ability to 'see' changes at different depths. 

Heres the full network

```
{
    D00 {
        E00 [O],
        D01 {
            E01 [O],
            D02 {
                E02 [O],
                D03 {
                    E03 [O],
                    D04 {
                        E04 [O],
                        D05 {
                            E06 [O],
                            D06 {
                                E07 [O],
                                D07 {
                                    A00 [O],
                                    D08 {
                                        A01 [O],
                                        D09 {
                                            A02 [O],
                                            B00 {
                                                A03 [O],
                                                B01 {
                                                    C01 [O],
                                                    B03 {
                                                        C02 [O],
                                                        B04 {
                                                            C03 [O],
                                                            B05 {
                                                                F00 [O],
                                                                B06 {
                                                                    G00 [O],
                                                                    B07 {
                                                                        G01 [O]
                                                                        B08 {
                                                                            B09 [O],
                                                                            B10 [O]
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
```