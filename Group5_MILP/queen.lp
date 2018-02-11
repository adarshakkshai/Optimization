\<queen>
\<exampleName>queen </exampleName>
\<category>milp</category>
\<solver>scip</solver>
\<inputMethod>CPLEX</inputMethod>
\<lp><![CDATA[
\\ This model allows to determine where and how many queens can be
\\ placed on the 8x8 chess board not beating each other.

Maximize
      x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 +
      x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 +
      x31 + x32 + x33 + x34 + x35 + x36 + x37 + x38 +
      x41 + x42 + x43 + x44 + x45 + x46 + x47 + x48 +
      x51 + x52 + x53 + x54 + x55 + x56 + x57 + x58 +
      x61 + x62 + x63 + x64 + x65 + x66 + x67 + x68 +
      x71 + x72 + x73 + x74 + x75 + x76 + x77 + x78 +
      x81 + x82 + x83 + x84 + x85 + x86 + x87 + x88

Subject To
h1:   x11 + x12 + x13 + x14 + x15 + x16 + x17 + x18 <= 1
h2:   x21 + x22 + x23 + x24 + x25 + x26 + x27 + x28 <= 1
h3:   x31 + x32 + x33 + x34 + x35 + x36 + x37 + x38 <= 1
h4:   x41 + x42 + x43 + x44 + x45 + x46 + x47 + x48 <= 1
h5:   x51 + x52 + x53 + x54 + x55 + x56 + x57 + x58 <= 1
h6:   x61 + x62 + x63 + x64 + x65 + x66 + x67 + x68 <= 1
h7:   x71 + x72 + x73 + x74 + x75 + x76 + x77 + x78 <= 1
h8:   x81 + x82 + x83 + x84 + x85 + x86 + x87 + x88 <= 1

v1:   x11 + x21 + x31 + x41 + x51 + x61 + x71 + x81 <= 1
v2:   x12 + x22 + x32 + x42 + x52 + x62 + x72 + x82 <= 1
v3:   x13 + x23 + x33 + x43 + x53 + x63 + x73 + x83 <= 1
v4:   x14 + x24 + x34 + x44 + x54 + x64 + x74 + x84 <= 1
v5:   x15 + x25 + x35 + x45 + x55 + x65 + x75 + x85 <= 1
v6:   x16 + x26 + x36 + x46 + x56 + x66 + x76 + x86 <= 1
v7:   x17 + x27 + x37 + x47 + x57 + x67 + x77 + x87 <= 1
v8:   x18 + x28 + x38 + x48 + x58 + x68 + x78 + x88 <= 1

e1:   x17 + x28                                     <= 1
e2:   x16 + x27 + x38                               <= 1
e3:   x15 + x26 + x37 + x48                         <= 1
e4:   x14 + x25 + x36 + x47 + x58                   <= 1
e5:   x13 + x24 + x35 + x46 + x57 + x68             <= 1
e6:   x12 + x23 + x34 + x45 + x56 + x67 + x78       <= 1
e7:   x11 + x22 + x33 + x44 + x55 + x66 + x77 + x88 <= 1
e8:   x21 + x32 + x43 + x54 + x65 + x76 + x87       <= 1
e9:   x31 + x42 + x53 + x64 + x75 + x86             <= 1
e10:  x41 + x52 + x63 + x74 + x85                   <= 1
e11:  x51 + x62 + x73 + x84                         <= 1
e12:  x61 + x72 + x83                               <= 1
e13:  x71 + x82                                     <= 1

d1:   x21 + x12                                     <= 1
d2:   x31 + x22 + x13                               <= 1
d3:   x41 + x32 + x23 + x14                         <= 1
d4:   x51 + x42 + x33 + x24 + x15                   <= 1
d5:   x61 + x52 + x43 + x34 + x25 + x16             <= 1
d6:   x71 + x62 + x53 + x44 + x35 + x26 + x17       <= 1
d7:   x81 + x72 + x63 + x54 + x45 + x36 + x27 + x18 <= 1
d8:   x82 + x73 + x64 + x55 + x46 + x37 + x28       <= 1
d9:   x83 + x74 + x65 + x56 + x47 + x38             <= 1
d10:  x84 + x75 + x66 + x57 + x48                   <= 1
d11:  x85 + x76 + x67 + x58                         <= 1
d12:  x86 + x77 + x68                               <= 1
d13:  x87 + x78                                     <= 1

Binaries
      x11 x12 x13 x14 x15 x16 x17 x18
      x21 x22 x23 x24 x25 x26 x27 x28
      x31 x32 x33 x34 x35 x36 x37 x38
      x41 x42 x43 x44 x45 x46 x47 x48
      x51 x52 x53 x54 x55 x56 x57 x58
      x61 x62 x63 x64 x65 x66 x67 x68
      x71 x72 x73 x74 x75 x76 x77 x78
      x81 x82 x83 x84 x85 x86 x87 x88

End

\* eof *\

\<comment>
\This is the example queen.lp in CPLEX-LP format from the GLPK distribution

\</comment>

\</queen>
