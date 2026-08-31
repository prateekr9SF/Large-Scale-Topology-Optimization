# Manual TAU Instrumentation

To manually instrument a C source file with TAU, first include the TAU header:

```c
#include <TAU.h>
```

Define the main profile and any manual timers:

```c
TAU_PROFILE("main()", "", TAU_DEFAULT);

TAU_PROFILE_TIMER(timer, "MyRegion", "", TAU_USER);
```

Place the timer around the region to be measured:

```c
TAU_PROFILE_START(timer);

myFunction();

TAU_PROFILE_STOP(timer);
```

Source files containing TAU probes must be compiled with the following additional flags:

```makefile
-DPROFILING_ON -I$(TAU_ROOT)/include
```

For example:

```makefile
gcc-12 ... -DPROFILING_ON -I$(TAU_ROOT)/include -c main.c -o main.o
```

Source files without TAU probes can be compiled normally.

The final executable should be linked using the TAU compiler wrapper:

```makefile
tau_cc.sh main.o other.o -o application.exe ...
```

Run the application normally. TAU will generate `profile.*` files that can be inspected using:

```bash
pprof
```

or:

```bash
paraprof
```

> **Important:** `-DPROFILING_ON` must be used when compiling source files containing manual TAU probes. Without it, the program may compile and run successfully but no TAU profile files will be generated.