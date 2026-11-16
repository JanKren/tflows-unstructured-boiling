# Building

T-Flows consists of four programs, each with its own makefile:

| program | source | binary | purpose |
|---|---|---|---|
| `Generate` | `Sources/Generate` | `Binaries/Generate` | block-structured hexahedral meshes from a `.dom` file |
| `Convert` | `Sources/Convert` | `Binaries/Convert` | meshes of other formats (Fluent, Gambit, Gmsh) to T-Flows `.cfn` and `.dim`; makes polyhedral duals |
| `Divide` | `Sources/Divide` | `Binaries/Divide` | domain decomposition with METIS |
| `Process` | `Sources/Process/Cpu` | `Binaries/Cpu/Process` | the flow solver |

## With the scripts

```bash
Cases/build_tools.sh               # Generate, Divide, Convert
Cases/build_case.sh Cases/Stefan   # Process for one case -> Cases/Stefan/Process
```

## By hand

```bash
cd Sources/Generate && make FORTRAN=gnu
cd ../Divide        && make FORTRAN=gnu
cd ../Convert       && make FORTRAN=gnu
cd ../Process/Cpu   && make FORTRAN=gnu MPI=yes DIR_CASE=/absolute/path/to/case
```

What to know before doing so:

- **The makefiles accept only these variables:** `FORTRAN FCOMP ASSERT DEBUG
  PROF REAL MPI OMP DIR_CASE SORT`. Anything else is an error.
- **Clear `FORTRAN` if the environment sets it.** Compiler modules often
  export `FORTRAN=/path/to/gfortran`, which the makefile takes for its own
  option and rejects. Pass `FORTRAN=gnu`, or `unset FORTRAN` first.
- **PETSc is not needed.** If `PETSC_DIR` is set in the environment and the
  build then fails, build with `PETSC_DIR= make ...`.
- **Do not use `make -j`.** The build races on the Fortran module files and
  stops with "Cannot open module file" on an unrelated module.
- **`Process` lands in `Binaries/Cpu/`**, the other three in `Binaries/`.
- **The MPI library of the build must be the one of the run.**

## Building for a case

`DIR_CASE` points to a folder with a `User_Mod/` of Fortran hooks
(initialisation, output at the end of a time step, sources). The makefile
links these files over the default hooks in `Sources/Process/Cpu/User_Mod`
and compiles them into `Process`. A case without `DIR_CASE` runs with empty
hooks: for the cases here that means no initial interface and no error.

Changing `DIR_CASE` in a tree that was already built for another case goes
wrong in two ways, and the build still reports success:

1. hooks that the previous case had and the new one lacks stay linked and
   are compiled in again;
2. a symbolic link carries the date of its target, so a case file written
   months ago looks older than the object of the last build, and `make`
   links the object of the previous case.

`Cases/build_case.sh` deletes the links, restores every default hook from
`Cases/default_hooks`, removes hooks that are not defaults, touches the case's
files and the modules that include hooks (`User_Mod.f90`, `Vof_Mod.f90`,
`Process_Mod.f90`), builds, and copies the binary into the case folder. After
a case build `git status` shows the hook files of `Sources/Process/Cpu/User_Mod`
as changed and the hooks that only the case has as new;
`git checkout -- Sources/Process/Cpu/User_Mod` and
`git clean -f Sources/Process/Cpu/User_Mod` restore the folder.

A check that costs nothing: the first rows of the case's own output. For the
Scriven case the first row of `bench-data.dat` must give a radius of 50 µm and
an interfacial area different from zero.

## Running

`Generate` and `Divide` read their answers from standard input:

```bash
printf 'stefan\nskip\nskip\nskip\n\n' | Binaries/Generate    # name, then prompts
printf 'bubble\n8\n'                  | Binaries/Divide      # name, number of parts
mpirun -np 8 ./Process > run.log
```

`Process` reads the file `control` of the folder it is started in. Start
parallel runs with as many ranks as `Divide` made parts. Large cases need
`ulimit -s unlimited`.
