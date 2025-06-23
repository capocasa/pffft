
from os import parentDir, `/`

{.compile: currentSourcePath.parentDir() & "/fftpack.c".}

when defined(FFTPACK_DOUBLE_PRECISION):
  type
    FftpackReal* = cdouble
    FftpackInt* = cint
else:
  type
    FftpackReal* = cfloat
    FftpackInt* = cint

{.push importc, cdecl.}

proc cffti*(n: FftpackInt; wsave: ptr UncheckedArray[FftpackReal])

proc cfftf*(n: FftpackInt; c: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

proc cfftb*(n: FftpackInt; c: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

proc rffti*(n: FftpackInt; wsave: ptr UncheckedArray[FftpackReal])
proc rfftf*(n: FftpackInt; r: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])
proc rfftb*(n: FftpackInt; r: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

proc cosqi*(n: FftpackInt; wsave: ptr UncheckedArray[FftpackReal])
proc cosqf*(n: FftpackInt; x: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])
proc cosqb*(n: FftpackInt; x: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

proc costi*(n: FftpackInt; wsave: ptr UncheckedArray[FftpackReal])
proc cost*(n: FftpackInt; x: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

proc sinqi*(n: FftpackInt; wsave: ptr UncheckedArray[FftpackReal])
proc sinqb*(n: FftpackInt; x: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])
proc sinqf*(n: FftpackInt; x: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

proc sinti*(n: FftpackInt; wsave: ptr UncheckedArray[FftpackReal])
proc sint*(n: FftpackInt; x: ptr UncheckedArray[FftpackReal]; wsave: ptr UncheckedArray[FftpackReal])

{.pop.}
