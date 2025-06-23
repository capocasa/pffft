## Nim wrapper for PFFFT - Pretty Fast FFT library
## 
## This is a Nim binding for the PFFFT library, which provides
## fast single-precision FFT transforms optimized for SIMD instructions.
##
## Original C library by Julien Pommier, based on FFTPACKv4.

from os import parentDir, `/`

{.compile: currentSourcePath.parentDir() / "pffft/pffft.c".}

type
  ## Opaque struct holding internal precomputed data
  PffftSetup* = ptr object
  
  ## Direction of the transform
  PffftDirection* {.size: sizeof(cint).} = enum
    PFFFT_FORWARD = 0
    PFFFT_BACKWARD = 1
  
  ## Type of transform
  PffftTransform* {.size: sizeof(cint).} = enum
    PFFFT_REAL = 0
    PFFFT_COMPLEX = 1

## Create a new PFFFT setup for transforms of size N
## Returns nil if N is not suitable for transformation
proc pffft_new_setup*(N: cint, transform: PffftTransform): PffftSetup {.
  cdecl, importc: "pffft_new_setup".}

## Destroy a PFFFT setup and free its memory
proc pffft_destroy_setup*(setup: PffftSetup) {.
  cdecl, importc: "pffft_destroy_setup".}

## Perform a Fourier transform
## The output is stored in the most efficient order for back-transforms or convolution
## Transforms are not scaled: BACKWARD(FORWARD(x)) = N*x
## 'work' should point to N floats (2*N for complex), or nil to use stack
proc pffft_transform*(setup: PffftSetup, input: ptr UncheckedArray[float32], output: ptr UncheckedArray[float32], 
                     work: ptr UncheckedArray[float32], direction: PffftDirection) {.
  cdecl, importc: "pffft_transform".}

## Similar to pffft_transform but ensures output is ordered as interleaved complex numbers
proc pffft_transform_ordered*(setup: PffftSetup, input: ptr UncheckedArray[float32], output: ptr UncheckedArray[float32],
                             work: ptr UncheckedArray[float32], direction: PffftDirection) {.
  cdecl, importc: "pffft_transform_ordered".}

## Reorder frequency components to canonical order (interleaved complex numbers)
## Call after pffft_transform(..., PFFFT_FORWARD) for proper frequency ordering
proc pffft_zreorder*(setup: PffftSetup, input: ptr UncheckedArray[float32], output: ptr UncheckedArray[float32],
                    direction: PffftDirection) {.
  cdecl, importc: "pffft_zreorder".}

## Perform multiplication of frequency components and accumulate: dft_ab += (dft_a * dft_b) * scaling
## Arrays should be obtained from pffft_transform(..., PFFFT_FORWARD) without reordering
proc pffft_zconvolve_accumulate*(setup: PffftSetup, dft_a: ptr UncheckedArray[float32], dft_b: ptr UncheckedArray[float32],
                                dft_ab: ptr UncheckedArray[float32], scaling: float32) {.
  cdecl, importc: "pffft_zconvolve_accumulate".}

## Allocate aligned memory buffer (16-byte boundary for SIMD)
proc pffft_aligned_malloc*(nb_bytes: csize_t): ptr UncheckedArray[float32] {.
  cdecl, importc: "pffft_aligned_malloc".}

## Free aligned memory buffer
proc pffft_aligned_free*(p: ptr UncheckedArray[float]) {.
  cdecl, importc: "pffft_aligned_free".}

## Returns 4 if SSE/Altivec support is enabled, 1 otherwise
proc pffft_simd_size*(): cint {.
  cdecl, importc: "pffft_simd_size".}

