## Test suite for PFFFT Nim implementation
## Based on test_pffft.c from the original PFFFT library
##
## Copyright (c) 2013 Julien Pommier.
##
## Small test & bench for PFFFT, comparing its performance with the scalar FFTPACK, FFTW, Intel MKL, and Apple vDSP

import math, random, times, strformat, sequtils
import pffft, pffft/fftpack

const
  MAX_OF_TEMPLATE = "max($1, $2)"

template MAX_OF(x, y: untyped): untyped =
  if x > y: x else: y

proc frand(): float32 =
  rand(1.0).float32

proc uclock_sec(): float64 =
  epochTime()

proc norm_inf_rel(v, w: ptr UncheckedArray[float32], N: cint): float32 =
  var max_w = 0.0f
  var max_diff = 0.0f
  
  for k in 0..<N:
    max_w = MAX_OF(max_w, abs(w[k]))
    max_diff = MAX_OF(max_diff, abs(w[k] - v[k]))
  
  assert(max_w > 0)
  return max_diff / max_w

proc pffft_validate_N(N: cint, cplx: bool) =
  let Nfloat = N * (if cplx: 2 else: 1)
  let Nbytes = (Nfloat * sizeof float32).cSizeT
  
  var reference = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  var ref2 = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  var input = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  var output = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  var tmp = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  var tmp2 = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  
  assert(Nbytes > 0)
  let s = pffftNewSetup(N, if cplx: PffftComplex else: PffftReal)
  
  if s == nil:
    echo &"Skipping N={N}, not supported"
    pffftAlignedFree(reference)
    pffftAlignedFree(ref2)
    pffftAlignedFree(input)
    pffftAlignedFree(output)
    pffftAlignedFree(tmp)
    pffftAlignedFree(tmp2)
    return
  
  let scratch = if N < 2000: nil else: cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  
  for pass in 0..1:
    var ref_max = 0.0f
    
    # Compute reference solution with FFTPACK
    if pass == 0:
      # Allocate workspace for FFTPACK
      var wrk = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc((2 * Nbytes + 25 * sizeof float32).cSizeT))
      
      for k in 0..<Nfloat:
        let val = frand() * 2 - 1
        reference[k] = val
        input[k] = val
        output[k] = 1e30
      
      # Compute forward transform with fftpack
      if not cplx:
        # Real FFT
        rffti(N, wrk)
        rfftf(N, reference, wrk)
        copyMem(ref2, reference, Nbytes)
        rfftb(N, ref2, wrk)
        
        # Use our ordering for real ffts instead of the one of fftpack
        let refN = reference[N - 1]
        for k in countdown(N - 2, 1):
          reference[k + 1] = reference[k]
        reference[1] = refN
      else:
        # Complex FFT
        cffti(N, wrk)
        cfftf(N, reference, wrk)
        copyMem(ref2, reference, Nbytes)
        cfftb(N, ref2, wrk)
      
      # Normalize ref2 for verification
      for k in 0..<Nfloat:
        ref2[k] /= N.float32
      
      # Verify FFTPACK back and forth error
      let fftpack_back_and_forth_error = norm_inf_rel(ref2, input, Nfloat)
      assert(fftpack_back_and_forth_error < 1e-3)
      
      pffftAlignedFree(wrk)
    
    for k in 0..<Nfloat:
      ref_max = MAX_OF(ref_max, abs(reference[k]))
    
    # Test forward transform, with different input / output
    if pass == 0:
      pffftTransform(s, input, tmp, scratch, PffftForward)
      copyMem(tmp2, tmp, Nbytes)
      copyMem(tmp, input, Nbytes)
      pffftTransform(s, tmp, tmp, scratch, PffftForward)
     
      for k in 0..<Nfloat:
        assert(tmp2[k] == tmp[k])
      
      # Test reordering
      pffftZreorder(s, tmp, output, PffftForward)
      pffftZreorder(s, output, tmp, PffftBackward)
      for k in 0..<Nfloat:
        assert(tmp2[k] == tmp[k])
      pffftZreorder(s, tmp, output, PffftForward)
    else:
      # Pass 1: canonical ordering of transform coeffs
      pffftTransformOrdered(s, input, tmp, scratch, PffftForward)
      copyMem(tmp2, tmp, Nbytes)
      copyMem(tmp, input, Nbytes)
      pffftTransformOrdered(s, tmp, tmp, scratch, PffftForward)
      for k in 0..<Nfloat:
        assert(tmp2[k] == tmp[k])
      copyMem(output, tmp, Nbytes)
    
    block:
      # Error for the forward transform when compared with fftpack
      let max_forward_transform_error = norm_inf_rel(output, reference, Nfloat)
      if not (max_forward_transform_error < 1e-3):
        echo &"{(if cplx: \"CPLX\" else: \"REAL\")} forward PFFFT mismatch found for N={N} relative error={max_forward_transform_error}"
        assert(false)
      
      if pass == 0:
        pffftTransform(s, tmp, output, scratch, PffftBackward)
      else:
        pffftTransformOrdered(s, tmp, output, scratch, PffftBackward)
      
      copyMem(tmp2, output, Nbytes)
      copyMem(output, tmp, Nbytes)
      
      if pass == 0:
        pffftTransform(s, output, output, scratch, PffftBackward)
      else:
        pffftTransformOrdered(s, output, output, scratch, PffftBackward)
      
      for k in 0..<Nfloat:
        assert(tmp2[k] == output[k])
        output[k] *= 1.0f / N.float32
      
      # Error when transformed back to the original vector
      let max_final_error_rel = norm_inf_rel(output, input, Nfloat)
      if max_final_error_rel > 1e-3:
        echo &"pass={pass}, {(if cplx: \"CPLX\" else: \"REAL\")} IFFFT does not match for N={N}, relative error={max_final_error_rel}"
        assert(false)
    
    # Quick test of the circular convolution in fft domain
    block:
      var conv_err = 0.0f
      var conv_max = 0.0f
      
      pffftZreorder(s, reference, tmp, PffftForward)
      for k in 0..< (Nbytes.int div sizeof float32):
        output[k] = 0
      pffftZconvolveAccumulate(s, reference, reference, output, 1.0)
      pffftZreorder(s, output, tmp2, PffftForward)
      
      for k in countup(0, Nfloat, 2):
        let ar = tmp[k]
        let ai = tmp[k + 1]
        if cplx or k > 0:
          tmp[k] = ar * ar - ai * ai
          tmp[k + 1] = 2 * ar * ai
        else:
          tmp[0] = ar * ar
          tmp[1] = ai * ai
      
      for k in 0..<Nfloat:
        let d = abs(tmp[k] - tmp2[k])
        let e = abs(tmp[k])
        if d > conv_err:
          conv_err = d
        if e > conv_max:
          conv_max = e
      
      if conv_err > 1e-5 * conv_max:
        echo &"zconvolve error ? {conv_err} {conv_max}"
        assert(false)
  
  echo &"{(if cplx: \"CPLX\" else: \"REAL\")} PFFFT is OK for N={N}"
  
  pffftDestroySetup(s)
  pffftAlignedFree(reference)
  pffftAlignedFree(ref2)
  pffftAlignedFree(input)
  pffftAlignedFree(output)
  pffftAlignedFree(tmp)
  pffftAlignedFree(tmp2)
  if scratch != nil:
    pffftAlignedFree(scratch)

proc pffft_validate(cplx: bool) =
  let Ntest = [16'i32, 32'i32, 64'i32, 96'i32, 128'i32, 160'i32, 192'i32, 256'i32, 288'i32, 384'i32, 5'i32*96'i32, 512'i32, 576'i32, 5'i32*128'i32, 800'i32, 864'i32, 1024'i32, 2048'i32, 2592'i32, 4000'i32, 4096'i32, 12000'i32, 36864'i32, 0'i32]
  
  for N in Ntest:
    if N == 0:
      break
    if N == 16 and not cplx:
      continue
    pffft_validate_N(N, cplx)

proc show_output(name: string, N: cint, cplx: bool, flops: float64, t0: float64, t1: float64, max_iter: int) =
  let mflops = flops / 1e6 / (t1 - t0 + 1e-16)
  if flops != -1:
    echo &"N={N:5}, {(if cplx: \"CPLX\" else: \"REAL\")} {name:>16} : {mflops:6.0f} MFlops [t={((t1-t0)/2/max_iter.float64 * 1e9):6.0f} ns, {max_iter} runs]"

proc benchmark_ffts(N: cint, cplx: bool) =
  let Nfloat = if cplx: N * 2 else: N
  let Nbytes = (Nfloat * cint(sizeof float32)).cSizeT
  let X = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  let Y = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  let Z = cast[ptr UncheckedArray[float32]](pffftAlignedMalloc(Nbytes))
  
  var max_iter = 5120000 div N * 4
  when defined(arm) or defined(arm64) or defined(aarch64):
    max_iter = max_iter div 8
  if max_iter == 0:
    max_iter = 1
  
  for k in 0..<Nfloat:
    X[k] = 0  # sqrtf(k+1)
  
  # PFFFT benchmark
  let s = pffftNewSetup(N, if cplx: PffftComplex else: PffftReal)
  if s != nil:
    let t0 = uclock_sec()
    for iter in 0..<max_iter:
      pffftTransform(s, X, Z, Y, PffftForward)
      pffftTransform(s, X, Z, Y, PffftBackward)
    let t1 = uclock_sec()
    pffftDestroySetup(s)
    
    let flops = (max_iter.float64 * 2) * ((if cplx: 5.0 else: 2.5) * N.float64 * ln(N.float64) / ln(2.0))
    show_output("PFFFT", N, cplx, flops, t0, t1, max_iter)
  
  pffftAlignedFree(X)
  pffftAlignedFree(Y)
  pffftAlignedFree(Z)

when isMainModule:
  randomize()
  
  echo "PFFFT validation..."
  pffft_validate(true)   # Complex
  pffft_validate(false)  # Real
  
  echo "\nPFFFFT benchmarks:"
  let Nvalues = [64'i32, 96'i32, 128'i32, 160'i32, 192'i32, 256'i32, 384'i32, 5'i32*96'i32, 512'i32, 5'i32*128'i32, 3'i32*256'i32, 800'i32, 1024'i32, 2048'i32, 2400'i32, 4096'i32]
  let Nmax = 1024'i32 * 1024'i32
  
  for N in Nvalues:
    if N < Nmax:
      benchmark_ffts(N, false)  # Real fft
  
  for N in Nvalues:
    if N < Nmax:
      benchmark_ffts(N, true)   # Complex fft
  
  echo "All tests passed!"
