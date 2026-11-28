<CsTest>
description = "Recursive UDO delay taps follow the caller's separate stereo readers across reinit"
[expect]
exit = 0
</CsTest>
<CsoundSynthesizer>
<CsOptions>
-n -d -m0
</CsOptions>
<CsInstruments>
sr = 8192
ksmps = 32
nchnls = 2
0dbfs = 1
gkSamplesChecked init 0

opcode RecursiveTap, a, iiii
  iTap, iCount, iBase, iInterpolate xin
  iDelay = (iBase + iTap * 30) / sr
  if iInterpolate == 0 then
    aCurrent deltap iDelay
  else
    aCurrent deltapi a(iDelay)
  endif
  aRest = 0
  if iTap + 1 < iCount then
    aRest RecursiveTap iTap + 1, iCount, iBase, iInterpolate
  endif
  xout aCurrent + aRest
endop

instr StereoTaps
  ; Each channel gets one impulse at sample 0 and another at sample 512.
  aLeft mpulse .5, .0625
  aRight mpulse -.25, .0625
  kTapCount init 4
  if gkSamplesChecked == 384 then
    kTapCount = 5
    reinit DELAYS
  endif
DELAYS:
  aReadL delayr .0625
  aTapL RecursiveTap 0, i(kTapCount), 37, 0
  delayw aLeft
  aReadR delayr .0625
  aTapR RecursiveTap 0, i(kTapCount), 45.5, 1
  delayw aRight
  rireturn

  kFrame = 0
  while kFrame < ksmps do
    kSinceImpulse = gkSamplesChecked % 512
    kExpectedL = 0
    kExpectedR = 0
    kTap = 0
    while kTap < kTapCount do
      if kSinceImpulse == 37 + kTap * 30 then
        kExpectedL = .5
      endif
      ; A half-sample delay divides the impulse equally across two samples.
      if kSinceImpulse == 45 + kTap * 30 || kSinceImpulse == 46 + kTap * 30 then
        kExpectedR = -.125
      endif
      kTap += 1
    od
    kLeft vaget kFrame, aTapL
    kRight vaget kFrame, aTapR
    if kLeft != kExpectedL || kRight != kExpectedR then
      printks "UDO delay taps at sample %g returned %g, %g instead of %g, %g\n", 0, gkSamplesChecked, kLeft, kRight, kExpectedL, kExpectedR
      exitnowk -1
    endif
    gkSamplesChecked += 1
    kFrame += 1
  od
endin

instr CheckCompletion
  if i(gkSamplesChecked) != 896 then
    prints "Expected 896 checked samples, got %g\n", i(gkSamplesChecked)
    exitnow -1
  endif
endin
</CsInstruments>
<CsScore>
i "StereoTaps" 0 .109375
i "CheckCompletion" .125 0
e
</CsScore>
</CsoundSynthesizer>
