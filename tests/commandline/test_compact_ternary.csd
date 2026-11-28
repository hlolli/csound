<CsTest>
description = "Compact ternary expressions preserve typed calls and declarations"

[expect]
exit = 0
output = ["Checked compact ternary expressions"]
</CsTest>
<CsoundSynthesizer>
<CsOptions>
-n -d -m0
</CsOptions>
<CsInstruments>
sr = 8192
ksmps = 32
nchnls = 1
gkChecks init 0

opcode Check(kActual:k, kExpected:k, SForm:S):void
  if kActual != kExpected then
    printks "%s returned %g, expected %g\n", 0, SForm, kActual, kExpected
    exitnowk -1
  endif
endop

instr 1
  kOdd = timeinstk() % 2
  kAmp = 1 + 2 * kOdd
  kMax = 2
  kOffset = 4
  kFactor = 5
  expected:k = 2 + kOdd

  kPlain = kAmp>kMax?kAmp:kMax
  Check(kPlain, expected, "compact max")

  kCompound = kAmp>kMax?kAmp+kOffset:kMax*kFactor
  Check(kCompound, 10 - 3*kOdd, "arithmetic in both arms")

  kNested = kAmp>kMax?(kMax<kOffset?kAmp:kMax):(kMax>kOffset?kOffset:kFactor)
  Check(kNested, 5 - 2*kOdd, "ternary in each arm")

  kChain = kAmp>kMax?kAmp:kMax<kOffset?kOffset:kFactor
  Check(kChain, 4 - kOdd, "right-associated ternary")

  kContinued = kAmp>kMax ? \
      kAmp:kMax
  Check(kContinued, expected, "continued expression")

  kUnescaped = kAmp>kMax?kAmp+
      kOffset:kMax*kFactor
  Check(kUnescaped, 10 - 3*kOdd, "newline after a binary operator")

  kFirst init 1
  kSecond init 0
  offsets:k[] fillarray 11, 22
  kIndexed = offsets[kAmp>kMax?kFirst:kSecond]
  Check(kIndexed, 11 + 11*kOdd, "ternary in an array index")

  iStart = 0
  iEnd = 1
  kSlice = kAmp>kMax?sumarray:k(offsets[iStart : iEnd])+kOffset:kFactor
  Check(kSlice, 5 + 32*kOdd, "slice in the true arm")

  picked:B = kAmp>kMax?truek:falsek
  kBoolean = picked ? kAmp : kMax
  Check(kBoolean, expected, "compact boolean arms")

  typed:k = kAmp>kMax?abs:k(-kAmp):abs:k(-kMax)
  Check(typed, expected, "typed calls in both arms")

  kMixed = kAmp>kMax?kAmp:abs:k(-kMax)
  Check(kMixed, expected, "compact separator before a typed call")

  gkChecks += 1
endin

instr 2
  if i(gkChecks) != 4 then
    prints "Expected four checks of each form, got %g\n", i(gkChecks)
    exitnow -1
  endif
  prints "Checked compact ternary expressions\n"
endin
</CsInstruments>
<CsScore>
i 1 0 .015625
i 2 .03125 0
e
</CsScore>
</CsoundSynthesizer>
