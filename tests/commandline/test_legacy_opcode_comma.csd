<CsTest>
description = "Legacy opcode separators preserve existing statement and function calls"

[expect]
exit = 0
output = ["Checked legacy opcode commas"]
</CsTest>
<CsoundSynthesizer>
<CsOptions>
-n -d -m0
</CsOptions>
<CsInstruments>
sr = 8192
ksmps = 32
nchnls = 1
0dbfs = 1
giCaptured init 0
giEchoed init 0
giRan init 0
gkOscillatorAmp init 0

opcode iCapture(iValue:i, iExpected:i):void
  if iValue != iExpected then
    prints "Expected %g from iCapture, got %g\n", iExpected, iValue
    exitnow -1
  endif
  giCaptured += 1
endop

opcode iEcho(iArg:i):i
  giEchoed += 1
  xout iArg
endop

opcode NamedOpcode(SText:S):void
  if strcmp(SText, "text") != 0 then
    exitnow -1
  endif
  giRan += 1
endop

opcode Pair(iValue:i):(i, i)
  xout iValue, iValue + 1
endop

instr 1
  icount = 0
  iVisits = 0
loop3:
  iVisits += 1
  loop_lt,icount,1,4,loop3
  if icount != 4 || iVisits != 4 then
    prints "loop_lt ran %g times and counted to %g\n", iVisits, icount
    exitnow -1
  endif

  iEcho = 37
  ; iEcho is also an opcode name. iCapture receives the variable.
  iCapture iEcho, 37
  iCapture(iEcho, 37)
  iExpected = 37
  iCapture,iEcho,iExpected
  iCapture,abs(iEcho),iExpected
  iValues[] fillarray 37
  iCapture,iValues[0],iExpected
  iCapture,(iEcho + 1),38
  if giCaptured != 6 || giEchoed != 0 then
    prints "Expected six captures and no echo calls, got %g and %g\n", giCaptured, giEchoed
    exitnow -1
  endif

  iEcho abs, -15
  if iEcho != 15 || giEchoed != 0 then
    prints "Expected abs to write 15 without an echo call, got %g and %g\n", iEcho, giEchoed
    exitnow -1
  endif

  NamedOpcode:Opcode = create(opcoderef("NamedOpcode"))
  run NamedOpcode,"text"
  if giRan != 1 then
    prints "Expected run to call the opcode object once, got %g calls\n", giRan
    exitnow -1
  endif

  iLeft,iRight Pair 5
  iFunctionLeft, iFunctionRight = Pair(8)
  if iLeft != 5 || iRight != 6 || iFunctionLeft != 8 || iFunctionRight != 9 then
    prints "Multi-output calls returned %g, %g, %g, %g\n", iLeft, iRight, iFunctionLeft, iFunctionRight
    exitnow -1
  endif
  iLegacyLeft,iLegacyRight Pair,11
  if iLegacyLeft != 11 || iLegacyRight != 12 then
    prints "Legacy multi-output call returned %g, %g\n", iLegacyLeft, iLegacyRight
    exitnow -1
  endif

  a1 oscili,0.5,440
  aFunction = oscili(0.5, 440)
  kLegacy downsamp a1
  kFunction downsamp aFunction
  if kLegacy != kFunction then
    printks "Legacy oscillator returned %g, expected %g\n", 0, kLegacy, kFunction
    exitnowk -1
  endif
  gkOscillatorAmp += abs(kLegacy)
endin

instr 2
  if i(gkOscillatorAmp) <= 0 then
    prints "The oscillator never produced a signal\n"
    exitnow -1
  endif
  prints "Checked legacy opcode commas\n"
endin
</CsInstruments>
<CsScore>
i 1 0 .015625
i 2 .03125 0
e
</CsScore>
</CsoundSynthesizer>
