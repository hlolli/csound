<CsTest>
description = "Caller delay taps require the same rate and one unambiguous reader"
[expect]
exit = "nonzero"
stderr = ["deltap: associated delayr not found", "4 errors in performance"]
stderr_regex = [
  'opcode EarlierPassTap\).*deltap: associated delayr not found',
  'opcode ReinitTap\).*deltap: associated delayr not found',
]
output = ["The child tap kept playing after reinit failed."]
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
gkReinitBlocksChecked init 0

opcode SmallerBlockTap, a, 0
  setksmps 16
  aTap deltap .015625
  xout aTap
endop

opcode IndexedTap, a, 0
  aTap deltapi a(.015625), -1
  xout aTap
endop

opcode DefaultTap, a, 0
  aTap deltapi a(.015625)
  xout aTap
endop

opcode LocalDelay, a, i
  iOpen xin
  if iOpen != 0 then
    aRead delayr .0625
  endif
  aTap DefaultTap
  if iOpen != 0 then
    delayw a(0)
  endif
  xout aTap
endop

opcode EarlierPassTap, a, 0
  aTap deltap .015625
  xout aTap
endop

opcode ReinitTap, a, 0
  kCycle init 0
  if kCycle == 8 then
    reinit TAP
  endif
TAP:
  aTap deltap .015625
  rireturn
  kCycle += 1
  xout aTap
endop

instr DifferentBlockSizes
  aRead delayr .0625
  aTap SmallerBlockTap
  delayw a(0)
endin

instr ExplicitIndex
  aRead delayr .0625
  aTap IndexedTap
  delayw a(0)
endin

instr MultipleReaders
  aFirst delayr .0625
  aSecond delayr .125
  aTap DefaultTap
  delayw a(0)
  delayw a(0)
endin

instr SkippedLocalReader
  ; A conditional reader in the UDO keeps its scope even when skipped.
  aRead delayr .0625
  aTap LocalDelay 0
  delayw a(0)
endin

instr EarlierPassReader
  aRead delayr .0625
  kCycle init 0
  if kCycle == 2 then
    reinit CALL_TAP
  endif
  ; Only reinit enters the call, after the reader's init pass has ended.
  goto DONE
CALL_TAP:
  aTap EarlierPassTap
  rireturn
DONE:
  kCycle += 1
endin

instr ChildOnlyReinit
  aRead delayr .0625
  aTap ReinitTap
  delayw a(.5)
  kCycle init 0
  if kCycle >= 8 then
    ; The failed child reinit must leave its first binding safe to use.
    if k(aTap) != .5 then
      printks "The child tap lost its delay reader at cycle %g\n", 0, kCycle
      exitnowk -1
    endif
    gkReinitBlocksChecked += 1
  endif
  kCycle += 1
endin

instr CheckChildReinit
  if i(gkReinitBlocksChecked) != 8 then
    prints "Expected eight blocks after reinit, got %g\n", i(gkReinitBlocksChecked)
    exitnow -1
  endif
  prints "The child tap kept playing after reinit failed.\n"
endin
</CsInstruments>
<CsScore>
i "DifferentBlockSizes" 0 .0625
i "ExplicitIndex" .125 .0625
i "MultipleReaders" .25 .0625
i "SkippedLocalReader" .375 .0625
i "EarlierPassReader" .5 .0625
i "ChildOnlyReinit" .625 .0625
i "CheckChildReinit" .75 0
e
</CsScore>
</CsoundSynthesizer>
