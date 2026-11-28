<CsTest>
description = "UDO recursion reports the depth limit and unwinds for the next call"

[expect]
exit = "nonzero"
stderr = ["UDO CountFrames max recursion depth 8 reached", "1 errors in performance"]
output = ["CountFrames returned 8 after the rejected call"]
</CsTest>
<CsoundSynthesizer>
<CsOptions>
-n --recursion-depth=8
</CsOptions>
<CsInstruments>
opcode CountFrames(depth:i):i
  iFrames = 1
  if depth > 1 then
    iFrames += CountFrames(depth - 1)
  endif
  xout iFrames
endop

instr 1
  iFrames = CountFrames(9)
endin

instr 2
  ; A leaked recursion count from the first call would reject this one.
  iFrames = CountFrames(8)
  if iFrames != 8 then
    prints "Expected 8 recursive frames, got %d\n", iFrames
    exitnow 1
  endif
  prints "CountFrames returned 8 after the rejected call\n"
endin
</CsInstruments>
<CsScore>
i 1 0 .01
i 2 .02 .01
e
</CsScore>
</CsoundSynthesizer>
