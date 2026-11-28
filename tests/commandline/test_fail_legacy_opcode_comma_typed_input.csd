<CsTest>
description = "A legacy opcode comma does not allow a type annotation on an input"

[expect]
exit = "nonzero"
stderr = ["Expected an input argument after the comma"]
</CsTest>
<CsoundSynthesizer>
<CsOptions>
-n -d -m0
</CsOptions>
<CsInstruments>
instr 1
  iValue = 7
  print,iValue:i
endin
</CsInstruments>
<CsScore>
i 1 0 .01
e
</CsScore>
</CsoundSynthesizer>
