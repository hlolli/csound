<CsoundSynthesizer>
<CsInstruments>

; Initialize the global variables.
sr = 44100
kr = 4410
ksmps = 10
nchnls = 1

; Instrument #1.
instr 1
  ; Create random numbers at k-rate in the range -1 to 1 
  ; with a uniform distribution, seed=10.
  k1 rnd31 1, 0, 10
        
  printks "k1=%f\\n", 0.1, k1
endin


</CsInstruments>
<CsScore>

; Play Instrument #1 for one second.
i 1 0 1
e


</CsScore>
</CsoundSynthesizer>