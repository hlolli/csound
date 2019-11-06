<CsoundSynthesizer>
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 2
0dbfs   =1

instr 1

ift  = ftlen(p4)
     print ift
aout loscil3 .8, 4, p4
     outs aout, aout

endin 
</CsInstruments> 
<CsScore> 
f 1 0 0 1 "fox.wav" 0 0 0 ;Csound computes tablesize
f 2 0 0 1 "beats.wav" 0 0 0 ;Csound computes tablesize

i 1 0 3 1 ;"fox.wav"
i 1 3 3 2 ;"beats.wav"

e
</CsScore>
</CsoundSynthesizer>