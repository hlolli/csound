<CsoundSynthesizer>
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 4	;quad
0dbfs  = 1

vbaplsinit 2, 4, 0, 90, 180, 270

instr 1

asig diskin2 "beats.wav", 1, 0, 1		;loop beats.wav

k0,k1,k2,k3 vbapgmove  p3, 1, 2, 310, 180	;change movement of soundsource in
     outq k0*asig,k1*asig,k2*asig,k3*asig	;the rear speakers

endin 
</CsInstruments>
<CsScore>

i 1 0 5

e
</CsScore>
</CsoundSynthesizer>
