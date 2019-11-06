<CsoundSynthesizer>
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 2
0dbfs  = 1

instr 1 

ims  = 100			;maximum delay time in msec
aout poscil .8, 220, 1		;make a signal
a2   poscil ims/2, 1/p3, 1	;make an LFO
a2   = a2 + ims/2		;offset the LFO so that it is positive
asig vdelay3 aout, a2, ims	;use the LFO to control delay time
     outs  asig, asig

endin

</CsInstruments>
<CsScore>
f1 0 8192 10 1

i 1 0 5 

e

</CsScore>
</CsoundSynthesizer> 