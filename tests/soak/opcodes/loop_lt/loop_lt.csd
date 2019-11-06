<CsoundSynthesizer>
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 2
0dbfs  = 1

seed 9864
gisine ftgen 0, 0, 2^10, 10, 1

instr 1 ;master instrument

ininstr = 5 ;number of called instances
indx = 0
loop:
      prints "play instance %d\\n", indx
ipan  random 0, 1
ifreq random 100, 1000
iamp  = 1/ininstr
event_i "i", 10, 0, p3, iamp, ifreq, ipan
loop_lt indx, 1, ininstr, loop

endin

instr 10

ipeak random 0, 1 ;where is the envelope peak
asig  poscil3 p4, p5, gisine
aenv  transeg 0, p3*ipeak, 6, 1, p3-p3*ipeak, -6, 0
aL,aR pan2 asig*aenv, p6
      outs aL, aR

endin

</CsInstruments>
<CsScore>
i1 0 10
e
</CsScore>
</CsoundSynthesizer>