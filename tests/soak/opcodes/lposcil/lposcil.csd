<CsoundSynthesizer>
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 2
0dbfs  = 1

instr 1

kcps  = 1.5	; a fifth up
kloop = 0	; loop start time in samples
kend  = 45000	; loop end time in samples

asig lposcil 1, kcps, kloop, kend, 1
     outs asig, asig

endin
</CsInstruments>
<CsScore>
; Its table size is deferred,
; and format taken from the soundfile header.
f 1 0 0 1 "beats.wav" 0 0 0

; Play Instrument #1 for 6 seconds.
; This will loop the drum pattern several times.
i 1 0 6

e
</CsScore>
</CsoundSynthesizer>