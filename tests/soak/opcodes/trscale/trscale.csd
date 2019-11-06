<CsoundSynthesizer>
<CsInstruments>

sr = 44100
ksmps = 32
nchnls = 2
0dbfs  = 1

instr 1

kpitch  = p4
ain	diskin2	"fox.wav", 1
fs1,fsi2 pvsifd	ain, 2048, 512, 1		; ifd analysis
fst	partials fs1, fsi2, .003, 1, 3, 500	; partial tracking
fscl	trscale	fst, kpitch			; frequency scale
aout	tradsyn	fscl, 1, 1, 500, 1		; resynthesis 
        outs    aout, aout

endin
</CsInstruments>
<CsScore>
f1 0 8192 10 1

i 1 0 3 1.5	;up a 5th
i 1 3 3 3	;two octaves higher
e
</CsScore>
</CsoundSynthesizer>