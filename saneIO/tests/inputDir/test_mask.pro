map = MRDFITS('polaris_plw.fits',0,header)

mask = map GT 1.0

dataPath = '/home/abeelen/sanepic/Workspace/saneIO/tests/data/PLW_only/'
READCOL, 'listFile.list',filenames, FORMAT='A'


CCR_tot = 0

iFile = 0
FOR iFile = 0, N_ELEMENTS(filenames)-1 DO BEGIN
   ra  = MRDFITS(dataPath+filenames[iFile], 'lon',/SILENT)
   dec = MRDFITS(dataPath+filenames[iFile], 'lat',/SILENT)
   
   adxy, header, ra, dec, x, y
   
   CCR = WHERE(mask[x,y] EQ 1, nCCR)
   print, iFile, nCCR
   CCR_tot += nCCR

ENDFOR

iFile = 0
ra  = MRDFITS(dataPath+filenames[iFile], 'lon',/SILENT)
dec = MRDFITS(dataPath+filenames[iFile], 'lat',/SILENT)
signal = MRDFITS(dataPath+filenames[iFile], 'signal',/SILENT)
adxy, header, ra, dec, x, y

CCR = WHERE(mask[x,y] EQ 1, nCCR)

FOR iChan = 0, 40 DO signal[*,iChan] -= MEDIAN(signal[*,iChan])

TRIANGULATE, x[ccr], y[ccr], tr, b 

imdisp, TRIGRID(x[ccr], y[ccr], signal[ccr], tr) 

END
