DEFSYSV, "!GD", GETDATA_CONSTANTS()

gd_pointer_ref = GD_OPEN('tempDir_ref/dirfile/OD162_0x500022F9L_SpirePacsParallel_polaris_notdf_all_sanepic/Noise_data',/VERBOSE)
gd_pointer     = GD_OPEN(    'tempDir/dirfile/OD162_0x500022F9L_SpirePacsParallel_polaris_notdf_all_sanepic/Noise_data',/VERBOSE)


list_ref = gd_field_list(gd_pointer_ref)
list     = gd_field_list(gd_pointer)

n_ref = gd_nframes(gd_pointer_ref)
n     = gd_nframes(gd_pointer)




det_ref = STRARR(N_ELEMENTS(list_ref))
FOR I=0, N_ELEMENTS(list_ref)-1 DO $
   det_ref[I] = (STRSPLIT(list_ref[I],"_",/EXTRACT,/REGEX))[0]

list = list[1:N_ELEMENTS(list)-2]
det = STRARR(N_ELEMENTS(list))
FOR I=0, N_ELEMENTS(list)-1 DO $
   det[I] = (STRSPLIT(list[I],"_",/EXTRACT,/REGEX))[8]


FOR index = 0, N_ELEMENTS(list)-1 DO BEGIN
;;   good = WHERE(STRMATCH(STRLOWCASE(list_ref),STRLOWCASE(list[index])) EQ 1,nGood)
;;    good = WHERE(STRMATCH(STRLOWCASE("Indexes_"+list_ref),STRLOWCASE(list[index])) EQ 1,nGood)
   good = WHERE(STRMATCH(STRLOWCASE(det_ref), STRLOWCASE(det[index])) EQ 1, nGood)

   IF nGood EQ 1 THEN BEGIN
      data = gd_getdata(gd_pointer, list[index],num_frames=n)
      data_ref = gd_getdata(gd_pointer_ref, list_ref[good[0]], num_frames=n_ref)
      print, minmax(data-data_ref)
   ENDIF ELSE BEGIN
      stop
   ENDELSE

ENDFOR

END
