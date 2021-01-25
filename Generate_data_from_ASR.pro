PRO Generate_data_from_ASR

  InPath = 'D:\Results\Onetime\20200706_1stevidence_cb_sub_updraft_GRL\1st_draft\'
;  OutPath = 'D:\Results\Onetime\20200706_1stevidence_cb_sub_updraft_GRL\1st_draft\Cases\'
  OutPath = 'D:\Results\Onetime\20200706_1stevidence_cb_sub_updraft_GRL\1st_draft\Cases_GRL_2ndsubmit\' 
  DLPath = 'D:\Data\SGP_DL\'
  ceilpath = 'D:\Data\SGP_VCEIL\'
  METPath = 'D:\Data\SGP_MET\'
  ecorPath = 'D:\Data\SGP_qcecor\'
  sondePath = 'D:\Data\SGP_sonde\'
  kazrPath = 'D:\Data\SGP_ARSCLKAZR\'

  READCOL, InPath + 'CaseDate_SGP.txt', $
    dateArray,overpassArray,$
    Comment = '*', Format = 'A,A'

  ;  READCOL, InPath + 'CaseDate_SGP_new_20200717.txt', $
  ;    dateArray,overpassArray,$
  ;    Comment = '*', Format = 'A,A'

  ;  dateArray = ['20120314']
  ;  overpassArray = ['1630']
  ;
  ncases = N_ELEMENTS(dateArray)

  ;  YZKAZRACMES, OutPath,$
  ;    dateArray, overpassArray,$
  ;    KAZRPath = kazrpath,$
  ;    CEILPath = ceilpath,$
  ;    METPath = metpath,$
  ;    sondepath = sondepath,$
  ;    Plotfig = 1,$
  ;    plotsonde = 1, $
  ;    timeinterval = 2,$
  ;    readkazr = 1,$
  ;    readmet = 1,$
  ;    readceil = 1,$
  ;    readsonde = 0, $
  ;    outputkazr = 1,$
  ;    outputmet = 0,$
  ;    outputceil = 0,$
  ;    outputsonde = 0,$
  ;    ind_st = 169,$
  ;    ind_nd = 207
  ;    ;ind_nd = ncases - 1
  ;
  ;  PRINT,''

  YZDLICMES, OutPath,$
    dateArray, overpassArray,$
    DLPath = dlpath,$
    CEILPath = ceilpath,$
    METPath = metpath,$
    sondepath = sondepath,$
    ecorpath = ecorpath,$
    Plotfig = 0,$
    plotsonde = 0, $
    timeinterval = 2,$
    readdl = 1,$
    readmet = 1,$
    readceil = 1,$
    readsonde = 1, $
    readecor = 1, $
    outputdl = 1,$
    outputmet = 0,$
    outputceil = 0,$
    outputsonde = 0,$
    outputecor = 0,$
    ind_st = 0,$
;    ind_nd = 6
    ind_nd = ncases - 1

  PRINT,''
END