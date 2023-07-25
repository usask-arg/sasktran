c**************************************
      BLOCK DATA BDMol_2016
c**************************************
c
c      23 November, 2016
c
      INCLUDE  'Species_2016.cmn'
      INCLUDE  'MOLEC.cmn'
c
      DATA (MOLID(I),I=0,NMOL)/ '   All',
c           1        2        3        4        5        6        7
     1'   H2O','   CO2','    O3','   N2O','    CO','   CH4','    O2',
c           8        9       10       11       12       13       14
     2'    NO','   SO2','   NO2','   NH3','  HNO3','    OH','    HF',
c          15       16       17       18       19       20       21
     3'   HCl','   HBr','    HI','   ClO','   OCS','  H2CO','  HOCl',
c          22       23       24       25       26       27       28
     4'    N2','   HCN',' CH3Cl','  H2O2','  C2H2','  C2H6','   PH3',
c          29       30       31       32       33       34       35
     5'  COF2','   SF6','   H2S',' HCOOH','   HO2','     O','ClONO2',
c          36       37       38       39       40       41       42
     6'   NO+','  HOBr','  C2H4',' CH3OH',' CH3Br',' CH3CN','   CF4',
c          43       44       45       46       47       48       49
     7'  C4H2','  HC3N','    H2','    CS','   SO3','  C2N2',' COCl2',
c          50       51       52       53       54       55       56
     8'    SO','  C3H4','   CH3','   CS2'/
c
      END
