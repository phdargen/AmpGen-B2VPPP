#!/bin/bash

condor_submit i=scans_v3 o=scan_psipi_P.txt d=scan_psipi_P n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipi_V.txt d=scan_psipi_V n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipi_A.txt d=scan_psipi_A n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipi_T.txt d=scan_psipi_T n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipi_PT.txt d=scan_psipi_PT n=90 submitScan.sub

condor_submit i=scans_v3 o=scan_psiK_P.txt d=scan_psiK_P n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psiK_V.txt d=scan_psiK_V n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psiK_A.txt d=scan_psiK_A n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psiK_T.txt d=scan_psiK_T n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psiK_PT.txt d=scan_psiK_PT n=90 submitScan.sub

condor_submit i=scans_v3 o=scan_psipipi_P.txt d=scan_psipipi_P n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi_S.txt d=scan_psipipi_S n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi_A.txt d=scan_psipipi_A n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi_V.txt d=scan_psipipi_V n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi_T.txt d=scan_psipipi_T n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi_PT.txt d=scan_psipipi_PT n=90 submitScan.sub

condor_submit i=scans_v3 o=scan_psipipi2_S.txt d=scan_psipipi2_S n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi2_A.txt d=scan_psipipi2_A n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi2_V.txt d=scan_psipipi2_V n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi2_T.txt d=scan_psipipi2_T n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_psipipi2_PT.txt d=scan_psipipi2_PT n=90 submitScan.sub

condor_submit i=scans_v3 o=scan_X2Z_SA.txt d=scan_X2Z_SA n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_X2Z_VA.txt d=scan_X2Z_VA n=90 submitScan.sub
condor_submit i=scans_v3 o=scan_X2Z_AA.txt d=scan_X2Z_AA n=90 submitScan.sub

condor_submit i=scans_v3 o=scan_XS2ZA.txt d=scan_XS2ZA n=90 submitScan.sub

#condor_submit i=scans_v3 o=scan_psiKpi_P.txt d=scan_psiKpi_P n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi_S.txt d=scan_psiKpi_S n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi_A.txt d=scan_psiKpi_A n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi_V.txt d=scan_psiKpi_V n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi_T.txt d=scan_psiKpi_T n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi_PT.txt d=scan_psiKpi_PT n=90 submitScan.sub

#condor_submit i=scans_v3 o=scan_psiKpi2_S.txt d=scan_psiKpi2_S n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi2_A.txt d=scan_psiKpi2_A n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi2_V.txt d=scan_psiKpi2_V n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi2_T.txt d=scan_psiKpi2_T n=90 submitScan.sub
#condor_submit i=scans_v3 o=scan_psiKpi2_PT.txt d=scan_psiKpi2_PT n=90 submitScan.sub

##
#condor_submit i=scans_v2 o=scan_X2Z_SA.txt d=scan_X2Z_SA n=90 submitScan.sub
#condor_submit i=scans_v2 o=scan_X2Z_VP.txt d=scan_X2Z_VP n=90 submitScan.sub
#condor_submit i=scans_v2 o=scan_X2Z_VV.txt d=scan_X2Z_VV n=90 submitScan.sub
#condor_submit i=scans_v2 o=scan_X2Z_VA.txt d=scan_X2Z_VA n=90 submitScan.sub
#condor_submit i=scans_v2 o=scan_X2Z_AV.txt d=scan_X2Z_AV n=90 submitScan.sub
#condor_submit i=scans_v2 o=scan_X2Z_AA.txt d=scan_X2Z_AA n=90 submitScan.sub

#condor_submit o=scan_Xs2Z_SP.txt d=scan_Xs2Z_SP2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Z_SA.txt d=scan_Xs2Z_SA2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Z_VP.txt d=scan_Xs2Z_VP2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Z_VV.txt d=scan_Xs2Z_VV2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Z_VA.txt d=scan_Xs2Z_VA2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Z_AV.txt d=scan_Xs2Z_AV2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Z_AA.txt d=scan_Xs2Z_AA2 n=90 submitScan.sub

#condor_submit o=scan_Xs2Zs_SP.txt d=scan_Xs2Zs_SP2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Zs_SA.txt d=scan_Xs2Zs_SA2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Zs_VP.txt d=scan_Xs2Zs_VP2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Zs_VV.txt d=scan_Xs2Zs_VV2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Zs_VA.txt d=scan_Xs2Zs_VA2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Zs_AV.txt d=scan_Xs2Zs_AV2 n=90 submitScan.sub
#condor_submit o=scan_Xs2Zs_AA.txt d=scan_Xs2Zs_AA2 n=90 submitScan.sub

