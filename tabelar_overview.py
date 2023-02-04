import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

columns = ['nr', 'geb-datum', 'behandlungsbeginn', 'behandlungsende', 'age',
           'sex_f1_m2', 'n_total', 'n_rechts', 'n_links', 'n_bilateral',
           'n_bilateral_lucentis', 'n_bilateral_eylea',
           'n_bilateral_lucentis_eylea', 'n_unilateral', 'n_unilateral_lucentis',
           'n_unilateral_eylea', 'n_lucentis_total', 'n_eylea_total',
           'ozurdex_nein_0__bds_1__rechts_2__links_3',
           'argon_laserbehandelt_nein_0__bds_1__rechts_2__links_3',
           'cat_surgery_nein_0__bds_1__rechts_2__links_3',
           'kapsulotomie_nein_0__bds_1__rechts_2__links_3',
           'oraya_srt_nein_0__bds_1___rechts_2__links_3',
           'eyesurgery_keine_0___ppv_mit_membranpeeling_1___lidkorrektur_(ektropium)_2___tumor_exzision_augenlid_3___dacryocystorhinostomie_4___exzision_fettprolaps_subkonjunktival_5___tpp_6__laseriridotomie,_slt,_te_7__vordere_vitrektomie_bei_kompl._cat-op_8___edta_abrasio_9___parazentese_zum_druckausgleich_10__bh-naht_nach_fk_11',
           'eysurgery_treatmentside_keine_0___bds_1___rechts_2___links_3_',
           'surgery_other_keine_0___diagnostische_eingriffe_zb_herzkatheter_1___kardiologische,_angiologische_eingriffe_2___orthopädische_eingriffe_3___onkologische_eingriffe_4__nierentransplantation_5___dermatologische_eingriffe_6__gastro_eingriffe_7___neurochirurgische_eingriffe_8___zahnärztliche,_kieferchirurgische_eingriffe_9',
           'op',
           'pretreated_(extern_oder_vor_01/11)_nein_0__rechts_2__links_3__bds_1',
           'bcva_baseline_r_bcva_(etdrs)_letters_\nbaseline_r',
           'bcva_baseline_l_bcva_(etdrs)_letters_\nbaseline_l',
           'bcva_endline_r_bcva_(etdrs)_letters_\nnach_letzter_ivi_im_beobachtungszeitraum_r',
           'bcva_endline_l_bcva_(etdrs)_letters_\nnach_letzter_ivi_im_beobachtungszeitraum_l',
           'hyposphagma_bilateral_keine_0___hyposphagma_unilateral_nicht_antikoag._lucentis_1___hyposphagma_unilateral_antikoag._lucentis_2___hyposp._bilat_nicht_antiko._luc_3__hyopsph._bilat._antiko._luc_4___hyposp._unilateral_nicht_antikoa._eylea_5___hyposphagma_unilat_antikoag._eylea_6___hyposph._bilat_nicht_antikoa._eylea_7___hyposphagma_bilat_antikoa._eylea_8___hyposphagma_unilat_nicht_antik._e&l_9',
           'sicca_bilateral_keine_0__sicca-sy_unilat_luc_1___sicca_sy_bilat_luc_2__sicca-sy_unilat_eylea_3__sicca_sy_bilat_eylea_4__sicca-sy_unilat_luc/eylea_5__sicca_sy_bilat_luc/eylea_6',
           'allergy_bilateral_keine_0___allerg._reaktion_unilat_luc_1__allerg_bilat_luc_2__allerg._reaktion_unilat_eylea_3___allergie_bilat_eylea_4',
           'allergy_cause_bilateral_keine_0___konservierungsmittel_1___jod_2___tobrex_3___octenisept_4__desinfektionsmittel_5___unklar_6___floxal_sdu_7___pollen_8',
           'iod_bilateral_keine_0___iod_dekompensation_luc_unilateral_ohne_vorb._glaukom_1__iod_dekomp_luc_unilateral_mit_glaukom_2__iod_dekomp_luc_bilateral_ohne_glaukom_3_iod_dekomp_luc_bilat_mit_glaukom_4_iod_dekomp_eylea_unilateral_ohne_glaukom_5_iod_dekomp_eylea_unilateral_mit_glaukom_6_iod_dekomp_eylea_bilat_ohne_glaukom_7_iod_dekomp_eylea_bilat_mit_glaukom_8',
           'hyposphagma_unilateral_keine_0___hyposphagma_unilateral_nicht_antikoag._lucentis_1___hyposphagma_unilateral_antikoag._lucentis_2___hyposp._bilat_nicht_antiko._luc_3__hyopsph._bilat._antiko._luc_4___hyposp._unilateral_nicht_antikoa._eylea_5___hyposphagma_unilat_antikoag._eylea_6___hyposph._bilat_nicht_antikoa._eylea_7___hyposphagma_bilat_antikoa._eylea_8___hyposphagma_unilat_nicht_antik._e&l_9__hyposphagma_unilat_antikoaguliert_e&l_10',
           'sicca_unilateral_keine_0__sicca-sy_unilat_luc_1___sicca_sy_bilat_luc_2__sicca-sy_unilat_eylea_3__sicca_sy_bilat_eylea_4__sicca-sy_unilat_luc/eylea_5__sicca_sy_bilat_luc/eylea_6',
           'allergy_unialteral_keine_0___allerg._reaktion_unilat_luc_1__allerg_bilat_luc_2__allerg._reaktion_unilat_eylea_3___allergie_bilat_eylea_4',
           'allergy_cause_unilateral_keine_0___konservierungsmittel_1___jod_2___tobrex_3___octenisept_4___desinfektionsmittel_5___unklar_6',
           'iod_unialteral_keine_0___iod_dekompensation_luc_unilateral_ohne_vorb._glaukom_1__iod_dekomp_luc_unilateral_mit_glaukom_2__iod_dekomp_luc_bilateral_ohne_glaukom_3_iod_dekomp_luc_bilat_mit_glaukom_4_iod_dekomp_eylea_unilateral_ohne_glaukom_5_iod_dekomp_unilateral_mit_glaukom_6_iod_dekomp_eylea_bilat_ohne_glaukom_7_iod_dekomp_eylea_bilat_mit_glaukom_8',
           'ocular_ae_bilateral_bilateral_(adverse_event/serious_adverse_event)_freitext_keine_0_lidprobleme_1___mouches,_gk_destruktion,_erf,_hga_2___vk_oder_gk_blutung_3___nh_foramen_und_makulaschichtforamen_4___zvv_oder_vav_unilateral_5___charles_bonnet_6___bh_pathologien_7___uveitis_anterior_8___hh_pathologien_9___doppelbilder_10___amotio_rhetinae_11__endophthalmitis_12___ah_nävus_13___erosio_cornea,_ulcus_cornea_14__iol_dezentriert_15__augenmigräne_16___iritis,_iridocyclitis,_irisdefekt_17___hkl_subluxiert_18_',
           'ocular_ae_bilateral_treatmentside_nicht_betroffen_0__bds_1___rechts_2___links_3',
           'ocular_ae_bilateral_medication_nicht_betroffen_0__bei_lucentis_1___bei_eylea_2___bei_luc_&_eylea_3',
           'ocular_ae_bilateral_time_time_nach_ivi_till_event_nicht_betroffen_0__innert_24h_1___24h-1_woche_2___innerhalb_1monat_3___>1monat_4___unbekannt_5',
           'ocular_ae_unilateral_unilateral_(adverse_event/serious_adverse_event)_freitext_keine_0_lidprobleme_1___mouches,_gk_destruktion,_hga_2___vk_oder_gk_blutung_3___nh_foramen_und_makulaschichtforamen_4___zvv_oder_vav_unilateral_5___charles_bonnet,_visuelle_halluzinationen_6___bh_pathologien_7___uveitis_anterior_8___hh_pathologien_9___doppelbilder_10___amotio_rhetinae_11__endophthalmitis_12___ah_nävus_13___erosio_cornea_14__iol_dezentriert_15__augenmigräne_16___iritis,_iridocyclitis_17__hkl_subluxiert_18_',
           'ocular_ae_unilateral_treatmentside_nicht_betroffen_0__bds_1___rechts_2___links_3___',
           'ocular_ae_unilateral_medication_nicht_betroffen_0__bei_lucentis_1___bei_eylea_2___bei_luc_&_eylea_3',
           'ocular_ae_unilateral_time_time_nach_ivi_till_event_nicht_betroffen_0__innert_24h_1___24h-1_woche_2___innerhalb_1monat_3___>1monat_4___unbekannt_5',
           'systemic_ae_bilateral_(adverse_event/serious_adverse_event)_freitext_keine_0_tod_1___hospitalisation_grund_unbekannt_2___hosp_sturzereignis,_fraktur_3___hosp_elektive_eingriffe_4___lagerungsbedingte_beschwerden,_rückensz_5___cvi_6__onkolog._path._7___hosp_stoffwechselerkrankungen,_infektionen,_dialyse,_fieber_8___hosp_lungenentz,_akute_dyspnoe._9__synkope_10__schwindel,_kopfschmerzen,_transiente_sehstörungen,_parästhesien_11__herzinfarkt_12___tvt_13__hosp_psychiatrisch_14__facialisparese_15__lungenembolie_16___dermatologische_probleme,_haarausfall_17___herzrhythmusstörungen,_herzrasen,_vhfli,_kreislaufbeschwerden_18___appendizitis_19___rheumatologische_pathologie_20___hosp_mediüberdosierung_21___epileptischer_anfall_22',
           'systemic_ae_bilateral_medication_nicht_betroffen_0__bei_lucentis_1___bei_eylea_2___bei_luc_&_eylea_3',
           'systemic_ae_bilateral_time_time_nach_ivi_till_event_nicht_betroffen_0__innert_24h_1___24h-1_woche_2___innerhalb_1monat_3___>1monat_4___unbekannt_5',
           'systemic_ae_unilateral_(adverse_event/serious_adverse_event)_freitext_keine_0_tod_1___hospitalisation_grund_unbekannt_2___hosp_sturzereignis,_fraktur_3___hosp_elektive_eingriffe_4___lagerungsbedingte_beschwerden,_rückensz_5___cvi_6__onkolog._path._7___hosp_stoffwechselerkrankungen,_infektionen,_gürtelrose_8___hosp_lungenentz.dyspnoe_9__synkope_10__schwindel,_kopfschmerzen,_transiente_sehstörungen,_paresen_11__herzinfarkt_12___tvt,_thrombosen_13__hosp_psychiatrisch_14__facialisparese_15__lungenembolie_16__hosp_dermatologisch_17__herzrhythmusstörungen,_vhfli,_herzprobleme_18___gi-pathologien,_appendizitis_19___rheumatologische_pathologie_20___hosp_mediüberdosierung_21___epileptischer_anfall_22___aneurysmaruptur_23___epistaxis_24',
           'systemic_ae_unilateral_medication_nicht_betroffen_0___bei_lucentis_1___bei_eylea_2___',
           'systemic_ae_unilateral_time_time_nach_ivi_till_event_nicht_betroffen_0__innert_24h_1___24h-1_woche_2___innerhalb_1monat_3___>1monat_4___unbekannt_5',
           'bemerkungen',
           'n_months_observationperiod_beobachtungszeitraum_t&e_in_anz_monate',
           'n_months_t&e_prn_anzahl_mt_im_t&e+prn_modus',
           'prn_eye_endline_prn_auge_ende_beobachtungszeitraum_keines_0__bds_1__rechts_2__links_3',
           'prn_eye_most_time_prn_auge_meiste_zeit_während_beobachtungszeitraum_keines_0__rechts_2__links_3',
           'n_months_bilateral_t&e_anzahl_mt_unter_bilateraler_t&e',
           'n_bilateral_t&e_anzahl_bilaterale_t&e_injektionen',
           'n_bilateral_t&e_prn_anzahl_bilaterale_t&e_+prn_injektionen',
           'startmodus_treatment_beginn_beobachtungszeitraum_bilat_t&e_1__t&e+prn_2__unilat_t&e_3',
           'bilateral_t&e_pathway_1._bilaterale_t&e:_im_verlauf_wechsel_auf_nicht_betroffen_1__t&e+prn_2___exit_t&e_3_ozurdex_bds_4_unilat_t&e_5',
           'bilateral_t&e_to_prn_cause_bilaterale_t&e_auf_t&e+prn:_begründung:_nicht_betroffen_1__rez._am_prn_auge_>_12-14_wo_2__weniger_aktivität_am_prn_auge,_bei_kurzen_intervallen_am_t&e_auge_3__flüssigkeit_am_prn_auge_toleriert_da_nicht_visusrelevant_4__flüssigkeit_am_prn_auge_toleriert_wegen_red_visuspotential_5__auf_pat_wunsch_6__st._n._cvi_7',
           'unilateral_pathway_3._unilat_t&e:_im_verlauf_wechsel_auf_nicht_betroffen_1__bilat_t&e_2__t&e+prn_3__',
           'prn_pathway_2._t&e+prn:_im_verlauf_wechsel_auf:_nicht_betroffen_1__bilat_t&e_2__unilat_t&e_(prn_stopp)_3__unverändert_t&e+prn_4__exit_t&e_5',
           'prn_t&e_cause_t&e+prn_beibehalten:_nicht_betroffen_1__flüssigkeit_toleriert_an_prn_auge_wegen_schlechtem_visus_2__rezidive_am_prn_auge_>12-14_wo_3__flüssigkeit_am_prn_auge_toleriert_da_nicht_visusrelevant_4__ergibt_sich_aus_unterschiedlicher_krankheitsaktivität_5__awdp_6__flüssigkeit_toleriert_an_prn_auge_da_stabil_7',
           'prn_t&e_to_bilateral_t&e_cause_t&e+prn:_grund_für_wechsel_auf_t&e_bds:_nicht_betroffen_1__rezidive_am_prn_auge_<12-14_wo_2___visus_prn_auge_gut_3___beides_4',
           'coordination_beginning_bilateral_t&e_synchronisierung_zu_beginn_bilat_t&e_nicht_erfolgt_1___beginn_synchronisiert_unterschiedl_intervall_2__synchronisiert_in_gleichem_intervall_3',
           'coordination_same_intervall_mehrheit_der_anz_bilat_injektionen:_synchronisierung_in_gleichen_intervallen_ja_1__nein_2',
           'coordination_different_intervall_mehrheit_der_anz_bilat_injektionen:_synchronisierung_in_unterschiedlichen_intervallen_ja_1__nein_2',
           'coordination_base_synchronisierungsbasis_ergibt_sich_aus_unterschiedl._erkrankungsaktivität_1__bei_besserem_auge_abstand_verkürzt_zu_gunsten_der_synchronisierung_2__flüssigkeit_toleriert_am_schlechteren_auge_zu_gunsten_synchronisierung_3__gleiche_krankheitsaktivität_bds_4__keine_synchronisierung_5',
           'coordination_cause_falls_2_oder_3_begründung_für_synchronisierung:_nicht_betroffen_0___situation_am_besseren_auge_nicht_stabil/vorgängig_rezidive_unter_längeren_intervallen_1___schlechter_az_/awdp_/_malcompliance_/bestmögliche_koordination_r/l_2___flüssigkeit_toleriert_am_schlechteren_auge_bei_schlechter_visusprognose_3__bei_schlechter_visusprognose_am_besseren_auge_anpassung_an_intervall_des_schlechteren,_besser_sehenden_auges_4__exsudat_stabil_/_sehr_wenig_flüssigkeit_regredient_uner_längeren_intervallen_5__i.r._erster_regulär_4_wöchentlicher_t&e_intervalle_6',
           'no_coordination_keine_synchronisierung_gefunden_nicht_betroffen_0___ja_1__nein_2',
           'no_coordination_cause_wieso_keine_synchronisierung_gefunden?_nicht_betroffen_0___intervallfindungsphase_/_unterschiedliche_krankheitsaktivität_1__sehr_unterschiedliche_behandlungsstadien_(behandlungsanfang_vs._ende)_2___awdp_/will_keine_bilateralen_injektionen_/_keine_koordination_nötig_bei_langen_intervallen_und_niedriger_terminfrequenz_3__st._n._cvi_4',
           'malcompliance_incompliant_ja_1__nein_2',
           'will_keine_bilateralen_injektionen_ja_1__nein_2',
           'treatment_ending_behandlungsabschluss_t&e_rezidivfrei_nach_3x_14_wo_1__behandlung_läuft_weiter_2__behandlung_abgebrochen_3__t&e_auge_rezidivfrei_oder_wechsel_zu_prn,_prn_auge_lohnt_sich_nicht_mehr_wegen_schlechtem_visus_4__weiterführende_therapie_mit_ozurdex_5',
           'treatment_abortion_cause_gründe_für_behandlungsabbruch_nicht_betroffen_1__arzt-/klinikwechsel_2__krankheit/hospitalisation_oder_tod_des_pat_3__lost_to_follow-up_4__pat_will_keine_ivi_mehr_5__kein_therapieansprechen_6',
           'bemerkung']

data = pd.read_excel("data_small.xlsx")
for col in data.columns:
       original_col = col
       # while '  ' in col:
       #        col = col.replace('  ', ' ')
       # col = col.replace(' ', '_')
       # while '__' in col:
       #        col = col.replace('__', '_')
       col = col.split(' ')[0]
       data.rename(columns={original_col: col.lower()}, inplace=True)
print(data.columns)
# data['new_col'] = data['Bemerkung']
# print(data)
data.hist("n_bilateral", bins=12)
plt.show()
# print(data.columns)
# create new column
print("Hello world")
n_tot = data.shape[0]
print(f'number of people {n_tot}')
print(f'number of IVI {sum(data["n_total"])}')
# data1 = data.query("")
# n_allergy1 = data.query(
# 'Allergy_cause_bilateral keine 0 __ Konservierungsmittel 1 __ Jod 2 __ Tobrex 3 __ Octenisept 4__ Desinfektionsmittel 5 __ unklar 6 __ Floxal SDU 7 __ Pollen 8')
# print(f"n allergy1 n people {n_allergy1.shape[0]}")
# for col in data.columns:
#     print(data[col])
medic1 = data.query("n_bilateral_lucentis")
medic2 = data.query("n_bilateral_eylea")
n_medic1 = medic1.shape[0]
n_medic2 = medic2.shape[0]
# n_medic2 = n_medic1
# adversarial = []
labels = ['adver_effect1']
means1 = [n_medic1]
means2 = [n_medic2]
std1 = [n_medic1 ** 0.5]
std2 = [n_medic2 ** 0.5]
width = 0.35  # the width of the bars: can also be len(x) sequence

plt.figure()
fig, ax = plt.subplots()

ax.bar(labels, means1, width, yerr=std1, label='medic1')
ax.bar(labels, means2, width, yerr=std2, bottom=means1, label='medic2')
plt.legend()
plt.show()

labels = ['G1']
men_means = [n_medic1]
women_means = [n_medic2 * 2]
women_means = men_means
# men_std = [2]
# women_std = [3]

data['bcva_diff'] = np.abs(data['bcva_baseline_r'] - data['bcva_baseline_l'])
diff_lt_10 = data.query('bcva_diff > 10')
print(diff_lt_10.shape[0])
print(diff_lt_10['n_bilateral'].sum())
plt.figure()
plt.scatter(data['bcva_diff'], data['n_bilateral'])
plt.title('aoeuaeo')
plt.show()
