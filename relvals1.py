from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# load config yaml
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

diagnoses = ['AMD', "DMÖ", 'ZVV_VAV', 'multiple_diagnoses']
# diagnoses.reverse()
datasets = [pd.read_excel("data_small.xlsx", sheet_name=name).reset_index(drop=True) for name in diagnoses]
pd.set_option('display.max_columns', None, 'display.max_rows', None, 'display.expand_frame_repr', False)
diagnoses_name_mapping = config['diagnosesnames']

########################################################################################################################
# Preprocessing of data
########################################################################################################################

# rename columns, shorten names
for data in datasets:
    for col in data.columns:
        original_col = col
        col = col.split(' ')[0]
        col = col.replace('&', 'x')
        data.rename(columns={original_col: col.lower()}, inplace=True)
print(datasets[0].columns)

all_cols = None
for data in datasets:
    if all_cols is None:
        all_cols = set(data.columns)
    else:
        all_cols.intersection_update(data.columns)
all_cols = list(all_cols)
datasets = [data[all_cols] for data in datasets]

# add new columns
for data in datasets:
    for side in ['r', 'l']:
        data.eval(f'bcva_{side}_diff = bcva_baseline_{side} - bcva_endline_{side}', inplace=True)
    data.eval('bcva_sum_diff = bcva_r_diff + bcva_l_diff', inplace=True)

# start analysis
medications = ['eylea', "lucentis"]
n_patients = [data.shape[0] for data in datasets]
name_data_iter = list(zip(diagnoses, datasets))


def str_mean_std(data: pd.Series, percent=False):
    factor = 100 if percent else 1
    out = f'{data.mean() * factor:.3g} +- {data.std() * factor / data.shape[0] ** 0.5:.2g}'
    if percent:
        out += '%'
    return out


n_total = sum(len(d) for d in datasets)
n_tot_inj = sum(d['n_total'].sum() for d in datasets)
linelength = 12
path_all_plots = Path('plots/all')
path_all_plots.mkdir(parents=True, exist_ok=True)
n_col_all = {}
n_col_std_all = {}
print(":" * linelength)
print("Medication")
print(":" * linelength)
for name, data in name_data_iter:
    n_collected = {}
    n_collected_std = {}
    for _ in [1]:
        df = data.copy()
        n = len(df)
        print("-" * linelength)
        print(f"Data {name}: {n} patients ({n / n_total * 100:.2g}% of all)")
        print("injections/person")
        print("-" * linelength)
        for medication in medications:
            n_treat_col = f"n_{medication}_total"
            message = f'{medication}: both {df[n_treat_col].mean():.3g} '
            for lat in ['uni', 'bi']:
                message += '\t'
                key = f'{medication[:3]}_{lat}'
                n_medic_col_lat = f"n_{lat}lateral_{medication}"
                n_medication = df[n_medic_col_lat]
                # n_collected[key] = n_medication.mean()
                # n_collected_std[key] = n_medication.std() / n_medication.shape[0] ** 0.5
                n_collected[key] = n_medication.sum()
                n_collected_std[key] = n_medication.sum() ** 0.5
                message += f" {n_medic_col_lat}: " + str_mean_std(n_medication)
            print(message)

    plt.figure(figsize=(10, 6))
    title = f"{name}: injections per person by med_lat"
    plt.title(title)
    plt.bar(n_collected.keys(), n_collected.values(), align='center')
    plt.xticks(rotation=45)
    filename = Path(title.replace(" ", "_")).with_suffix('.png')
    plt.savefig(path_all_plots / filename)
    output_file = Path(f'plots/injections/{name}') / filename
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file)
    plt.close()

    n_col_all[name] = n_collected
    n_col_std_all[name] = n_collected_std
n_col_all = {k: {k2: f'{v:.4g} +- {n_col_std_all[k][k2]:.2g}' for k2, v in v.items()} for k, v in n_col_all.items()}

print("Injections by med_lat")
injections_med_lat = pd.DataFrame({k: v.values() for k, v in n_col_all.items()},
                                  index=list(n_col_all.values())[0].keys())
print(injections_med_lat)

print("\n\n\n")

########################################################################################################################
# treatments
########################################################################################################################
n_col_all = {}
n_col_std_all = {}

print(":" * linelength)
print("Treatments")
print(":" * linelength)
for name, data in name_data_iter:
    print("=" * linelength)
    print(f"Data {diagnoses_name_mapping[name]} treatments")
    print("=" * linelength)
    n_collected = {}
    n_collected_std = {}
    for _ in [1]:
        df = data.copy()
        n = len(df)
        treatments = ['txe', 'txe_prn']
        message = ""
        for treatment in treatments:
            message += f'{treatment}:'
            n_treat_col = f"n_bilateral_{treatment}"
            # for lat in ['uni', 'bi']:
            key = f'{treatment[-3:]}'
            # n_medic_col_lat = f"n_{lat}lateral_{medication}"
            n_treatment = df[n_treat_col]
            n_collected[key] = n_treatment.mean()
            n_collected_std[key] = n_treatment.std() / n_treatment.shape[0] ** 0.5
            message += ' ' + str_mean_std(n_treatment) + ' inj/pp '
            message += "\t"
        print(message)

    plt.figure(figsize=(10, 6))
    title = f"{diagnoses_name_mapping[name]}: treatment injections per person"
    plt.title(title)
    plt.bar(n_collected.keys(), n_collected.values(), align='center')
    plt.xticks(rotation=45)
    filename = Path(f"{name}: treatment injections per person".replace(" ", "_")).with_suffix('.png')
    plt.savefig(path_all_plots / filename)
    output_file = Path(f'plots/treatments/{name}') / filename
    output_file.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_file)
    plt.close()
    n_col_all[name] = n_collected
    n_col_std_all[name] = n_collected_std
n_col_all = {k: {k2: f'{v:.3g} +- {n_col_std_all[k][k2]:.2g}' for k2, v in v.items()} for k, v in n_col_all.items()}
print(pd.DataFrame({k: v.values() for k, v in n_col_all.items()}, index=list(n_col_all.values())[0].keys()))

print("\n\n\n")

print(":" * linelength)
print("Treatment reasons")
print(":" * linelength)
n_collected_all = {}
changes_nmax = {"bilateral_txe_pathway": 4, "bilateral_txe_to_prn_cause": 7,
                "prn_txe_cause": 6, "prn_txe_to_bilateral_txe_cause": 4}
for changes, nmax in changes_nmax.items():
    n_collected_all[changes] = {}
    for name, data in name_data_iter:
        # print("=" * linelength)
        # print(f"Data {name} reasons")
        # print("=" * linelength)

        for _ in [1]:
            df = data.copy()
            n = len(df)
            # print("-" * linelength)
            # print("-" * linelength)

            message = ""
            # for lat in ['uni', 'bi']:
            key = f'{name[:3]}'
            # n_medic_col_lat = f"n_{lat}lateral_{medication}"
            message += "\t"
            n_collected_all[changes][key] = []
            for i in range(nmax):
                n = (df[changes] == i + 1).sum()
                std = n ** 0.5
                n_rel = n / df.shape[0] * 100
                std_rel = std / df.shape[0] * 100
                n_formatted = f'{n_rel:.0f} +- {std_rel:.0f}'
                n_collected_all[changes][key].append(n_formatted)
            # print(message)

for change, values in n_collected_all.items():
    print(f"{change}: reasons per patients in percent")
    print(f"mul etc are datasets, index is the encoding")
    print(pd.DataFrame(values, columns=list(values.keys())))
    print()

    # plt.figure(figsize=(10, 6))
    # plt.title(title)
    # plt.bar(n_collected.keys(), n_collected.values(), align='center')
    # plt.xticks(rotation=45)
    # filename = Path(title.replace(" ", "_")).with_suffix('.png')
    # plt.savefig(path_all_plots / filename)
    # output_file = Path(f'plots/treatments/{name}') / filename
    # output_file.parent.mkdir(parents=True, exist_ok=True)
    # plt.savefig(output_file)
    # plt.close()
    # n_col_all[name] = n_collected
    # n_col_std_all[name] = n_collected_std
# n_col_all = {k: {k2: f'{v:.3g} +- {n_col_std_all[k][k2]:.2g}' for k2, v in v.items()} for k, v in n_col_all.items()}
# print(pd.DataFrame({k: v.values() for k, v in n_col_all.items()}, index=list(n_col_all.values())[0].keys()))


data_all = pd.concat(datasets, ignore_index=True)
data_all.reset_index(drop=True, inplace=True)

print("\n\n\n")
print("Baseline values")
print()
print("Number of patients total: ", data_all.shape[0])
print(f"n_total={data_all['n_total'].sum():.0f}")
print(f"n_bilateral={(data_all['n_bilateral'] * 2).sum():.0f}")
print(f"n_unilateral={data_all['n_unilateral'].sum():.0f}")
print(f"n women={data_all.query('sex==1').shape[0]}")
print(f"n men={data_all.query('sex==2').shape[0]}")

print("\n\n\n")

replace_keys = ["lat", "med"]
adverse_effects_cut = {"uni":
                           {"eylea": {},
                            "lucentis": {}
                            },
                       "bi":
                           {"eylea": {},
                            "lucentis": {}
                            }
                       }

print("\n\n\n")

lateral = ['uni', 'bi']
index = pd.MultiIndex.from_product([diagnoses, medications, lateral])
adverse_eff = ['hyposphagma', 'sicca', 'allergy']
adverse_hypo = {
    1: "unilat_non-antikoag",
    2: "unilat_antikoag",
    3: "bilat_non-antikoag",
    4: "bilat_antikoag"}
shifts_med = {'eylea': 4, 'lucentis': 0}
# shifts_lat = {'uni': 0, 'bi': 2}
results = pd.DataFrame(index=index, columns=adverse_hypo.values())
results_abs = pd.DataFrame(index=index, columns=adverse_hypo.values())
for name, data in name_data_iter:
    print("=" * linelength)
    print(f"Data {diagnoses_name_mapping[name]} adverse effects")
    print("=" * linelength)
    n_collected = {}
    n_collected_std = {}

    df = data.copy()
    for medication in medications:
        for lat in lateral:
            key = f'_{medication[:3]}_{lat}'
            # shift_lat = shifts_lat[lat]
            shift = shifts_med[medication]
            for i, adverse_name in adverse_hypo.items():
                series_sel = df[f'hyposphagma_{lat}lateral']
                series_sel = series_sel.apply(lambda x: x[0] if isinstance(x, (list, tuple)) else x)
                if not series_sel.dtype == 'int64':
                    pass
                have_it = (series_sel == i + shift)
                n_have_it = have_it.sum()
                ntot = df["n_total"].sum()
                percent_val = n_have_it / ntot
                results.loc[
                    (name, medication, lat), adverse_name] = f'{percent_val:.1%} +- {n_have_it ** 0.5 / ntot:.1%}'
                results_abs.loc[(name, medication, lat), adverse_name] = f'{n_have_it} +- {n_have_it ** 0.5:.1f}'

print("Adverse effects per data in percent of total injections")
print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
      "\n(medication, uni/bilateral injection taken from table name/encoding)")
print(results)

print("Adverse effects per data in absolute numbers")
print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
      "\n(medication, uni/bilateral injection taken from table name/encoding)")
print(results_abs)

print("\n\n\n")

lateral = ['uni', 'bi']
index = pd.MultiIndex.from_product([diagnoses, medications, lateral])
adverse_eff = {'hyposphagma': 'hyposphagma_{lat}lateral > 0',
               'sicca': 'sicca_{lat}lateral > 0', 'allergy': 'allergy_{lat}lateral > 0'}
# adverse_hypo = {
#     1: "unilat_non-antikoag",
#     2: "unilat_antikoag",
#     3: "bilat_non-antikoag",
#     4: "bilat_antikoag"}
shifts_med = {'eylea': 4, 'lucentis': 0}
# shifts_lat = {'uni': 0, 'bi': 2}
results = pd.DataFrame(index=index, columns=adverse_hypo.values())
results_abs = pd.DataFrame(index=index, columns=adverse_hypo.values())
for name, data in name_data_iter:
    print("=" * linelength)
    print(f"Data {diagnoses_name_mapping[name]} adverse effects")
    print("=" * linelength)
    n_collected = {}
    n_collected_std = {}

    df = data.copy()
    for medication in medications:
        for lat in lateral:
            key = f'_{medication[:3]}_{lat}'
            # shift_lat = shifts_lat[lat]
            shift = shifts_med[medication]
            for i, adverse_name in adverse_hypo.items():
                series_sel = df[f'hyposphagma_{lat}lateral']
                series_sel = series_sel.apply(lambda x: x[0] if isinstance(x, (list, tuple)) else x)
                if not series_sel.dtype == 'int64':
                    pass
                have_it = (series_sel == i + shift)
                n_have_it = have_it.sum()
                ntot = df["n_total"].sum()
                percent_val = n_have_it / ntot
                results.loc[
                    (name, medication, lat), adverse_name] = f'{percent_val:.1%} +- {n_have_it ** 0.5 / ntot:.1%}'
                results_abs.loc[(name, medication, lat), adverse_name] = f'{n_have_it} +- {n_have_it ** 0.5:.1f}'

print("Adverse effects per data in percent of total injections")
print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
      "\n(medication, uni/bilateral injection taken from table name/encoding)")
print(results)

print("Adverse effects per data in absolute numbers")
print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
      "\n(medication, uni/bilateral injection taken from table name/encoding)")
print(results_abs)

# plt.figure(figsize=(10, 5))
# df.plot.scatter(x='bcva_sum_diff', y=n_medic_col, cmap='coolwarm',)
#
#
#
#         # plt.figure(figsize=(10, 5))
#         # df.plot.scatter(x='bcva_sum_diff', y='hyposphagma_bilateral', cmap='coolwarm', )
#
#     # print(data.columns)
#     # create new column
#     # print("hello world")
#     # n_tot = data.shape[0]
#     # print(f'number of people {n_tot}')
#     # print(f'number of ivi {sum(data["n_total"])}')
plt.show()


def pie_diagram(data, name: str):
    """Pie diagram of the different treatment (coord, bilat t&e, coord diff and no coord)."""
    n_bilat_te = data['n_bilateral_txe'].sum() * 2  # 2 injections per person
    n_bilat_teprn = data['n_bilateral_txe_prn'].sum() * 2  # 2 injections per person

    coord_same = data.query('coordination_same_intervall == 1')['n_bilateral_txe'].sum() * 2
    coord_diff = data.query('coordination_different_intervall == 1')['n_bilateral_txe'].sum() * 2
    coord_none = data.query('no_coordination == 1')['n_bilateral_txe'].sum() * 2
    # TODO: horizontal bar diagrams below that show the causes of the different coordination types


    all_coord = coord_same + coord_diff + coord_none
    # hack, what's the right number?
    n_bilat_te = all_coord
    assert all_coord == n_bilat_te, f'{all_coord} != {n_bilat_te}'

    # plt.figure(figsize=(10, 10))

    kwargs = dict(autopct='%1.1f%%', textprops={'fontsize': 10},
                  counterclock=True, startangle=90)
    double_plot = False
    if double_plot:

        plt.pie([n_bilat_te,
                 n_bilat_teprn],
                labels=['bilat t&e', 'bilat t&e prn'],
                **kwargs)
        pie = plt.pie([coord_same,
                       coord_diff,
                       coord_none,
                       n_bilat_teprn],
                      labels=['coord same', 'coord diff', 'no coord', ''],
                      **kwargs)
        # change color of text/bkg in pie for readability
        for wedge in pie[0]:
            wedge.set_alpha(0.2)
        pie[0][-1].set_visible(False)
    else:
        _ = plt.pie([coord_same,
                     coord_diff,
                     coord_none,
                     n_bilat_teprn],
                    labels=['bilateral T&E coordination same interval', 'bilateral T&E coordination different interval',
                            'bilateral T&E no coordination', 'T&E PRN'],
                    # colors=['deepskyblue', 'turquoise', 'green', 'red'],
                    colors=[

                        "#14185c",
                        "#ffa600",
                        "#e9484a",
                        "#951269", ],
                    labeldistance=None,
                    **kwargs)

    plt.title(f'{name}')


plt.subplots(2, 2, figsize=(10, 10))
for i, (name, data) in enumerate(name_data_iter):
    print("=" * linelength)
    print(f"Data {diagnoses_name_mapping[name]} pie diagram")
    print("=" * linelength)
    plt.subplot(2, 2, i + 1)
    pie_diagram(data, diagnoses_name_mapping[name])
# plt.suptitle('Pie diagrams of different treatment')
plt.legend(loc='lower center', bbox_to_anchor=(0.0, -0.2), ncol=2, frameon=False)
filename = f'pie_diagrams.png'
output_file = Path(f'plots/pie1') / filename
output_file.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_file)
plt.show()

# split by coordination same or different interval
# and count n_bilateral t&e (times two) for coordination base
coord_bases = [1, 2, 3, 4, 5]
coord_reasons_2o3 = [1, 2, 3, 4, 5, 6]

coords_all_cols = {str(cbase) if coord_bases not in [2, 3] else f'{cbase}_{creason}'
                   for cbase in coord_bases for creason in coord_reasons_2o3}
coords_all_cols = sorted(coords_all_cols)

df_coord = pd.DataFrame()
for name, data in name_data_iter:
    for coord_type in ['same', 'different']:
        rowname = f'{name}_{coord_type}_interval'
        for coord_base in coord_bases:
            colname = f'coordbase{coord_base}'
            if coord_base in [2, 3]:
                for coord_cause in coord_reasons_2o3:
                    colname_add = f'{colname}_{coord_cause}'
                    df_coord.loc[rowname, colname_add] = 0.
            else:
                df_coord.loc[rowname, colname] = 0.

for name, data in name_data_iter:
    for coord_base in coord_bases:
        colname = f'coordbase{coord_base}'
        base_query = f'coordination_base == {coord_base}'

        for coord_type in ['same', 'different']:
            df_sel = data.query(f"coordination_{coord_type}_intervall == 1 & {base_query}")
            rowname = f'{name}_{coord_type}_interval'
            if coord_base in [2, 3]:
                for row in df_sel.itertuples():
                    coord_causes = []
                    for x in str(row.coordination_cause).split(','):
                        try:
                            x_int = int(x)
                        except ValueError:
                            x_int = int(float(x))
                        coord_causes.append(x_int)
                    n_coord_causes = len(coord_causes)
                    for coord_cause in coord_causes:
                        colname_add = f'{colname}_{coord_cause}'

                        # if colname_add not in df_coord.columns:
                        #     df_coord[colname_add] = 0
                        # if rowname not in df_coord.index:
                        #     df_coord.loc[rowname] = 0

                        df_coord.loc[rowname, colname_add] += row.n_bilateral_txe * 2 / n_coord_causes
                    # print(type(row.coordination_cause), row.coordination_cause)
            else:
                df_coord.loc[rowname, colname] = df_sel['n_bilateral_txe'].sum() * 2  # two eyes

        # for coord_type in ['same', 'different']:
        #     queried_data = data.query(f"coordination_{coord_type}_intervall == 1 & {base_query}")

df_coord_out = "Coordination types, n_bilateral_txe * 2 \n" \
               "============================================\n"
df_coord_out += str(df_coord)
outfile = Path('outputs/df_coord')
outfile.parent.mkdir(parents=True, exist_ok=True)
with open(outfile.with_suffix('.txt'), 'w') as f:
    f.write(df_coord_out)
df_coord.to_excel(outfile.with_suffix('.xlsx'), sheet_name='n_bilateral_txe')
print(df_coord_out)

# get number of bilateral, unilateral and n_patients for all four diagnoses and the two medication types
cfg_general = config['data_general']
df_general = pd.DataFrame(index=['bilat', 'unilat', 'n_patients'])
medications = ['eylea', 'lucentis']
for medication in medications:
    df_general[medication] = 0
for name, data in name_data_iter:
    df_general.loc['bilat', name] = data['n_bilateral'].sum() * 2
    df_general.loc['unilat', name] = data['n_unilateral'].sum()
    df_general.loc['n_patients', name] = data.shape[0]
    for i, medication in enumerate(medications):
        df_general.loc['bilat', medication] += data[f'n_bilateral_{medication}'].sum() * 2
        df_general.loc['bilat', medication] += data[f'n_bilateral_lucentis_eylea'].sum()  # each on for each eye
        df_general.loc['unilat', medication] += data[f'n_unilateral_{medication}'].sum()
        df_general.loc['n_patients', medication] += \
            data.query(f'n_{medication}_total >= n_{medications[(i + 1) % 2]}_total').shape[0]

df_general.rename(index=cfg_general['indexnames'], columns=cfg_general['columnnames'], inplace=True)
df_general.rename(columns=diagnoses_name_mapping, inplace=True)
df_general_out = str(df_general)
print(df_general_out)
outfile = Path('outputs/df_general')
outfile.parent.mkdir(parents=True, exist_ok=True)
with open(outfile.with_suffix('.txt'), 'w') as f:
    f.write(df_general_out)
df_general.to_excel(outfile.with_suffix('.xlsx'), sheet_name=cfg_general['sheetname'])

# adversarial effects
lateralities = ['bilateral', 'unilateral']

for laterality in lateralities:
    for name, data in name_data_iter:
        for i, adverse_col in enumerate(['hyposphagma', 'sicca', 'allergy', 'iod']):
            has_adverse = (data[f'{adverse_col}_{laterality}'].astype(str) != "0")
            data[f'ocular_ae_{laterality}'] = data[f'ocular_ae_{laterality}'].astype(str) + " " + (
                    has_adverse * (19 + i)).astype(str)

        # print(data[f'ocular_ae_{laterality}'])
        # preprocess data, add the two adverse effects

dfs = {'ocular': None, 'systemic': None}
neffects = {'ocular': 22, 'systemic': 24}

# make number to effect mapping
effect_mapping = {}

for adverse_type in ['ocular', 'systemic']:
    df_adverse = pd.DataFrame(index=lateralities, columns=[f'effect_{i + 1}' for i in range(neffects[adverse_type])])
    df_adverse.fillna(0, inplace=True)
    dfs[adverse_type] = df_adverse


    for name, data in name_data_iter:
        for laterality in lateralities:
            adv_counter = Counter()
            for row in data.itertuples():
                effects = str(getattr(row, f'{adverse_type}_ae_{laterality}'))
                effects = effects.replace(',', ' ').split(' ')
                effects = [int(x) for x in effects if x not in ['', '0']]
                row_counter = Counter(effects)
                adv_counter += row_counter

            for effect, count in adv_counter.items():
                df_adverse.loc[laterality, f'effect_{effect}'] += count

rename_mapping_ocular = config['adversenames']['ocular']
dfs['ocular'].rename(columns=rename_mapping_ocular, inplace=True)
rename_mapping_systemic = config['adversenames']['systemic']
dfs['systemic'].rename(columns=rename_mapping_systemic, inplace=True)
for df in dfs.values():
    df['total'] = df.sum(axis=1)

dfs_normalized = {}
for adverse_type in ['ocular', 'systemic']:
    df = dfs[adverse_type].copy()
    # df = df.div(df.sum(axis=1), axis=0).fillna(0) * 100  # percentage
    df = df.div(n_tot_inj, axis=0).fillna(0) * 100  # percentage
    dfs_normalized[adverse_type] = df

for dfs_use, normalized in zip([dfs, dfs_normalized], [False, True]):

    out = f"Adverse effects {'percent' if normalized else ''}\n" \
          "============================\n"
    for adverse_type in ['ocular', 'systemic']:
        out += f"{adverse_type} effects (number of occurences) \n"
        out += str(dfs_use[adverse_type])
        out += "\n____________________\n"
    print(out)

    outfile = Path(f'outputs/df_adverse{"_percent" if normalized else ""}')
    outfile.parent.mkdir(parents=True, exist_ok=True)
    with open(outfile.with_suffix('.txt'), 'w') as f:
        f.write(out)

    for adv_type, df in dfs_use.items():
        for transpose in True, False:
            file = Path(str(outfile) + f'_{adv_type}{"_transposed" if transpose else ""}')
            df_tmp = df
            if transpose:
                df_tmp = df.transpose()
            df_tmp.to_excel(file.with_suffix('.xlsx'), sheet_name=adv_type)

# dump summary statistics of all four datasets

columns = [
    'avg_age',
    'n_patients',
    'n_injections',
    'male_female_ratio_percent',
    'n_unilateral',
    'n_bilateral',
    'delta_bcva',
    'n_lucentis_total',
    'n_lucentis_unilateral',
    'n_lucentis_bilateral',
    'n_eylea_total',
    'n_eylea_unilateral',
    'n_eylea_bilateral',
]

df_summary = pd.DataFrame(index=diagnoses, columns=columns)
df_summary.fillna(0, inplace=True)
for name, data in name_data_iter:
    df_summary.loc[name, 'avg_age'] = round(data['age'].mean(), 1)
    n_patients = data.shape[0]
    df_summary.loc[name, 'n_patients'] = n_patients
    df_summary.loc[name, 'n_injections'] = data['n_total'].sum()
    # female is 1, male is 2
    # TODO: change to nmale (perc_nmale)/ nfemale (perc_nfemale)
    df_summary.loc[name, 'male_female_ratio_percent'] = round(np.sum(data['sex'] == 2) / np.sum(data['sex'] == 1) * 100,
                                                              0)
    df_summary.loc[name, 'n_unilateral'] = data['n_unilateral'].sum()
    df_summary.loc[name, 'n_bilateral'] = data['n_bilateral'].sum()
    for medication in medications:
        df_summary.loc[name, f'n_{medication}_total'] = data[f'n_{medication}_total'].sum()
        df_summary.loc[name, f'n_{medication}_unilateral'] = data[f'n_unilateral_{medication}'].sum()
        df_summary.loc[name, f'n_{medication}_bilateral'] = data[f'n_bilateral_{medication}'].sum()
    df_summary.loc[name, 'delta_bcva'] = round(data['bcva_sum_diff'].mean(), 3)

df_summary.rename(columns=config['data_general']['columnnames'], index=diagnoses_name_mapping, inplace=True)
# df_summary.rename(, inplace=True)
df_summary_out = str(df_summary)
print(df_summary_out)
outfile = Path('outputs/df_summary')
outfile.parent.mkdir(parents=True, exist_ok=True)
with open(outfile.with_suffix('.txt'), 'w') as f:
    f.write(df_summary_out)

df_summary.to_excel(outfile.with_suffix('.xlsx'), sheet_name='summary')
