from collections import Counter
from pathlib import Path

import hist
import matplotlib.pyplot as plt
import mplhep
import numpy as np
import pandas as pd
import scipy.stats
import yaml

# load config yaml
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)
# manual config, change digits of percent if needed
perc_formatstr = '{:#.2g}'.format

diagnoses = ['AMD', "DMÖ", 'ZVV_VAV', 'multiple_diagnoses']
lateralities = ['bilateral', 'unilateral']
# diagnoses.reverse()
datasets = [pd.read_excel("data_small.xlsx", sheet_name=name).reset_index(drop=True) for name in diagnoses]
pd.set_option('display.max_columns', None, 'display.max_rows', None, 'display.expand_frame_repr', False)
diagnoses_name_mapping = config['diagnosesnames']
coordnames = config['coordnames']

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
    out = f'{data.mean() * factor:#.3g} +- {data.std() * factor / data.shape[0] ** 0.5:#.2g}'
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
n_col_all = {k: {k2: f'{v:#.4g} +- {n_col_std_all[k][k2]:#.2g}' for k2, v in v.items()} for k, v in n_col_all.items()}

print("Injections by med_lat")
injections_med_lat = pd.DataFrame({k: v.values()
                                   for k, v in n_col_all.items()}, index=list(n_col_all.values())[0].keys())
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
n_col_all = {k: {k2: f'{v:#.3g} +- {n_col_std_all[k][k2]:#.2g}' for k2, v in v.items()} for k, v in n_col_all.items()}
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
    print("mul etc are datasets, index is the encoding")
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

# lateral = ['uni', 'bi']
# index = pd.MultiIndex.from_product([diagnoses, medications, lateral])
# adverse_eff = ['hyposphagma', 'sicca', 'allergy']
# adverse_hypo = {
#     1: "unilat_non-antikoag",
#     2: "unilat_antikoag",
#     3: "bilat_non-antikoag",
#     4: "bilat_antikoag"}
# shifts_med = {'eylea': 4, 'lucentis': 0}
# # shifts_lat = {'uni': 0, 'bi': 2}
# results = pd.DataFrame(index=index, columns=adverse_hypo.values())
# results_abs = pd.DataFrame(index=index, columns=adverse_hypo.values())
# for name, data in name_data_iter:
#     print("=" * linelength)
#     print(f"Data {diagnoses_name_mapping[name]} adverse effects")
#     print("=" * linelength)
#     n_collected = {}
#     n_collected_std = {}
#
#     df = data.copy()
#     for medication in medications:
#         for lat in lateral:
#             key = f'_{medication[:3]}_{lat}'
#             # shift_lat = shifts_lat[lat]
#             shift = shifts_med[medication]
#             for i, adverse_name in adverse_hypo.items():
#                 series_sel = df[f'hyposphagma_{lat}lateral']
#                 series_sel = series_sel.apply(lambda x: x[0] if isinstance(x, (list, tuple)) else x)
#                 if not series_sel.dtype == 'int64':
#                     pass
#                 have_it = (series_sel == i + shift)
#                 n_have_it = have_it.sum()
#                 ntot = df["n_total"].sum()
#                 percent_val = n_have_it / ntot
#                 results.loc[
#                     (name, medication, lat), adverse_name] = f'{percent_val:#.1%} +- {n_have_it ** 0.5 / ntot:.1%}'
#                 results_abs.loc[(name, medication, lat), adverse_name] = f'{n_have_it} +- {n_have_it ** 0.5:.1f}'
#
# print("Adverse effects per data in percent of total injections")
# print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
#       "\n(medication, uni/bilateral injection taken from table name/encoding)")
# print(results)
#
# print("Adverse effects per data in absolute numbers")
# print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
#       "\n(medication, uni/bilateral injection taken from table name/encoding)")
# print(results_abs)

print("\n\n\n")

lateral = ['uni', 'bi']
index = pd.MultiIndex.from_product([diagnoses, medications, lateral])
adverse_eff = {'hyposphagma': 'hyposphagma_{lat}lateral > 0',
               'sicca': 'sicca_{lat}lateral > 0',
               # 'allergy': 'allergy_{lat}lateral > 0'
               }
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
                series_sel = series_sel.apply(lambda x: i if (isinstance(x, (list, tuple)) and i in x) else x)
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
# save results
adversedir = Path('outputs/adverse_effects')
adversedir.mkdir(parents=True, exist_ok=True)
outfile_hyposphagma = adversedir / 'hyposphagma'
out = "Adverse effects, per data in percent of total injections \n" \
        "====================================================\n"
out += str(results_abs)

outfile_hyposphagma.with_suffix('.txt').write_text(out)

results_abs.to_excel(outfile_hyposphagma.with_suffix('.xlsx'), sheet_name='hyposphagma')


# add subcategories for hyposhagma and IOD. Looks similar to above, but is different
lateral = ['uni', 'bi']
index = pd.MultiIndex.from_product([diagnoses, medications, lateral])
adverse_eff = {'hyposphagma': 'hyposphagma_{lat}lateral > 0',
               # 'sicca': 'sicca_{lat}lateral > 0',
               # 'allergy': 'allergy_{lat}lateral > 0',
               'iod': 'iod_{lat}lateral > 0',
               }

hypo_iods = ['hypo_antikoag', 'hypo_nonanti', 'IOD_no_Glaukom', 'IOD_with_Glaukom']
allergy_causes =  [f'allergy_{i}' for i in range(1, 7)]
columns = [f'{eff}_{lat}' for eff in hypo_iods + allergy_causes for lat in lateralities]
index = diagnoses
results2 = pd.DataFrame(index=index, columns=columns).fillna(0.)
results_abs2 = pd.DataFrame(index=index, columns=columns).fillna(0.)
for name, data in name_data_iter:
    print("=" * linelength)
    print(f"Data {diagnoses_name_mapping[name]} adverse effects")
    print("=" * linelength)
    n_collected = {}
    n_collected_std = {}

    df = data.copy()
    colname = '{adverse}_{lat}'
    for i, (eff, query) in enumerate(adverse_eff.items()):
        # df = data.query(query)
        evenname = hypo_iods[(i * 2)]
        oddname = hypo_iods[(i * 2) + 1]

        for j in range(10):
            if j == 0:
                # for lat in lateralities:
                #     results_abs2.loc[name, colname.format(adverse=evenname, lat=lat)] = 0
                #     results_abs2.loc[name, colname.format(adverse=oddname, lat=lat)] = 0
                #     results2.loc[name, colname.format(adverse=evenname, lat=lat)] = 0
                #     results2.loc[name, colname.format(adverse=oddname, lat=lat)] = 0
                continue
            if j % 2 == 0:
                adverse_name = evenname
            else:
                adverse_name = oddname



            for lat in lateralities:
                series_sel = df[f"{eff}_{lat}"]
                have_it = series_sel.apply(lambda x: j in x if isinstance(x, (list, tuple)) else j == x)
                # if not series_sel.dtype == 'int64':
                #     pass

                # have_it = (series_sel == i + shift)
                n_have_it = have_it.sum()
                ntot = df["n_total"].sum()
                percent_val = n_have_it / ntot

                # results2.loc[name, adverse_name] += f'{percent_val:.1%} +- {n_have_it ** 0.5 / ntot:.1%}'
                results_abs2.loc[name, colname.format(adverse=adverse_name, lat=lat)] += n_have_it
                results2.loc[name, colname.format(adverse=adverse_name, lat=lat)] += percent_val  * 100
    for i, adverse_name in enumerate(allergy_causes):
        i += 1

        for lat in lateralities:
            series_sel = df[f'allergy_cause_{lat}']
            have_it = series_sel.apply(lambda x: i in x if isinstance(x, (list, tuple)) else i == x)
            n_have_it = have_it.sum()
            ntot = df["n_total"].sum()
            percent_val = n_have_it / ntot
            results_abs2.loc[name, colname.format(adverse=adverse_name, lat=lat)] += n_have_it
            results2.loc[name, colname.format(adverse=adverse_name, lat=lat)] += percent_val * 100



rename_mapping = config['adversenames']['ocular']
rename_mapping_allergy = config['adversenames']['allergy']
rename_mapping.update(rename_mapping_allergy)
results_abs2.rename(columns=rename_mapping, inplace=True)
results2.rename(columns=rename_mapping, inplace=True)


print("Adverse effects per data in percent of total injections, percent, antikoag and stuff together")
print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
      "\n(medication, uni/bilateral injection taken from table name/encoding)")
print(results2)

print("Adverse effects per data in absolute numbers")
print("Columns: adverse effects.\nIndex: split by data, medication, uni/bilateral injection"
      "\n(medication, uni/bilateral injection taken from table name/encoding)")
print(results_abs2)
# save results
adversedir = Path('outputs/adverse_effects')
adversedir.mkdir(parents=True, exist_ok=True)
outfile_hyposphagma = adversedir / 'hypo_iod_summary'
out = "Adverse effects, per data in percent of total injections \n" \
      "====================================================\n"
out += str(results_abs2)

outfile_hyposphagma.with_suffix('.txt').write_text(out)

results_abs2.to_excel(outfile_hyposphagma.with_suffix('.xlsx'), sheet_name='hyposphagma')


# sys.exit(0)







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

# split by coordination same or different interval
# and count n_bilateral t&e (times two) for coordination base
coord_bases = [1, 2, 3, 4, 5]
coord_reasons_2o3 = [1, 2, 3, 4, 5, 6]

coords_all_cols = {str(cbase) if coord_bases not in [2, 3] else f'{cbase}_{creason}'
                   for cbase in coord_bases for creason in coord_reasons_2o3}
coords_all_cols = sorted(coords_all_cols)

dfs_coordbase = [pd.DataFrame() for _ in name_data_iter]

df_coord = pd.DataFrame()
df_coord_all = pd.DataFrame()
for name, data in name_data_iter:
    for coord_type in ['same', 'different']:
        rowname_nodata = f'{coord_type}_interval'
        rowname = f'{name}_{rowname_nodata}'
        for coord_base in coord_bases:
            colname = f'coordbase{coord_base}'
            if coord_base in [2, 3]:
                for coord_cause in coord_reasons_2o3:
                    colname_add = f'{colname}_{coord_cause}'
                    df_coord.loc[rowname, colname_add] = 0.
                    df_coord_all.loc[rowname_nodata, colname_add] = 0.
            else:
                df_coord.loc[rowname, colname] = 0.
                df_coord_all.loc[rowname_nodata, colname] = 0.

for name, data in name_data_iter:
    for coord_base in coord_bases:
        colname = f'coordbase{coord_base}'
        base_query = f'coordination_base == {coord_base}'

        for coord_type in ['same', 'different', 'none']:
            if coord_type == 'none':
                df_sel = data.query(f"no_coordination == 1 & {base_query}")
            else:
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
df_coord.loc[:, 'coordbase4'] += df_coord.pop('coordbase2_6')

for name, _ in name_data_iter:
    for inter in ['same', 'different']:
        df_coord_all.loc[f'{inter}_interval'] = 0.

for name, _ in name_data_iter:
    for inter in ['same', 'different']:
        df_coord_all.loc[f'{inter}_interval'] += df_coord.loc[f'{name}_{inter}_interval']

df_coord_all

# for coord_type in ['same', 'different']:
#     queried_data = data.query(f"coordination_{coord_type}_intervall == 1 & {base_query}")

df_coord_out = "Coordination types, n_bilateral_txe * 2 \n" \
               "============================================\n"
df_coord_out += str(df_coord)
adversedir = Path('outputs/df_coord')
adversedir.parent.mkdir(parents=True, exist_ok=True)
with open(adversedir.with_suffix('.txt'), 'w') as f:
    f.write(df_coord_out)
df_coord.to_excel(adversedir.with_suffix('.xlsx'), sheet_name='n_bilateral_txe')
print(df_coord_out)

df_coord_out = "Coordination types merged, n_bilateral_txe * 2 \n" \
               "============================================\n"
df_coord_out += str(df_coord_all)
adversedir = Path('outputs/df_coord_merged')
adversedir.parent.mkdir(parents=True, exist_ok=True)
with open(adversedir.with_suffix('.txt'), 'w') as f:
    f.write(df_coord_out)
df_coord_all.to_excel(adversedir.with_suffix('.xlsx'), sheet_name='n_bilateral_txe')
print(df_coord_out)


# colorsdiag = [
#     # "#14185c",
#     "#555eff",
#     "#ffa600",
#     "#e9484a",
#     "#951269", ]
colorsdiag = plt.cm.gray([ 0.4, 0.7, 0.55, 0.9,])

# hatches = ['/', '\\', '|', '-', '+', 'x', 'o', 'O', '.', '*']
# hatches = [
#     '/',
#     '\\',
#     '+',
#     'o',
# ]
hatchesbar = [
    '/',
    '\\',
    '+',
    'o',
]
hatches = None

def pie_diagram(data, name: str, *,beta=0):
    """Pie diagram of the different treatment (coord, bilat t&e, coord diff and no coord)."""
    n_bilat_te = data['n_bilateral_txe'].sum() * 2  # 2 injections per person
    n_bilat_teprn = data['n_bilateral_txe_prn'].sum() * 2  # 2 injections per person

    coord_same = data.query('coordination_same_intervall == 1')['n_bilateral_txe'].sum() * 2
    coord_diff = data.query('coordination_different_intervall == 1')['n_bilateral_txe'].sum() * 2
    coord_none = data.query('no_coordination == 1')['n_bilateral_txe'].sum() * 2
    # done: horizontal bar diagrams below that show the causes of the different coordination types

    all_coord = coord_same + coord_diff + coord_none
    # hack, what's the right number?
    n_bilat_te = all_coord
    assert all_coord == n_bilat_te, f'{all_coord} != {n_bilat_te}'

    # plt.figure(figsize=(10, 10))

    kwargs = dict(autopct=lambda x: f'{int(x * all_coord / 100)} ({x / 100:.1%})', textprops={'fontsize': 15},
                  counterclock=True, startangle=0, # 90 + beta,
                  pctdistance=0.7)
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

        datatoplot = [
            coord_same,
                  coord_diff,
            n_bilat_teprn,
                      coord_none,
                      ]
        labelstoplot = [
            'Equal TER',
            'Coordinated TER',
            'Mixed TER and PRN',
            'Async TER',
        ]
        colorstoplot = colorsdiag
        datacleaned = []
        labelscleaned = []
        colorscleaned = []
        for i, d in enumerate(datatoplot):
            if d > 0:
                datacleaned.append(d)
                labelscleaned.append(labelstoplot[i])
                colorscleaned.append(colorstoplot[i])
        _ = plt.pie(datacleaned,
                    labels=labelscleaned,
                    # colors=['deepskyblue', 'turquoise', 'green', 'red'],
                    colors=colorscleaned,
                    labeldistance=None,
                    hatch=hatches,
                    radius=1,
                    **kwargs)

    plt.title(f'{name}', y=0.95, fontsize = 20, pad=0.0)


cols = 2
rows = 4
plotnrpie = {0: 1, 1: 2, 2: 5, 3: 6}
plotnrbar = {0: 3, 1: 4, 2: 7, 3: 8}

plt.subplots(rows, cols, figsize=(20, 20), height_ratios=[7, 8 / 3, 7, 1])



def plot_barh(data, coordcolnames: list[str], name: str, invert=False):
    # n_bilat_teprn = data['n_bilateral_txe_prn'].sum() * 2  # 2 injections per person

    coordlabels = [coordnames[col] for col in coordcolnames]
    left = np.zeros((len(coordcolnames),))

    coords = ['same', 'different', 'none']
    for coord, color, hatch in zip(coords, colorsdiag[:-1], hatchesbar):
        diff = [data.loc[f'{name}_{coord}_interval', col] for col in coordcolnames]

        plt.barh(coordlabels, left + diff, color=color, left=left, hatch=hatch)
        if invert:
            # plt.gca().yaxis.tick_right()
            plt.gca().tick_params(axis='both', which='major', labelsize=6)
            plt.gca().invert_xaxis()
        left += diff


coordcolnames = [
    "coordbase4",
    "coordbase2_1",
    "coordbase2_2",
    "coordbase2_4",
    "coordbase3_2",
    "coordbase3_3",
    "coordbase3_5",
    "coordbase1",
]
coordcolnames.reverse()
coordcols = coordcolnames

# Bar plot with pia diagrams, uncomment to plot
# for i, (name, data) in enumerate(name_data_iter):
#     if i >= 2:
#         coordcols = [n for i, n in enumerate(coordcolnames) if i in (0, 5, 7)]
#     print("=" * linelength)
#     print(f"Data {diagnoses_name_mapping[name]} pie diagram")
#     print("=" * linelength)
#     plt.subplot(rows, cols, plotnrpie[i])
#     pie_diagram(data, diagnoses_name_mapping[name])
#
#     plt.subplot(rows, cols, (ploti := plotnrbar[i]))
#     if ploti in (4, 8):  # right hbar, no labels
#         plt.gca().set_yticklabels([])
#     plot_barh(df_coord, coordcols, name, invert=i % 2 == 0)
# # plt.suptitle('Pie diagrams of different treatment')
#
# plt.subplot(rows, cols, 6)
# plt.legend(
#     loc='lower center',
#     bbox_to_anchor=(-0.15, -0.5),
#     ncol=2, frameon=False,)
# filename = 'pie_diagrams_bars.pdf'
# output_file = Path('plots/pie1') / filename
# output_file.parent.mkdir(parents=True, exist_ok=True)
# plt.savefig(output_file)

# only bar plot
plt.subplots(2, 2, figsize=(20, 20))
for i in range(4):  # TODO
    # if i >= 2:
    coordcols = [n for i, n in enumerate(coordcolnames) if i in (1, 3)]
    plt.subplot(rows, cols, i + 1)

    if i in (1, 3):  # right hbar, no labels
        plt.gca().set_yticklabels([])
    plot_barh(df_coord, coordcols, name, invert=i % 2 == 0)
# plt.suptitle('Pie diagrams of different treatment')

plt.subplot(rows, cols, 6)
plt.legend(
    loc='lower center',
    bbox_to_anchor=(-0.15, -0.5),
    ncol=2, frameon=False,)
filename = 'coord_reason_bars.pdf'
output_file = Path('plots/bars') / filename
output_file.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_file)

# only pie diagrams
ncols = 2
nrows = 2
plt.subplots(nrows, ncols, figsize=(20, 20))
for i, (name, data) in enumerate(name_data_iter):
    plt.subplot(nrows, ncols, i + 1)
    pie_diagram(data, diagnoses_name_mapping[name], beta=-27 if i == 1 else 0)

    # plt.subplot(rows, cols, (ploti := plotnrbar[i]))
    # if ploti in (4, 8):  # right hbar, no labels
    #     plt.gca().set_yticklabels([])
    # plot_barh(df_coord, coordcols, name, invert=i % 2 == 0)
# plt.suptitle('Pie diagrams of different treatment')

plt.subplot(nrows, ncols, 3)
plt.tight_layout()
plt.legend(
    loc='upper right',
    bbox_to_anchor=(1.2, 1.1),
    # ncol=2,
    frameon=False,
    fontsize=20,
    borderaxespad=0.0,
    borderpad=0.0,
)
filename = 'pie_diagrams.pdf'
output_file = Path('plots/pie1') / filename
output_file.parent.mkdir(parents=True, exist_ok=True)
plt.savefig(output_file)
# plt.show()

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
        df_general.loc['bilat', medication] += data['n_bilateral_lucentis_eylea'].sum()  # each on for each eye
        df_general.loc['unilat', medication] += data[f'n_unilateral_{medication}'].sum()
        df_general.loc['n_patients', medication] += \
            data.query(f'n_{medication}_total >= n_{medications[(i + 1) % 2]}_total').shape[0]

df_general.rename(index=cfg_general['indexnames'], columns=cfg_general['columnnames'], inplace=True)
df_general.rename(columns=diagnoses_name_mapping, inplace=True)
df_general_out = str(df_general)
print(df_general_out)
adversedir = Path('outputs/df_general')
adversedir.parent.mkdir(parents=True, exist_ok=True)
with open(adversedir.with_suffix('.txt'), 'w') as f:
    f.write(df_general_out)
df_general.to_excel(adversedir.with_suffix('.xlsx'), sheet_name=cfg_general['sheetname'])

# adversarial effects


for laterality in lateralities:
    for name, data in name_data_iter:
        for i, adverse_col in enumerate(['hyposphagma', 'sicca', 'allergy', 'iod']):
            has_adverse = (data[f'{adverse_col}_{laterality}'].astype(str) != "0")
            data[f'ocular_ae_{laterality}_encoded'] = data[f'ocular_ae_{laterality}'].astype(str) + " " + (
                    has_adverse * (19 + i)).astype(str)

        print(data[f'ocular_ae_{laterality}'])
        # preprocess data, add the two adverse effects

dfs = {'ocular': None, 'systemic': None}
neffects = {'ocular': 22, 'systemic': 24}

# make number to effect mapping
effect_mapping = {}

def n_adverse_from_list_of_types(x):
    print(x)
    return sum([int(y) > 0 if y != '' else 0 for y in str(x).replace(' ', '').split(',')])


for adverse_type in ['ocular', 'systemic']:
    df_adverse = pd.DataFrame(
        index=lateralities,
        columns=[f"effect_{i + 1}" for i in range(neffects[adverse_type])]
        + [f"effect_{i}" for i in range(100, 104)]  # hyposphagma, IOD
        + allergy_causes
    )
    df_adverse.fillna(0., inplace=True)

    dfs[adverse_type] = df_adverse


    for name, data in name_data_iter:
        for lat in lateralities:
            # count the number of adverse events
            data[f'n_{adverse_type}_ae_{lat}'] = data[f'{adverse_type}_ae_{lat}'].apply(lambda x: n_adverse_from_list_of_types(x))
            adv_counter = Counter()
            for row in data.itertuples():
                effects = str(getattr(row, f'{adverse_type}_ae_{lat}'))
                effects = effects.replace(' ', '').rstrip(',').split(',')
                effects = [int(x) for x in effects if x not in ['', '0']]
                row_counter = Counter(effects)
                adv_counter += row_counter

            for effect, count in adv_counter.items():
                df_adverse.loc[lat, f'effect_{effect}'] += count
            # add 'IOL dezentriert/HKL subluxiert', effect 15 and 18
            if adverse_type == 'ocular':
                for i, effect in enumerate(hypo_iods):
                    df_adverse.loc[lat, f"effect_10{i}"] += results_abs2.loc[name, f'{effect}_{lat}']
                for i, cause in enumerate(allergy_causes):
                    df_adverse.loc[lat, f"allergy_{i + 1}"] += results_abs2.loc[name, f'{cause}_{lat}']
    if adverse_type == 'systemic':
        df_adverse.drop(columns=[f'effect_10{i}' for i in range(4)], inplace=True)
        df_adverse.drop(columns=[f'allergy_{i}' for i in range(1,7)], inplace=True)
    # add hyposphagma non-antikoag, antikoag and iod glaukom, non-glaukom
    df_adverse.loc[:, 'effect_15'] += df_adverse.loc[:, 'effect_18']
    df_adverse.drop(columns='effect_18', inplace=True)
    # remove effect 6, Charles Bonnet
    df_adverse.drop(columns='effect_6', inplace=True)



for adverse_type in ['ocular', 'systemic']:
    rename_mapping = config['adversenames'][adverse_type]
    dfs[adverse_type].rename(columns=rename_mapping, inplace=True)
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

    adversedir = Path(f'outputs/df_adverse{"_percent" if normalized else ""}')
    adversedir.parent.mkdir(parents=True, exist_ok=True)
    with open(adversedir.with_suffix('.txt'), 'w') as f:
        f.write(out)

    for adv_type, df in dfs_use.items():
        for transpose in True, False:
            file = Path(str(adversedir) + f'_{adv_type}{"_transposed" if transpose else ""}')
            df_tmp = df
            if normalized:

                df_tmp = df_tmp.map(perc_formatstr)
            if transpose:
                df_tmp = df.transpose()
            df_tmp.to_excel(file.with_suffix('.xlsx'), sheet_name=adv_type)


# dump summary statistics of all four datasets

columns = [
    'avg_age',
    'n_patients',
    'n_injections',
    'male_female_percent',
    'n_unilateral',
    'n_bilateral',
    'bcva_baseline_r_mean',
    'bcva_baseline_l_mean',
    'delta_bcva',
    'n_lucentis_total',
    'n_lucentis_unilateral',
    'n_lucentis_bilateral',
    'n_eylea_total',
    'n_eylea_unilateral',
    'n_eylea_bilateral',
    'n_cat_surgery',
    'n_argon_surgery',
    'n_kapsulotomie_surgery',
]

df_summary = pd.DataFrame(index=diagnoses, columns=columns)
df_summary.fillna(0, inplace=True)
for name, data in name_data_iter:
    df_summary.loc[name, 'avg_age'] = round(data['age'].mean(), 1)
    n_patients = data.shape[0]
    df_summary.loc[name, 'n_patients'] = n_patients
    df_summary.loc[name, 'n_injections'] = data['n_total'].sum()
    # female is 1, male is 2
    # done: change to nmale (perc_nmale)/ nfemale (perc_nfemale)
    nmale = np.sum(data['sex'] == 2)
    nfemale = np.sum(data['sex'] == 1)
    ntot = nmale + nfemale
    female_male_string = f'{nmale} ({nmale / ntot:.1%}) / {nfemale} ({nfemale / ntot:.1%})'
    df_summary.loc[name, 'male_female_percent'] = female_male_string
    df_summary.loc[name, 'n_unilateral'] = data['n_unilateral'].sum()
    df_summary.loc[name, 'n_bilateral'] = data['n_bilateral'].sum()
    df_summary.loc[name, 'bcva_baseline_r_mean'] = round(data['bcva_baseline_r'].mean(), 3)
    df_summary.loc[name, 'bcva_baseline_l_mean'] = round(data['bcva_baseline_l'].mean(), 3)
    df_summary.loc[name, 'bcva_baseline_mean'] = round((data['bcva_baseline_r'].mean() + data['bcva_baseline_l'].mean())/2, 3)

    # cat surgery
    cat_surgery_map = {0: 0, 1: 2, 2: 1, 3: 1}
    df_summary.loc[name, 'n_cat_surgery'] = data['cat_surgery'].apply(lambda x: cat_surgery_map[x]).sum()
    df_summary.loc[name, 'n_argon_surgery'] = data['argon'].apply(lambda x: cat_surgery_map[x]).sum()
    df_summary.loc[name, 'n_kapsulotomie_surgery'] = data['kapsulotomie'].apply(lambda x: cat_surgery_map[x]).sum()
    for medication in medications:
        df_summary.loc[name, f'n_{medication}_total'] = data[f'n_{medication}_total'].sum()
        df_summary.loc[name, f'n_{medication}_unilateral'] = data[f'n_unilateral_{medication}'].sum()
        df_summary.loc[name, f'n_{medication}_bilateral'] = data[f'n_bilateral_{medication}'].sum()
    df_summary.loc[name, 'delta_bcva'] = round(data['bcva_sum_diff'].mean(), 3)

n_injections = pd.Series({'unilateral': df_summary['n_unilateral'].sum(),
                          'bilateral': df_summary['n_bilateral'].sum()})
df_summary.rename(columns=config['data_general']['columnnames'], index=diagnoses_name_mapping, inplace=True)
# df_summary.rename(, inplace=True)
df_summary_out = str(df_summary)
print(df_summary_out)
adversedir = Path('outputs/df_summary')
adversedir.parent.mkdir(parents=True, exist_ok=True)
with open(adversedir.with_suffix('.txt'), 'w') as f:
    f.write(df_summary_out)

df_summary.to_excel(adversedir.with_suffix('.xlsx'), sheet_name='summary')

from scipy.stats import binom as binomsp

name_n_unilat = '# unilateral'
name_n_bilat = '# bilateral'

# compare
def calc_pval_fisher(row):

    n_unilat = row[name_n_unilat]
    n_bilat = row[name_n_bilat]
    n_unilat_ae = row['unilateral']
    n_bilat_ae = row['bilateral']
    pval = scipy.stats.fisher_exact([[n_unilat_ae, n_bilat_ae], [n_unilat - n_unilat_ae, n_bilat - n_bilat_ae]], alternative='two-sided')[1]
    # pval = scipy.stats.fisher_exact([[n_unilat_ae, n_unilat - n_unilat_ae], [n_bilat_ae, n_bilat - n_bilat_ae]])[1]
    return pval


for advtype in ['ocular', 'systemic']:
    df = dfs[advtype]
    n_unilat = n_injections['unilateral']
    p_unilat = df.loc['unilateral'] / n_unilat
    n_bilat = n_injections['bilateral'] * 2  # two eyes
    p_bilat = df.loc['bilateral'] / n_bilat
    ptot = (n_unilat * p_unilat + n_bilat * p_bilat) / (n_unilat + n_bilat)
    z = (p_unilat - p_bilat) / (ptot * (1 - ptot) * (1 / n_unilat + 1 / n_bilat)) ** 0.5
    z = pd.DataFrame({'p value': scipy.stats.norm.sf(np.abs(z)) *2}, #twosided
                     index=z.index)
    nadv = ptot.shape[0]
    zbilat_larger = np.where((p_bilat > p_unilat).array, z.values[:, 0] * 0.5,  # one sided test
                             1)

    dftmp = pd.DataFrame({'zstat p value bilat > unilat': zbilat_larger,
                          'zstat p value': z.values[:, 0],
                          'Prob. unilateral': p_unilat, 'Prob. bilateral': p_bilat, name_n_unilat: [n_unilat] * nadv,
                          name_n_bilat: [n_bilat] * nadv, 'Prob. total': ptot})


    dfadvtest = pd.concat([dftmp, df.loc['unilateral'], df.loc['bilateral']], axis=1)
    dfadvtest['p value fisher'] = dfadvtest.apply(calc_pval_fisher, axis='columns')
    dfadvtest['rel. incr. fisher to zstat in %'] = np.round((dfadvtest['p value fisher'] - dfadvtest['zstat p value']) / dfadvtest['p value fisher'], 2) * 100
    print(f"Adverse events {advtype} test statistics")
    print(dfadvtest)
    with open(Path('outputs') / f'df_adverse_{advtype}_teststatistic.txt', 'w') as f:
        f.write(str(dfadvtest))
    dfadvtest.to_excel(Path('outputs') / f'df_adverse_{advtype}_teststatistic.xlsx', sheet_name='teststatistic')


dataall = pd.concat(datasets, ignore_index=True).reset_index(drop=True)
nmax = 10  # TODO, how many?
ntot_ae = None
nmin = 0
npoints = nmax + 1
x = np.linspace(nmin, nmax, npoints)
def plot_data_pdf(ps, title, data):
    if not isinstance(ps, list):
        ps = [ps]
    global x, bins
    plt.figure()
    plt.title(title)
    nbins = 50
    h = hist.Hist.new.Reg(npoints, nmin -0.5, nmax + 0.5).Double()
    h.fill(data)
    # bins, edges, fig = plt.hist(data, bins=50, density=True, alpha=0.5)
    mplhep.histplot(h, histtype='errorbar', color='k', label='data', alpha=0.5, yerr=True)
    ntot_ae = np.sum(h.counts())
    for i, p in enumerate(ps):
        if i == 0:
            label = 'binom pmf'
            fmt = 'b-'
        elif i == 1:
            label = r'$ + 1 \sigma $'
            fmt = 'r--'
        elif i == 2:
            label = r'$ - 1 \sigma $'
            fmt = 'r--'
        elif i == 3:
            label = r'$ + 2 \sigma $'
            fmt = 'y-.'
        elif i == 4:
            label = r'$ - 2 \sigma $'
            fmt = 'y-.'
        else:
            raise ValueError("Too many ps")
        plt.plot(x, ntot_ae * binomsp.pmf(x, nmax, p), fmt, label=label)
    plt.xlabel('Number of adverse events')
    plt.ylabel('# patients')
    plt.legend()
    outpath = Path('plots') / 'stats' / 'fits'
    outpath.mkdir(parents=True, exist_ok=True)
    plt.savefig(outpath / f'{title.replace(" ", "_")}.png')
    plt.show()




def create_loss(data, n):
    def pdf(x, p):
        return np.maximum(binomsp.pmf(x, n, p), 1e-10)

    def nll(p):  # binned poisson loss
        return -np.sum(np.log(pdf(data, p)))

    nll.errordef = 0.5

    return nll

# TODO: mix ocular & systemic adverse events (one dataset)
datafit_notseldf = dataall[['n_ocular_ae_bilateral']]
datafit_notsel = dataall['n_ocular_ae_bilateral']

for sel in ['n_ocular_ae_bilateral < 30', ' n_ocular_ae_bilateral < 3']:
    datafit = datafit_notseldf.query(sel)['n_ocular_ae_bilateral']
    ntot_ae = datafit.sum()
    nll = create_loss(datafit, nmax)

    import zfit
    zfit.run.set_autograd_mode(False)
    zfit.run.set_graph_mode(False)

    p = zfit.Parameter("p", 0.02, 0, 1)
    plot_data_pdf(p, title="Before fit" + sel, data =datafit)
    minimizer = zfit.minimize.Minuit()
    result = minimizer.minimize(nll, p).update_params()
    result.hesse()
    result.errors()
    print(result)
    title_afterfit = "After fit" + sel
    plot_data_pdf(p, title=title_afterfit, data=datafit_notsel)
    pval = result.params[p]['value']
    upper = result.params[p]['errors']['upper']
    lower = result.params[p]['errors']['lower']
    plot_data_pdf([pval, pval + upper, pval + lower, pval + 2 * upper, pval + 2 * lower],
                  title=title_afterfit, data=datafit_notsel)

    asimov = binomsp.rvs(size=int(ntot_ae) *1, n=nmax, p=p)  # TODO: asimov?
    nll_asimov = create_loss(asimov[asimov < 3], nmax)
    result_asimov = minimizer.minimize(nll_asimov, p).update_params()
    result_asimov.hesse()
    result_asimov.errors()
    print(result_asimov)
    pval = result_asimov.params[p]['value']
    plot_data_pdf([pval, pval + result_asimov.params[p]['errors']['upper'],
                   pval + result_asimov.params[p]['errors']['lower']],
                  title="Asimov fit" + sel, data=asimov)

# TODO: check if higher chance that multiple adverse events (poisson?)

# TODO: correlation of diabetes + 2 more mit adverse events per patient


# Done below: unilateral vs bilateral: adverse events per injections (or per two injections?)
occular_occ = dfs['ocular'].sum(axis=1).to_frame("oculartotal")
nunilat = df_summary['# unilateral'].sum()
nbilat = df_summary['# bilateral'].sum() * 2
occular_occ.loc['bilateral', '# injections'] = nbilat
occular_occ.loc['unilateral', '# injections'] = nunilat

occular_occ['relsig_adv'] = occular_occ['oculartotal'] ** 0.5 / occular_occ['oculartotal']
occular_occ['rel_adv'] = occular_occ['oculartotal'] / occular_occ['# injections']
occular_occ['sig_adv'] = occular_occ['relsig_adv'] * occular_occ['rel_adv']

for lat in ['bilateral', 'unilateral']:
    rel_adv_per_bilat = occular_occ.loc[lat, 'rel_adv']
    sig_adv_per_bilat = occular_occ.loc[lat, 'sig_adv']
    print(f"Ocular adv per {lat} injection {rel_adv_per_bilat:.3f} +- {sig_adv_per_bilat:.3f}")

import scipy.stats

# Done below: correlation of number of bilateral/unilateral with number of adversarial events
# TODO: ocular correlation with systemic

corr_and_sig = {}
corr_and_sig_tot = {}
alldata_n_adverse = pd.DataFrame()
for lat in ['bilateral', 'unilateral']:  # no *2 needed for bilateral as we're only looking at the rank correlation
    for adv_type in ['ocular', 'systemic']:
        alldata_n_adverse_adv_lat = pd.DataFrame()
        alldata_n_injections = pd.DataFrame()
        for name, data in name_data_iter:
            n_adverse = data[f'n_{adv_type}_ae_{lat}']
            n_injections = data[f'n_{lat}']
            hasinj = n_injections > -1  # no selection
            n_adverse = n_adverse[hasinj]
            n_injections = n_injections[hasinj]
            alldata_n_injections = pd.concat([alldata_n_injections, n_injections])
            alldata_n_adverse_adv_lat = pd.concat([alldata_n_adverse_adv_lat, n_adverse])
            spearmancorr = scipy.stats.spearmanr(n_adverse, n_injections)


            corr_and_sig[(name, adv_type, lat)] = spearmancorr
        spearmancorrtot = scipy.stats.spearmanr(alldata_n_adverse_adv_lat, alldata_n_injections)
        corr_and_sig_tot[(adv_type, lat)] = spearmancorrtot
        alldata_n_adverse[f'n_{adv_type}_ae_{lat}'] = alldata_n_adverse_adv_lat[f'n_{adv_type}_ae_{lat}']

# TODO: ocular correlation with systemic
for lat in ['bilateral', 'unilateral']:  # TODO
    spearmanncorr = scipy.stats.spearmanr(alldata_n_adverse[f'n_ocular_ae_{lat}'], alldata_n_adverse[f'n_systemic_ae_{lat}'])
    print(f"Total ocular/systemic {lat} (spearman) correlation: {spearmanncorr.statistic:.2f} p-value: {spearmanncorr.pvalue:.2g}")

for name in diagnoses:
    for adv_type in ['ocular', 'systemic']:
        for lat in ['bilateral', 'unilateral']:
            spearmancorr = corr_and_sig[(name, adv_type, lat)]
            print(f"{name} {adv_type} {lat} (spearman) correlation: {spearmancorr.statistic:.2f} "
                  f"p-value: {spearmancorr.pvalue:.2g}")

for adv_type in ['ocular', 'systemic']:
    for lat in ['bilateral', 'unilateral']:
        spearmancorrtot = corr_and_sig_tot[(adv_type, lat)]
        print(f"Total {adv_type} {lat} (spearman) correlation: {spearmancorrtot.statistic:.2f} "
              f"p-value: {spearmancorrtot.pvalue:.2g}")

print("End")
# TODO: more columns in table, christina will provide more

# Done: how many digits? round df_*_percent to 3 digits

# TODO: number and percent, christina will provide

# Done: hatch? plots in grey, no colors
