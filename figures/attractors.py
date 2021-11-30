import json
from os.path import join
import sqlite3
import numpy as np

from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import patches
from sklearn.manifold import MDS
import pandas as pd

cytokine_storm = {
    13: [2, 0, 0, 0, 1, 1, 0, 0, 2, 2, 0, 0, 2, 2, 0, 1, 1, 1, 2],
    15: [2, 0, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0, 2, 1, 0, 1, 0, 1, 2],
    18: [2, 0, 0, 0, 1, 1, 0, 0, 2, 1, 0, 0, 2, 0, 0, 1, 0, 1, 2],
}

rest_states = {
    13: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    15: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    18: [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
}


def plot_attractors(df):
    fig, axs = plt.subplots(
        len(df['drug'].unique()),
        len(df['model_id'].unique()),
        sharex=True,
        sharey=True,
        figsize=(12, 14),
        squeeze=False,
        gridspec_kw={'wspace': 0.05},
        constrained_layout=True,
    )
    style = 'Simple, tail_width=0.5, head_width=4, head_length=8'
    kw = dict(arrowstyle=style, color='k', alpha=0.4)
    jitter = 0.1
    num_points = 25
    colors = cm.plasma(np.linspace(0, 1, 3))
    for i, model_id in enumerate(reversed(df['model_id'].unique())):
        for j, drug in enumerate(df['drug'].unique()):
            df_sub = df[(df['model_id'] == model_id) & (df['drug'] == drug)]
            periods_attractors_cv_values = zip(
                df_sub['period'],
                df_sub['states_2d'],
                df_sub['cv_value'],
            )
            for period, attractor, cv_values in periods_attractors_cv_values:
                ax = axs[j, i]
                attractor = np.array(attractor)
                if period == 'rest':
                    for k in range(num_points):
                        ax.scatter(
                            attractor[:, 0],
                            attractor[:, 1],
                            s=k ** 3.2,
                            alpha=((num_points - k) / num_points) / 10,
                            color='green',
                        )
                    ax.annotate(
                        'Immune\nQuiescence',
                        (attractor[:, 0], attractor[:, 1] + 1.0),
                        ha='center',
                        fontsize=16,
                    )
                if period == 'cytokine_storm':
                    for k in range(num_points):
                        ax.scatter(
                            attractor[:, 0] + jitter,
                            attractor[:, 1] - jitter,
                            s=k ** 3.2,
                            alpha=((num_points - k) / num_points) / 10,
                            color='red',
                        ).set_zorder(1)
                    ax.annotate(
                        'Cytokine\nStorm',
                        (attractor[:, 0], attractor[:, 1] + 0.5),
                        ha='center',
                        fontsize=16,
                    )
                if period not in {'rest', 'cytokine_storm'}:
                    ax.scatter(
                        attractor[:, 0],
                        attractor[:, 1],
                        c=list(map(lambda x: colors[x], cv_values)),
                    ).set_zorder(3)
                    for k, state in enumerate(attractor[:-1]):
                        arrow = patches.FancyArrowPatch(
                            state,
                            attractor[k + 1, :],
                            connectionstyle='arc3,rad=.45',
                            **kw,
                        )
                        ax.add_patch(arrow)
                ax.set_xticks([])
                ax.set_yticks([])
                if model_id == 18:
                    if drug == 'None':
                        ax.set_title(
                            f"""$\\bf{{a.}}$ Model {model_id}: Attractors""",
                            fontsize=16,
                        )
                    elif drug == 'antiviral':
                        ax.set_title(
                            f"""$\\bf{{e.}}$ Model {model_id}: Attractors with Idealized Antiviral""",
                            fontsize=16,
                        )
                    elif drug == 'hydroxychloroquine':
                        ax.set_title(
                            f"""$\\bf{{c.}}$ Model {model_id}: Attractors with Hydroxychloroquine""",
                            fontsize=16,
                        )
                    elif drug == 'hydroxychloroquine_antiviral':
                        ax.set_title(
                            f"""$\\bf{{g.}}$ Model {model_id}: Attractors with
Hydroxychloroquine & Idealized Antiviral""",
                            fontsize=16,
                        )
                elif model_id == 15:
                    if drug == 'None':
                        ax.set_title(
                            f"""$\\bf{{b.}}$ Model {model_id}: Attractors""",
                            fontsize=16,
                        )
                    elif drug == 'antiviral':
                        ax.set_title(
                            f"""$\\bf{{f.}}$ Model {model_id}: Attractors with Idealized Antiviral""",
                            fontsize=16,
                        )
                    elif drug == 'hydroxychloroquine':
                        ax.set_title(
                            f"""$\\bf{{d.}}$ Model {model_id}: Attractors with Hydroxychloroquine""",
                            fontsize=16,
                        )
                    elif drug == 'hydroxychloroquine_antiviral':
                        ax.set_title(
                            f"""$\\bf{{h.}}$ Model {model_id}: Attractors with
Hydroxychloroquine & Idealized Antiviral""",
                            fontsize=16,
                        )
    fake_points = [
        plt.Line2D([0], [1], color=c, marker='o', linestyle='')
        for c in colors
    ]
    fig.legend(
        fake_points,
        [
            'No coronavirus\n(viral titer $\\leq 10^3$ pfu/ml)',
            'Low coronavirus\n($10^3 <$ viral titer $\\leq 10^7$ pfu/ml)',
            'High coronavirus\n(viral titer $> 10^7$ pfu/ml)',
        ],
        loc='lower center',
        fontsize=16,
    )
    # fig.tight_layout()
    plt.savefig(join('.', 'figures', 'attractors.pdf'))


def reduce_attractor_dimensions(df):
    model_dfs = []
    for model_id in df['model_id'].unique():
        model_df = df[df['model_id'] == model_id]
        attractor_states = [
            i for s in model_df['states'] for i in s
        ] + [cytokine_storm[model_id]] + [rest_states[model_id]]
        mds = MDS(random_state=850732)
        reduced_attractor_states = mds.fit_transform(attractor_states)
        mapping = dict(zip(
            [tuple(s) for s in attractor_states],
            reduced_attractor_states,
        ))
        all_states_2d = []
        for states in model_df['states']:
            states_2d = []
            for state in states:
                states_2d += [mapping[tuple(state)]]
            all_states_2d += [states_2d]
        model_df = model_df.assign(states_2d=all_states_2d)
        for drug in df['drug'].unique():
            model_df = model_df.append({
                'model_id': model_id,
                'period': 'rest',
                'drug': drug,
                'states': [rest_states[model_id]],
                'states_2d': [mapping[tuple(rest_states[model_id])]],
            }, ignore_index=True)
            model_df = model_df.append({
                'model_id': model_id,
                'period': 'cytokine_storm',
                'drug': drug,
                'states': [cytokine_storm[model_id]],
                'states_2d': [mapping[tuple(cytokine_storm[model_id])]],
            }, ignore_index=True)
        model_dfs += [model_df]
    return pd.concat(model_dfs)


def parse_steady_state_col(df):
    steady_state_dicts = list(map(json.loads, df['steady_state']))
    steady_state = pd.json_normalize(steady_state_dicts)
    return pd.merge(
        df,
        steady_state,
        left_index=True,
        right_index=True,
    ).drop('steady_state', axis=1)


if __name__ == '__main__':
    conn = sqlite3.connect(join('.', 'models', '20200415_cv_MCM_Guimera.db'))
    df = pd.read_sql_query(
        """
SELECT steady_state, model_id FROM steady_state
WHERE (model_id = 18 OR model_id = 15)
GROUP BY model_id, JSON_EXTRACT(steady_state, "$.states"),
JSON_EXTRACT(steady_state, "$.drug");
        """,
        conn,
    )
    df = parse_steady_state_col(df)
    df = reduce_attractor_dimensions(df)
    df.loc[df['drug'].isna(), 'drug'] = 'None'
    df['cv_value'] = df.apply(
        lambda row: [state[-1] for _, state in enumerate(row['states'])],
        axis=1,
    )
    plot_attractors(df)
