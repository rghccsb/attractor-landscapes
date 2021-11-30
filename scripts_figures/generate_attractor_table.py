from collections import defaultdict
import json
from os.path import join
import sqlite3

import numpy as np


def generate_attractor_table(cursor, table_path):
    with open(table_path, 'w') as fh:
        fh.write('| {0} |\n'.format(
            ' | '.join([
                'Model ID',
                '# of $n = 1$ Attractors',
                '# of $n = 2$ Attractors',
                '# of $n = 3$ Attractors',
                '# of $n = 4$ Attractors',
                'Total # of Attractors',
                'Minimum Rest Distance',
                'Minimum Cytokine Storm Distance',
            ]),
        ))
        fh.write(
            '|----------+-----------------------+-------------------------+-------------------------+-------------------------|\n',
        )
        model_attractors = defaultdict(dict)
        model_id_period_pairs = cursor.execute(
            """
SELECT DISTINCT model_id, JSON_EXTRACT(steady_state, "$.period")
FROM steady_state;
            """,
        ).fetchall()
        for model_id_period in model_id_period_pairs:
            attractors = cursor.execute(
                """
SELECT steady_state FROM steady_state WHERE
model_id = ? AND
JSON_EXTRACT(steady_state, "$.period") = ? AND
JSON_EXTRACT(steady_state, "$.drug") IS NULL;
                """,
                model_id_period,
            ).fetchall()
            model_attractors[model_id_period[0]].update({
                model_id_period[1]: {
                    'num': len(attractors),
                },
            })
        model_ids = cursor.execute(
            """
SELECT DISTINCT model_id FROM steady_state;
            """,
        ).fetchall()
        for model_id in model_ids:
            attractors = cursor.execute(
                """
SELECT steady_state FROM steady_state WHERE model_id = ? AND
JSON_EXTRACT(steady_state, "$.drug") IS NULL;
                """,
                model_id,
            ).fetchall()
            cytokine_storm = np.array(json.loads(
                cursor.execute(
                    """
SELECT JSON_EXTRACT(model, "$.trajectory_outputs[0][10]") FROM model WHERE
rowid = ?;
                    """,
                    (model_id[0],),
                ).fetchone()[0],
            ))

            attractor_states = [
                np.array(json.loads(a[0])['states']) for a in attractors
            ]
            cytokine_storm_dists = [
                np.sum(np.abs(
                    a[:-1, :] - cytokine_storm,
                ))
                for a in attractor_states if a.shape[0] == 2
            ]
            rest_dists = [
                np.sum(
                    a[:-1, :],
                )
                for a in attractor_states if a.shape[0] == 2
            ]
            model_attractors[model_id[0]].update({
                'min_cytokine_storm_dist': min(cytokine_storm_dists),
                'min_rest_dist': min(rest_dists),
            })

        for model_id, attractors in model_attractors.items():
            if 1 in attractors:
                num1attractors = attractors.get(1).get('num')
            else:
                num1attractors = 0
            if 2 in attractors:
                num2attractors = attractors.get(2).get('num')
            else:
                num2attractors = 0
            if 3 in attractors:
                num3attractors = attractors.get(3).get('num')
            else:
                num3attractors = 0
            if 4 in attractors:
                num4attractors = attractors.get(4).get('num')
            else:
                num4attractors = 0
            num_total_attractors = num1attractors + num2attractors \
                + num3attractors + num4attractors
            fh.write('| {0} |\n'.format(
                ' | '.join(
                    map(
                        str,
                        [
                            model_id,
                            num1attractors,
                            num2attractors,
                            num3attractors,
                            num4attractors,
                            num_total_attractors,
                            attractors['min_rest_dist'],
                            attractors['min_cytokine_storm_dist'],
                        ],
                    ),
                ),
            ))


if __name__ == '__main__':
    db = sqlite3.connect(join('models', '20200415_cv_MCM_Guimera.db'))
    cursor = db.cursor()
    generate_attractor_table(cursor, join('results', 'attractor_table.md'))
