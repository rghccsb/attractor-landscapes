import sqlite3
from os.path import join

from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator
import pandas as pd
import seaborn as sns


def plot_performance(df):
    fig, axs = plt.subplots(2, 1)
    ax_line = sns.lineplot(
        data=df,
        x='Period, $n$',
        y='Time (s)',
        err_style='bars',
        ci='sd',
        ax=axs[0],
    )
    axs[0].set_title('$\\bf{a.}$ Mean Time Taken to Search Attractors\nper Period for Model 18', fontsize=14)
    axs[0].xaxis.label.set_size(12)
    axs[0].yaxis.label.set_size(12)
    ax_line.set_yscale('log')

    ax_bar = sns.barplot(
        data=df,
        x='Period, $n$',
        y='Number of Attractors',
        ax=axs[1],
    )
    axs[1].set_title('$\\bf{b.}$ Mean Number of Attractors per Period for Model 18', fontsize=14)
    axs[1].set_ylabel('Number of\nAttractors')
    axs[1].xaxis.label.set_size(12)
    axs[1].yaxis.label.set_size(12)
    fig.tight_layout()
    plt.savefig(join('.', 'figures', 'performance.pdf'))


if __name__ == '__main__':
    conn = sqlite3.connect(join('.', 'models', '20200415_cv_MCM_Guimera.db'))
    df = pd.read_sql_query(
        """
SELECT * FROM steady_state_performance;
        """,
        conn,
    ).rename(columns={
        'period': 'Period, $n$',
        'time_delta': 'Time (s)',
        'num_attractors': 'Number of Attractors',
    })
    plot_performance(df)
