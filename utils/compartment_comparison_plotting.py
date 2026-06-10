
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def compute_compartment_advantage(df):
    comp = df[df["compartment_scenario"] == "comp"].copy()
    nocomp = df[df["compartment_scenario"] == "nocomp"].copy()

    comp = comp.rename(columns={"ATP_Consumption": "ATP_comp"})
    nocomp = nocomp.rename(columns={"ATP_Consumption": "ATP_nocomp"})

    merged = pd.merge(
        comp,
        nocomp,
        on=["growth", "transport_dG0", "atp_scenario"],
        suffixes=("_comp", "_nocomp"),
    )

    merged = merged.dropna(subset=["ATP_comp", "ATP_nocomp"])

    merged["advantage"] = merged["ATP_comp"] - merged["ATP_nocomp"]

    return merged


def plot_protein_cost_vs_growth(df, title="Protein Cost vs Growth Rate"):
    """Plot protein pool delivery vs growth rate, with separate lines for each dG0 value."""
    plt.figure(figsize=(10, 6))

    for scenario in sorted(df['scenario'].unique()):
        subset = df[df['scenario'] == scenario].sort_values('growth')
        for dg0_val in sorted(subset['transport_dG0'].unique()):
            dg0_subset = subset[subset['transport_dG0'] == dg0_val]
            label = f"{scenario.replace('_', ' ')} dG0={dg0_val}"
            plt.plot(
                dg0_subset['growth'],
                dg0_subset['ATP_Consumption'],
                marker='o',
                markersize=4,
                label=label,
                alpha=0.7,
            )

    plt.xlabel('Growth rate')
    plt.ylabel('Protein pool delivery (g/gDW)')
    plt.title(title)
    plt.legend(fontsize=7, ncol=3, loc='upper left', bbox_to_anchor=(1.02, 1))
    plt.tight_layout()
    plt.show()


def plot_cost_vs_growth_by_dg0(df, title="ATP Consumption vs Growth Rate"):
    """One subplot per transport_dG0, comparing scenarios within each."""

    dg0_values = sorted(df['transport_dG0'].unique())
    scenarios = sorted(df['scenario'].unique())

    fig, axes = plt.subplots(
        len(dg0_values),
        1,
        figsize=(10, 4 * len(dg0_values)),
        sharex=True,
        sharey=True
    )

    if len(dg0_values) == 1:
        axes = [axes]  

    for ax, dg0_val in zip(axes, dg0_values):
        subset_dg0 = df[df['transport_dG0'] == dg0_val]

        for scenario in scenarios:
            scenario_data = subset_dg0[subset_dg0['scenario'] == scenario] \
                .sort_values('growth')

            ax.plot(
                scenario_data['growth'],
                scenario_data['ATP_Consumption'],
                marker='o',
                markersize=3,
                label=scenario
            )

        ax.set_title(f"dG0 = {dg0_val}")
        ax.set_ylabel("ATP Consumption")
        ax.grid(alpha=0.3)

    axes[-1].set_xlabel("Growth rate")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", fontsize=8)

    plt.suptitle(title, y=1.02)
    plt.tight_layout()
    plt.show()


def plot_advantage_heatmap(adv_df, title="Compartment Advantage Heatmap"):
    """Plot compartment advantage as a heatmap with growth rate on x-axis and dG0 on y-axis."""
    for scenario in sorted(adv_df['atp_scenario'].unique()):
        subset = adv_df[adv_df['atp_scenario'] == scenario].sort_values('growth')

        pivot = subset.pivot_table(
            values='advantage',
            index='transport_dG0',
            columns='growth',
            aggfunc='mean',
        )

        pivot = pivot.reindex(index=sorted(pivot.index))

        plt.figure(figsize=(10, 6))
        sns.heatmap(
            pivot,
            annot=True,
            fmt='.4f',
            cmap='RdBu_r',
            center=0,
            linewidths=0.5,
            cbar_kws={'label': 'Advantage', 'shrink': 0.8},
           # vmin=-0.5,
           # vmax=0.5,
        )
        plt.xlabel('Growth rate')
        plt.ylabel('Transporter dG0')
        plt.title(f"{title}\nScenario: {scenario}")
        plt.tight_layout()
        plt.show()


def plot_dg0_effect_by_growth(adv_df):
    """For each growth rate, show how advantage changes with dG0."""
    growth_values = sorted(adv_df['growth'].unique())
    n_growth = len(growth_values)

    for scenario in ['comp_noATP', 'comp_ATP']:
        subset = adv_df[adv_df['scenario_comp'] == scenario].sort_values('growth')

        fig, axes = plt.subplots(1, n_growth, figsize=(5 * n_growth, 4))
        if n_growth == 1:
            axes = [axes]

        for gi, growth_val in enumerate(growth_values):
            g_subset = subset[subset['growth'] == growth_val].sort_values('transport_dG0')
            bars = axes[gi].bar(
                [str(int(v)) for v in g_subset['transport_dG0']],
                g_subset['advantage'].values,
                color=['green' if v > 0 else 'red' for v in g_subset['advantage'].values],
                edgecolor='white',
            )
            for bar, val in zip(bars, g_subset['advantage'].values):
                axes[gi].text(
                    bar.get_x() + bar.get_width()/2.,
                    bar.get_height() + 0.001,
                    f'{val:+.3f}',
                    ha='center', va='bottom', fontsize=8,
                )
            axes[gi].set_ylabel('Advantage')
            axes[gi].set_title(f'growth = {growth_val}')
            axes[gi].axhline(y=0, color='gray', linestyle=':', alpha=0.5)

        fig.suptitle(f'Advantage vs dG0 at Each Growth Rate\nScenario: {scenario.replace("_", " ")}', fontsize=13)
        plt.tight_layout()
        plt.show()



def plot_enzyme_heatmaps(df, variables):
    df = df.copy()

    for var in variables:

        pivot = df.pivot_table(
            values=var,
            index="transport_dG0",
            columns="growth",
            aggfunc="mean"
        )

        pivot = pivot.sort_index().sort_index(axis=1)

        plt.figure(figsize=(6, 4))

        sns.heatmap(
            pivot,
            cmap="viridis",
            linewidths=0.5,
            cbar_kws={"label": var},
            mask=pivot.isna()  
        )

        plt.title(var)
        plt.xlabel("Growth rate")
        plt.ylabel("Transport dG0")
        plt.tight_layout()
        plt.show()


def plot_enzyme_lines(df, variables):
    df = df.copy()

    dG0_values = sorted(df["transport_dG0"].unique())

    for var in variables:

        plt.figure(figsize=(7, 5))

        for dG0 in dG0_values:

            sub = df[df["transport_dG0"] == dG0]

            grouped = sub.groupby("growth")[var].mean().sort_index()

            plt.plot(
                grouped.index,
                grouped.values,
                marker="o",
                label=f"dG0 = {dG0}"
            )

        plt.title(var)
        plt.xlabel("Growth rate")
        plt.ylabel(var)
        plt.legend(title="Transport dG0")
        plt.tight_layout()
        plt.show()


def plot_atp_adp_ratio_vs_growth(df, title="ATP/ADP Ratio vs Growth Rate"):
    """Plot ATP/ADP ratio vs growth rate, with separate lines for each dG0 value."""
    # Compute ATP/ADP ratio
    df_plot = df.copy()
    atp = np.exp(df_plot['x_ATP'].fillna(0))
    adp = np.exp(df_plot['x_ADP'].fillna(0))
    df_plot['ATP_ADP_ratio'] = atp / (adp + 1e-15)

    for scenario in sorted(df_plot['scenario'].unique()):
        subset = df_plot[df_plot['scenario'] == scenario].sort_values('growth')

        plt.figure(figsize=(10, 6))
        for dg0_val in sorted(subset['transport_dG0'].unique()):
            dg0_subset = subset[subset['transport_dG0'] == dg0_val]
            label = f"{scenario.replace('_', ' ')} dG0={dg0_val}"
            plt.plot(
                dg0_subset['growth'],
                dg0_subset['ATP_ADP_ratio'],
                marker='o',
                markersize=4,
                label=label,
                alpha=0.7,
            )

        plt.xlabel('Growth rate')
        plt.ylabel('ATP/ADP ratio')
        plt.title(f"{title}\nScenario: {scenario.replace('_', ' ')}")
        plt.legend(fontsize=7, ncol=3, loc='upper left', bbox_to_anchor=(1.02, 1))
        plt.tight_layout()
        plt.show()