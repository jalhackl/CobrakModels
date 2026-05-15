import matplotlib.pyplot as plt
import pandas as pd


def plot_enzyme_analysis(
    df,
    stack_order=None,
    columns_to_plot=None,
    x_col_ex="EX_B",
    x_col_growth="growth",
    supertitle=None
):


    df = df.copy()

    def make_title(title):
        if supertitle:
            return f"{supertitle}\n{title}"
        return title

    # defaults
    if stack_order is None:
        stack_order = ["Transporters", "Respiration", "Glycolysis", "Biomass", "Other"]

    if columns_to_plot is None:
        columns_to_plot = [
            'Glycolysis', 'Respiration',
            'Transport_M', 'Transport_C',
            'Transport_ATP', 'Transport_ADP',
            'EX_S', 'EX_C', 'EX_P', 'EX_B',
            'Biomass_Rea'
        ]


    group_df = pd.DataFrame()

    if "group_totals" in df.columns:
        try:
            group_df = df["group_totals"].apply(pd.Series).astype(float)
        except Exception:
            group_df = pd.DataFrame(index=df.index)

    if x_col_ex in df.columns:
        df_plot = df.sort_values(x_col_ex).copy()
    else:
        df_plot = df.copy()

    if not group_df.empty:
        group_df = group_df.reindex(df_plot.index)

    x_ex = df_plot[x_col_ex].to_numpy() if x_col_ex in df_plot.columns else None

    # stacked plot of enzyme mass vs growth
    if x_ex is not None and not group_df.empty:

        plt.figure(figsize=(8, 5))

        y_stack = []
        for col in stack_order:
            if col not in group_df.columns:
                group_df[col] = 0.0
            y_stack.append(group_df[col].to_numpy())

        try:
            plt.stackplot(x_ex, *y_stack, labels=stack_order, alpha=0.8)
        except Exception:
            pass

        if "total_enzyme_weight" in df_plot.columns:
            plt.plot(x_ex, df_plot["total_enzyme_weight"], "k-", lw=2, label="Total enzyme")
            plt.scatter(x_ex, df_plot["total_enzyme_weight"], c="black", s=10)

        if "prot_pool_delivery" in df_plot.columns:
            plt.plot(x_ex, df_plot["prot_pool_delivery"], "k")
            plt.scatter(x_ex, df_plot["prot_pool_delivery"], c="black", s=10)

        plt.xlabel(x_col_ex)
        plt.ylabel("Enzyme mass")
        plt.title(make_title("Enzyme usage across reaction groups"))
        plt.legend()
        plt.tight_layout()
        plt.show()

    # plot of enzyme mass vs growth
    if x_col_growth in df_plot.columns and not group_df.empty:

        plt.figure()

        for col in stack_order:
            if col in group_df.columns:
                plt.plot(df_plot[x_col_growth], group_df[col], marker="o", label=col)

        plt.xlabel(x_col_growth)
        plt.ylabel("Enzyme mass")
        plt.title(make_title("Enzyme groups vs growth"))
        plt.legend()
        plt.tight_layout()
        plt.show()

    # plot of fluxes
    if x_col_growth in df.columns:

        def style(col):
            if col.startswith("Transport"):
                return {"linestyle": "--", "linewidth": 2}
            if col.startswith("EX"):
                return {"linestyle": ":", "linewidth": 2}
            return {"linestyle": "-", "linewidth": 1.5}

        plt.figure()

        for col in columns_to_plot:
            if col in df.columns:
                plt.plot(
                    df[x_col_growth],
                    df[col],
                    marker="o",
                    label=col,
                    **style(col)
                )

        plt.xlabel(x_col_growth)
        plt.ylabel("Value")
        plt.title(make_title("Flux / transport / biomass variables"))
        plt.legend(ncol=2)
        plt.tight_layout()
        plt.show()

    # plot of concentration ratios and concentrations vs growth


    if x_col_ex not in df.columns:
        return

    def plot_concentrations(base_name, ylabel, title, color1="green", color2="red"):

        plt.figure()

        col_c = f"{base_name}, c"
        col_m = f"{base_name}, m"
        col_single = base_name

        plotted = False

        # if compartment-specific data exists
        if col_c in df.columns and col_m in df.columns:

            plt.plot(df[x_col_ex], df[col_c], c=color1, label=f"{base_name} (c)")
            plt.scatter(df[x_col_ex], df[col_c], c=color1)

            plt.plot(df[x_col_ex], df[col_m], c=color2, label=f"{base_name} (m)")
            plt.scatter(df[x_col_ex], df[col_m], c=color2)

            plotted = True

        # if no compartment
        elif col_single in df.columns:
            plt.plot(df[x_col_ex], df[col_single], c=color1, label=col_single)
            plt.scatter(df[x_col_ex], df[col_single], c=color1)

            plotted = True

        if not plotted:
            plt.close()
            return

        plt.xlabel(x_col_ex)
        plt.ylabel(ylabel)
        plt.title(make_title(title))
        plt.legend()
        plt.tight_layout()
        plt.show()

    # ATP / ADP concentrations
    plot_concentrations(
        "ATP conc.",
        "concentration",
        "ATP concentration"
    )

    plot_concentrations(
        "ADP conc.",
        "concentration",
        "ADP concentration"
    )

    # ATP/ADP ratio (in case of no compartmentalization, or if ratio column is directly available)
    if "ATP/ADP" in df.columns:
        plt.figure()
        plt.plot(df[x_col_ex], df["ATP/ADP"].round(6), c="purple", label="ATP/ADP")
        plt.scatter(df[x_col_ex], df["ATP/ADP"].round(6), c="purple")

        plt.xlabel(x_col_ex)
        plt.ylabel("ratio")
        plt.title(make_title("ATP/ADP ratio"))
        plt.legend()
        plt.tight_layout()
        plt.show()