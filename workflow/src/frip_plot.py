import pandas as pd
import plotly.graph_objects as go


def main():
    inputs = list(snakemake.input)
    output_html = str(snakemake.output[0])

    pd.options.plotting.backend = "plotly"
    dataframes = []
    for path in sorted(inputs):
        parts = path.split("_")
        # expecting .../frip_{sample}.tsv -> use last 2 underscore parts of filename without extension
        # fallback to whole stem if unexpected
        try:
            cond_marker = "_".join(parts[-2:]).replace(".tsv", "")
        except Exception:
            cond_marker = path.rsplit("/", 1)[-1].replace(".tsv", "")
        temp_df = (
            pd.read_csv(path, sep="\t", usecols=["percent"]).rename(columns={"percent": cond_marker})
        )
        dataframes.append(temp_df)

    frip_df = pd.concat(dataframes, axis=1)
    frip_df = frip_df.rename(index={0: "inside"})
    frip_df.loc["outside"] = 100 - frip_df.loc["inside"]

    fig = go.Figure(
        data=[
            go.Bar(
                name="inside_peaks",
                x=frip_df.columns,
                y=frip_df.loc["inside"],
                marker_color="rgb(255, 201, 57)",
            ),
            go.Bar(
                name="outside_peaks",
                x=frip_df.columns,
                y=frip_df.loc["outside"],
                marker_color="rgb(0,39, 118)",
            ),
        ]
    )
    fig.update_layout(
        barmode="stack",
        title="Fraction of Reads in Peaks by Sample",
        xaxis=dict(title="Samples", tickfont=dict(size=14)),
        yaxis=dict(title="Fraction of reads in peaks", tickfont=dict(size=14)),
    )
    fig.write_html(output_html)


if __name__ == "__main__":
    main()


