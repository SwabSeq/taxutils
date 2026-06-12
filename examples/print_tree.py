import pandas as pd

from taxutils import taxutils


tu = taxutils()

taxon_counts = pd.Series(
    {
        12059: 500,
        147711: 120,
        147712: 80,
    },
    name="count",
)

tree = tu.format_tree(taxon_counts.index)

report = pd.DataFrame(
    {
        "count": taxon_counts.reindex(tree.index).fillna(0).astype(int),
        "name": tree,
    }
)

with open("examples/example_kreport.txt", "w") as f:
    for taxon, row in report.iterrows():
        f.write(f"{taxon}\t{row['count']}\t{row['name']}\n")

print(report.to_string())
