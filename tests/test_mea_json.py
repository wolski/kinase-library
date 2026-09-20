import json

import gseapy as gp
import pandas as pd
import pytest

from kinase_library.enrichment.mea import MeaEnrichmentResults


def make_mea_result():
    ranks = pd.Series(
        [8, 7, 7, 5, 4, 3, 2, 1, -1, -2, -3, -4, -5, -6, -7, -8],
        index=[f"g{index}" for index in range(1, 17)],
        name="statistic",
    )
    gene_sets = {
        "positive": ["g1", "g2", "g4", "g6"],
        "negative": ["g11", "g13", "g15", "g16"],
    }
    fitted = gp.prerank(
        rnk=ranks,
        gene_sets=gene_sets,
        weight=1.5,
        min_size=2,
        max_size=10,
        permutation_num=100,
        seed=42,
        verbose=False,
    )
    table = fitted.res2d.rename(
        columns={
            "Term": "Kinase",
            "NOM p-val": "p-value",
            "FDR q-val": "FDR",
            "Tag %": "Subs fraction",
            "Lead_genes": "Leading substrates",
        }
    ).set_index("Kinase")
    return MeaEnrichmentResults(
        enrichment_results=table,
        pps_data=None,
        kin_sub_sets=gene_sets,
        gseapy_obj=fitted,
        kin_type="custom",
        kl_method="custom",
        kl_thresh=None,
        tested_kins=list(gene_sets),
        data_att="custom",
        kl_comp_direction=None,
    )


def test_mea_uses_the_shared_native_gsea_structure(tmp_path):
    result = make_mea_result()
    document = result.to_gsea_result_data("A_vs_B")
    contrast = document["data"]["A_vs_B"]
    category = contrast["categories"]["MEA"]
    native = category["gsea_result"]

    assert list(document["rank_lists"]["A_vs_B"]["entries"]) == [
        f"g{index}" for index in range(1, 17)
    ]
    assert list(contrast["gene_pool"]) == [f"G{index}" for index in range(1, 17)]
    assert native["result"]["columns"]["ID"] == ["negative", "positive"]
    assert set(native["gene_sets"]) == {"negative", "positive"}
    assert native["params"]["exponent"] == 1.5

    for term, source in result.gseapy_obj.results.items():
        assert native["running_scores"][term] == pytest.approx(source["RES"])
        assert native["hit_indices"][term] == [index + 1 for index in source["hits"]]
    assert {term["method"] for term in category["terms"]} == {"gseapy"}

    path = tmp_path / "mea.json"
    result.write_gsea_result_json(path, "A_vs_B")
    assert json.loads(path.read_text()) == document
