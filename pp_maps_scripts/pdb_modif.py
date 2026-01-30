from biopandas.pdb import PandasPdb


def modify_pdb(biopdb: PandasPdb) -> PandasPdb:

    for atom in range(
        len(biopdb.df["ATOM"]) - 1, -1, -1
    ):  # Lit le df de la fin au début, avec un pas de -1 pour remonter les itérations
        res_nb_init = biopdb.df["ATOM"]["residue_number"][atom]
        biopdb.df["ATOM"].at[atom, "residue_number"] = 1
        biopdb.df["ATOM"].at[atom, "residue_name"] = "LIG"
        if biopdb.df["ATOM"]["residue_number"][atom - 1] > res_nb_init:
            break

    return biopdb
