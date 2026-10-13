#!/usr/bin/env python3
"""
Reduce a Hipparcos-like FITS catalog to the columns used by the Julia project:
RAdeg, DEdeg, Vmag

It filters rows with:
- Vmag <= mag_lim
- finite Vmag

Output FITS layout:
- HDU 0: PrimaryHDU
- HDU 1: BinTableHDU with columns RAdeg, DEdeg, Vmag

This keeps compatibility with Julia code that reads f[2].
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
from astropy.io import fits


REQUIRED_COLUMNS = ("RAdeg", "DEdeg", "Vmag")
PREFERRED_EXTENSIONS = ("hip_main.dat", "HIP_MAIN.DAT")


def _normalize_names(names):
    """Map normalized column names to original names."""
    if not names:
        return {}
    return {str(name).strip().lower(): str(name) for name in names}


def _has_required_columns(hdu) -> bool:
    """Return True if the HDU contains the required columns."""
    columns = getattr(hdu, "columns", None)
    names = getattr(columns, "names", None)
    normalized = _normalize_names(names)
    return all(col.lower() in normalized for col in REQUIRED_COLUMNS)


def _resolve_column_names(hdu):
    """Return the exact column names present in the HDU."""
    columns = getattr(hdu, "columns", None)
    names = getattr(columns, "names", None)
    normalized = _normalize_names(names)
    try:
        return tuple(normalized[col.lower()] for col in REQUIRED_COLUMNS)
    except KeyError as exc:
        raise ValueError(
            f"La extensión seleccionada no contiene las columnas requeridas {REQUIRED_COLUMNS}."
        ) from exc


def find_table_hdu(hdul: fits.HDUList, preferred_ext: str | None = "hip_main.dat"):
    """
    Return the most appropriate table HDU containing RAdeg, DEdeg and Vmag.

    The input catalog may store the table as TableHDU instead of BinTableHDU,
    so both types are accepted.
    """
    if preferred_ext:
        for hdu in hdul:
            if str(getattr(hdu, "name", "")).strip().lower() == preferred_ext.strip().lower():
                if _has_required_columns(hdu):
                    return hdu
                raise ValueError(
                    f"La extensión '{preferred_ext}' existe, pero no contiene {REQUIRED_COLUMNS}."
                )

    for ext_name in PREFERRED_EXTENSIONS:
        for hdu in hdul:
            if str(getattr(hdu, "name", "")).strip().lower() == ext_name.lower() and _has_required_columns(hdu):
                return hdu

    for hdu in hdul:
        if isinstance(hdu, (fits.TableHDU, fits.BinTableHDU)) and _has_required_columns(hdu):
            return hdu

    raise ValueError(
        "No se encontró ninguna extensión de tabla con las columnas "
        f"{REQUIRED_COLUMNS}."
    )


def reduce_catalog(
    input_path: str | Path,
    output_path: str | Path,
    mag_lim: float = 10.0,
    preferred_ext: str | None = "hip_main.dat",
) -> None:
    """Read, filter and write the reduced catalog."""
    input_path = Path(input_path)
    output_path = Path(output_path)

    with fits.open(input_path, memmap=True) as hdul:
        table_hdu = find_table_hdu(hdul, preferred_ext=preferred_ext)
        ra_name, dec_name, vmag_name = _resolve_column_names(table_hdu)
        data = table_hdu.data

        input_columns = list(table_hdu.columns.names)
        input_num_columns = len(input_columns)
        input_num_stars = len(data)

        ra = np.asarray(data[ra_name], dtype=np.float64)
        dec = np.asarray(data[dec_name], dtype=np.float64)
        vmag = np.asarray(data[vmag_name], dtype=np.float64)

        mask = (vmag <= mag_lim) & np.isfinite(vmag)

        ra_out = ra[mask]
        dec_out = dec[mask]
        vmag_out = vmag[mask]

        output_columns = list(REQUIRED_COLUMNS)
        output_num_columns = len(output_columns)
        output_num_stars = len(vmag_out)

        cols = fits.ColDefs([
            fits.Column(name="RAdeg", format="D", array=ra_out),
            fits.Column(name="DEdeg", format="D", array=dec_out),
            fits.Column(name="Vmag", format="D", array=vmag_out),
        ])

        primary_hdu = fits.PrimaryHDU()
        table_out = fits.BinTableHDU.from_columns(cols, name="HIPPARCOS_FILTERED")

        table_out.header["MAGLIM"] = (float(mag_lim), "Maximum Vmag kept")
        table_out.header["NINPUT"] = (int(input_num_stars), "Rows/stars in input table")
        table_out.header["NOUTPUT"] = (int(output_num_stars), "Rows/stars after filtering")
        table_out.header["NCOLIN"] = (int(input_num_columns), "Columns in input table")
        table_out.header["NCOLOUT"] = (int(output_num_columns), "Columns in output table")
        table_out.header["NREJECT"] = (int(input_num_stars - output_num_stars), "Rows/stars rejected")
        table_out.header["SRCFILE"] = (input_path.name[:68], "Original FITS file")
        table_out.header["SRCEXT"] = (str(getattr(table_hdu, "name", ""))[:68], "Source extension")

        fits.HDUList([primary_hdu, table_out]).writeto(output_path, overwrite=True)

    print(f"Catálogo reducido guardado en: {output_path}")
    print(f"Extensión usada: {getattr(table_hdu, 'name', '<sin nombre>')}")
    print("")
    print("Resumen de reducción:")
    print(f"- Columnas de entrada: {input_num_columns}")
    print(f"- Columnas de salida:  {output_num_columns}")
    print(f"- Estrellas de entrada: {input_num_stars}")
    print(f"- Estrellas de salida:  {output_num_stars}")
    print(f"- Estrellas descartadas: {input_num_stars - output_num_stars}")
    print(f"- Límite de magnitud aplicado: Vmag <= {mag_lim}")


def build_parser() -> argparse.ArgumentParser:
    """Create the CLI parser."""
    parser = argparse.ArgumentParser(
        description="Reduce un catálogo FITS a RAdeg, DEdeg y Vmag, filtrando por Vmag."
    )
    parser.add_argument("input_fits", help="Ruta al FITS original")
    parser.add_argument("output_fits", help="Ruta del FITS reducido de salida")
    parser.add_argument(
        "--mag-lim",
        type=float,
        default=10.0,
        help="Límite de magnitud V a conservar (por defecto: 10.0)",
    )
    parser.add_argument(
        "--ext",
        default="hip_main.dat",
        help="Nombre de la extensión a leer (por defecto: hip_main.dat)",
    )
    return parser


def main() -> None:
    """CLI entry point."""
    parser = build_parser()
    args = parser.parse_args()
    reduce_catalog(
        args.input_fits,
        args.output_fits,
        mag_lim=args.mag_lim,
        preferred_ext=args.ext,
    )


if __name__ == "__main__":
    main()
