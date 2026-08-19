"""
Module description:

Convert Dutch Emission Registration (ER) Excel files
to semicolon-separated CSV files.
The Excel worksheet is expected to contain an 'Emissies'
sheet and an 'Emissie' column.

Author: Dr. Arseni Doyennel
Email: a.doyennel@vu.nl
"""

from pathlib import Path
import pandas as pd

from emission_preparation_setting import (
    point_ER_excel,
    map_ER_excel,
    datadir_point_source_plume_processing,
    point_source_harm_file,
    mapped_data,
)

class ERDataConverter:

    def __init__(
        self,
        point_excel,
        map_excel,
        point_csv,
        map_csv,
    ):
        self.point_excel = Path(point_excel)
        self.map_excel = Path(map_excel)

        self.point_csv = Path(point_csv)
        self.map_csv = Path(map_csv)

        self.sheet_name = "Emissies"
        self.emission_column = "Emissie"

    # =========================================================
    # Convert one Excel file
    # =========================================================

    def convert_excel_to_csv(self, input_file, output_file):

        # -----------------------------------------------------
        # IMPORTANT:
        # Check output BEFORE reading Excel.
        # -----------------------------------------------------

        if output_file.exists():
            print(
                f"WARNING: output file already exists, "
                f"conversion is unnecessary:\n"
                f"    {output_file}"
            )
            return None

        print()
        print("=" * 60)
        print("ER DATA CONVERSION")
        print("=" * 60)
        print(f"Input : {input_file}")
        print(f"Output: {output_file}")

        # -----------------------------------------------------
        # Check input
        # -----------------------------------------------------

        if not input_file.exists():
            raise FileNotFoundError(
                f"ER Excel file not found:\n{input_file}"
            )

        # -----------------------------------------------------
        # Read Excel
        # -----------------------------------------------------

        df = pd.read_excel(
            input_file,
            sheet_name=self.sheet_name,
            dtype=str
        )

        print(f"Rows read: {len(df)}")

        # -----------------------------------------------------
        # Check emission column
        # -----------------------------------------------------

        if self.emission_column not in df.columns:
            raise ValueError(
                f"Column '{self.emission_column}' was not found "
                f"in sheet '{self.sheet_name}'.\n\n"
                f"Available columns:\n{list(df.columns)}"
            )

        # -----------------------------------------------------
        # Convert decimal comma -> decimal point
        # -----------------------------------------------------

        df[self.emission_column] = (
            df[self.emission_column]
            .fillna("")
            .str.strip()
            .str.replace(",", ".", regex=False)
        )

        # -----------------------------------------------------
        # Create output directory if necessary
        # -----------------------------------------------------

        output_file.parent.mkdir(
            parents=True,
            exist_ok=True
        )

        # -----------------------------------------------------
        # Write CSV
        # -----------------------------------------------------

        df.to_csv(
            output_file,
            sep=";",
            index=False,
            encoding="utf-8",
            lineterminator="\n"
        )

        print(f"CSV written successfully:")
        print(f"    {output_file}")

        return df

    # =========================================================
    # Convert point source
    # =========================================================

    def convert_point(self):

        return self.convert_excel_to_csv(
            input_file=self.point_excel,
            output_file=self.point_csv
        )

    # =========================================================
    # Convert gridded map
    # =========================================================

    def convert_map(self):

        return self.convert_excel_to_csv(
            input_file=self.map_excel,
            output_file=self.map_csv
        )

    # =========================================================
    # Convert both
    # =========================================================

    def convert_all(self):

        print()
        print("#" * 60)
        print("# Checking ER conversion files")
        print("#" * 60)

        df_point = self.convert_point()
        df_map = self.convert_map()

        print()
        print("#" * 60)
        print("# ER data conversion check finished")
        print("#" * 60)

        return df_point, df_map
        
if __name__ == "__main__":

    converter = ERDataConverter(
        point_excel=(
            datadir_point_source_plume_processing
            + point_ER_excel
        ),
        map_excel=(
            datadir_point_source_plume_processing
            + map_ER_excel
        ),
        point_csv=(datadir_point_source_plume_processing + point_source_harm_file),
        map_csv=(datadir_point_source_plume_processing + mapped_data),
    )

    converter.convert_all()
