import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

from inputs import RESULTS, WORK, unpack_topic
from counties import read_rows, run as fit_counties

SCRIPTS = Path(__file__).resolve().parent


def check_stations():
    import json
    import numpy as np

    reference = read_rows(RESULTS / "stations.csv")
    for model in json.loads((RESULTS / "models.json").read_text(encoding="utf-8")):
        variable = model["variable"]
        cv = read_rows(WORK / "tables" / (model["id"] + "_cv.csv"))
        by_id = {int(row["station_id"]): row for row in cv}
        assert len(cv) == len(by_id) == len(reference) == 53
        prediction = [float(by_id[int(r["id"])]["predicted"]) for r in reference]
        observed = np.array([float(r[variable]) for r in reference])
        np.testing.assert_allclose(prediction, [float(r[variable + "_loo"]) for r in reference],
                                   atol=1e-4, rtol=1e-7)
        print(f"{variable} RMSE = {np.sqrt(np.mean((prediction - observed) ** 2)):.2f} mm")


def main():
    parser = argparse.ArgumentParser(description="甘肃地统计实习")
    parser.add_argument("--full", action="store_true", help="运行空间验证、变程敏感性和台站回归")
    parser.add_argument("--figures", action="store_true", help="用 R 重绘报告图")
    parser.add_argument("--arcmap", metavar="PYTHON", help="使用 ArcMap 10.8 的 Python 从原始栅格重算县域统计")
    options = parser.parse_args()
    subprocess.run([sys.executable, SCRIPTS / "kriging.py"], check=True)
    check_stations()
    if options.arcmap or options.figures:
        unpack_topic()
    zonal_table = None
    if options.arcmap:
        subprocess.run([options.arcmap, SCRIPTS / "zonal.py"], check=True)
        zonal_table = WORK / "tables" / "county_zones.csv"
    fit_counties(zonal_table)
    if options.full:
        subprocess.run([sys.executable, SCRIPTS / "spatial.py"], check=True)
    if options.figures:
        rscript = shutil.which("Rscript")
        env = os.environ.copy()
        if os.name == "nt":
            env["LC_ALL"] = "Chinese (Simplified)_China.utf8"
            if not rscript:
                versions = sorted((Path(env["ProgramFiles"]) / "R").glob("R-*/bin/Rscript.exe"), reverse=True)
                rscript = next(iter(versions), None)
        if not rscript:
            parser.exit(1, "Rscript is required to draw figures.\n")
        subprocess.run([rscript, "--vanilla", "--encoding=UTF-8", SCRIPTS / "figures.R"],
                       env=env, check=True)


if __name__ == "__main__":
    main()
