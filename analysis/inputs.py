"""Read the course archive and the recorded station coordinates."""
import csv
import os
import shutil
import subprocess
from pathlib import Path
from zipfile import ZipFile

import pyogrio
import shapely
import xlrd
from pyproj import Transformer

PROJECT = Path(__file__).resolve().parents[1]
ARCHIVE = PROJECT / "实习-第一部分.zip"
RESULTS = PROJECT / "results"
WORK = PROJECT / "run"
SOURCE = WORK / "source"


def station_inputs():
    province = SOURCE / "province"
    province.mkdir(parents=True, exist_ok=True)
    with ZipFile(ARCHIVE, metadata_encoding="gbk") as archive:
        sheet_name = next(name for name in archive.namelist() if name.endswith(".xls"))
        book = xlrd.open_workbook(file_contents=archive.read(sheet_name))
        sheet = book.sheet_by_name("空间分析")
        for name in archive.namelist():
            if Path(name).stem == "Gansu":
                (province / Path(name).name).write_bytes(archive.read(name))

    with (RESULTS / "stations.csv").open(encoding="utf-8", newline="") as stream:
        saved = {int(row["id"]): row for row in csv.DictReader(stream)}
    project = Transformer.from_crs(4326, 2381, always_xy=True).transform
    stations = []
    for number in range(1, sheet.nrows):
        station_id, name, lon, lat, elev, rain, evap = sheet.row_values(number)
        record = saved[int(station_id)]
        x, y = float(record["x_m"]), float(record["y_m"])
        px, py = project(lon, lat)
        if name != record["name"] or abs(x - px) > 1 or abs(y - py) > 1:
            raise ValueError(f"Station coordinates disagree for {name}")
        stations.append(dict(id=int(station_id), name=name, lon=lon, lat=lat,
                             elev=elev, precip=rain, evap=evap, xy=[x, y]))
    _, _, geometries, _ = pyogrio.raw.read(province / "Gansu.shp")
    boundary = shapely.union_all(shapely.from_wkb(geometries))
    return stations, shapely.transform(boundary, project, interleaved=False)


def unpack_topic():
    """The land-cover raster and its county boundaries are inside the nested RAR."""
    target = SOURCE / "topic"
    target.mkdir(parents=True, exist_ok=True)
    exe = shutil.which("unrar") or shutil.which("7z")
    if not exe:
        for path in ("WinRAR/UnRAR.exe", "7-Zip/7z.exe"):
            candidate = Path(os.environ.get("ProgramFiles", "/opt")) / path
            if candidate.is_file():
                exe = str(candidate)
                break
    if not exe:
        raise SystemExit("Install UnRAR or 7-Zip to extract the course's thematic data.")
    with ZipFile(ARCHIVE, metadata_encoding="gbk") as archive:
        rar = next(name for name in archive.namelist() if name.endswith(".rar"))
        temp = target / "topic.rar"
        temp.write_bytes(archive.read(rar))
    try:
        patterns = ["Gansu-County.*", "Gansu-lucc-2015.*"]
        if Path(exe).stem.lower() == "unrar":
            args = [exe, "e", "-o+", "-inul", "-p-", str(temp), *patterns]
        else:
            args = [exe, "e", "-y", "-bso0", str(temp), *patterns]
        subprocess.run(args, cwd=target, check=True)
    finally:
        temp.unlink(missing_ok=True)
