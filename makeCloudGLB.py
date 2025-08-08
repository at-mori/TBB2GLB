import bz2
import numpy as np
import xarray as xr
import trimesh
import os
import urllib.request
import geopandas as gpd
from shapely.geometry import box, LineString, MultiLineString
from math import radians, sin, cos, tan, pi, log
from trimesh.visual.color import to_rgba
from trimesh.visual.material import SimpleMaterial

def tbb_to_color(tbb_value):
    C2K = 273.15
    if np.isnan(tbb_value):
        return [0.0, 0.0, 0.0, 1.0]  # NaNは黒で示す。
    elif tbb_value <=  0 + C2K:
        return [1.0, 1.0, 1.0, 1.0]
    elif tbb_value <=  5 + C2K:
        return [0.8, 0.8, 0.8, 1.0]
    elif tbb_value <= 10 + C2K:
        return [0.6, 0.6, 0.6, 1.0]
    elif tbb_value <= 15 + C2K:
        return [0.4, 0.4, 0.4, 1.0]
    else:
        return [0.2, 0.2, 0.2, 1.0]
    
def loadFile(year,month,day,hour,minutes):
    server = "ftp://hmwr829gr.cr.chiba-u.ac.jp/gridded/FD/V20190123"
    yyyymm = str(year)+str(month).zfill(2)
    dd     = str(day).zfill(2)
    hhnn   = str(hour).zfill(2)+str(minutes).zfill(2)
    band   = "tir"
    ch     = "01"
    ext    = "tbb.fld.4km.bin.bz2"
    fname  = f"{yyyymm}{dd}{hhnn}.{band}.{ch}.{ext}"
    url    = f"{server}/{yyyymm}/4KM/{yyyymm}{dd}/{fname}"

    if not os.path.exists(fname):
        urllib.request.urlretrieve(url, fname)
        print(f"✅ ダウンロードしました : {fname}")
    else:
        print(f"✅ 既にあります。 : {fname}")
        
    return fname

def load_tbb(filename, west=116, east=146, south=27, north=47, downsample=4):
    """圧縮されたTBBファイルを読み込み、指定範囲を抽出"""
    """3000 x 3000, 0.04 degree (4 km 相当)"""
    n = 3000  # global size
    lon_min, lon_max = 85, 205
    lat_min, lat_max = -60, 60
    resolution = 0.04
    lon = np.arange(lon_min, lon_max, resolution)
    lat = np.arange(lat_max, lat_min, -resolution)

    print("✅ Bz2 file start reading")

    with bz2.BZ2File(filename) as bz2file:
        dataTBB = np.frombuffer(bz2file.read(), dtype=">f4").reshape(n, n)

    print("✅ Bz2 file has read")
        
    xr_tbb = xr.DataArray(np.float32(dataTBB), name="tbb",
                          coords={'lat': ('lat', lat), 'lon': ('lon', lon)},
                          dims=['lat', 'lon'])

    sub = xr_tbb.sel(lat=slice(north, south), lon=slice(west, east))

    if downsample > 1:
        sub = sub.isel(
            lat=slice(0, None, downsample),
            lon=slice(0, None, downsample)
        )

    return sub

def tbb_to_height(tbb, surface_temp=300.0, lapse_rate=6.5):
    """等価黒体温度TBBを仮定した温度減率で高度[km]に変換"""
    return (surface_temp - tbb) / lapse_rate

def generate_cloud_and_coastline_glb(xr_tbb, out_path="combined_model.glb", threshold_diff=0.5):
    """雲と海岸線を統合してGLB出力"""

    heights = tbb_to_height(xr_tbb)
    max_tbb = np.nanmax(xr_tbb.values)
    base_height = tbb_to_height(max_tbb)
    dz = (heights - base_height).values

    lat_arr = xr_tbb['lat'].values
    lon_arr = xr_tbb['lon'].values
    ny, nx = dz.shape

    # === 正角円錐図法 ===
    lat_min, lat_max = lat_arr.min(), lat_arr.max()
    lon_min, lon_max = lon_arr.min(), lon_arr.max()
    lat1 = radians(lat_min + (lat_max - lat_min) * 3 / 10)
    lat2 = radians(lat_min + (lat_max - lat_min) * 7 / 10)
    lat0 = radians((lat_min + lat_max) / 2)
    lon0 = radians((lon_min + lon_max) / 2)
    n = (sin(lat1) + sin(lat2)) / 2
    F = cos(lat1) * (tan(pi/4 + lat1/2) ** n) / n
    rho0 = F / (tan(pi/4 + lat0/2) ** n)
    R = 6371  # 地球半径 [km]

    def project(phi_deg, lambda_deg):
        phi = radians(phi_deg)
        lam = radians(lambda_deg)
        rho = F / (tan(pi/4 + phi/2) ** n)
        theta = n * (lam - lon0)
        x = R * rho * sin(theta)
        y = R * (rho0 - rho * cos(theta))
        return x, y

    # === 雲メッシュ ===
    cloud_vertices, cloud_faces, cloud_colors = [], [], []

    def is_valid(i, j):
        return not np.isnan(dz[i, j])

    def smooth(i1, j1, i2, j2):
        return abs(dz[i1, j1] - dz[i2, j2]) < threshold_diff

    for i in range(ny - 1):
        for j in range(nx - 1):
            if all([is_valid(i+a, j+b) for a in [0, 1] for b in [0, 1]]):
                if (smooth(i,   j, i+1, j)   and smooth(i, j,   i, j+1) and
                    smooth(i+1, j, i+1, j+1) and smooth(i, j+1, i+1, j+1)):

                    x0, y0 = project(lat_arr[i],     lon_arr[j])
                    x1, y1 = project(lat_arr[i],     lon_arr[j+1])
                    x2, y2 = project(lat_arr[i+1],   lon_arr[j])
                    x3, y3 = project(lat_arr[i+1],   lon_arr[j+1])

                    v0 = [x0, y0, dz[i, j]]
                    v1 = [x1, y1, dz[i, j+1]]
                    v2 = [x2, y2, dz[i+1, j]]
                    v3 = [x3, y3, dz[i+1, j+1]]

                    idx = len(cloud_vertices)
                    cloud_vertices.extend([v0, v1, v2, v3])
                    cloud_faces.append([idx, idx+1, idx+2])
                    cloud_faces.append([idx+1, idx+3, idx+2])

                    color0 = tbb_to_color(xr_tbb.values[i, j])
                    color1 = tbb_to_color(xr_tbb.values[i, j+1])
                    color2 = tbb_to_color(xr_tbb.values[i+1, j])
                    color3 = tbb_to_color(xr_tbb.values[i+1, j+1])
                    cloud_colors.extend([color0, color1, color2, color3])
                    
    cloud_vertices = np.array(cloud_vertices)[:, [0, 2, 1]]
    cloud_vertices[:, 2] *= -1  # 方向を反転
    cloud_faces = np.array(cloud_faces)[:, ::-1] # 裏表反転
                    
    # === 海岸線メッシュ ===
    gdf = gpd.read_file("ne_10m_coastline.geojson")
    clip_box = box(lon_min, lat_min, lon_max, lat_max)
    gdf = gdf.to_crs("EPSG:4326").clip(clip_box)
    gdf['geometry'] = gdf['geometry'].simplify(0.01, preserve_topology=True)

    coast_vertices, coast_faces, coast_colors = [], [], []
    cream = [1.0, 0.95, 0.7, 1.0]

    half_width = 2.5  # [km]
    half_thick = 0.5  # [km]
    # z_offset = -base_height
    z_offset = 1.5

    for geom in gdf.geometry:
        if isinstance(geom, (LineString, MultiLineString)):
            lines = [geom] if isinstance(geom, LineString) else geom.geoms
            for line in lines:
                coords = list(line.coords)
                projected = [project(lat, lon) for lon, lat in coords]
                for i in range(len(projected) - 1):
                    p1 = np.array(projected[i])
                    p2 = np.array(projected[i + 1])
                    vec = p2 - p1
                    norm = np.linalg.norm(vec)
                    if norm == 0:
                        continue
                    vec /= norm
                    normal = np.array([-vec[1], vec[0]])

                    base = []
                    for dw in [-half_width, half_width]:
                        for dh in [-half_thick, half_thick]:
                            offset = normal * dw
                            base.append([p1[0] + offset[0], p1[1] + offset[1], z_offset + dh])
                            base.append([p2[0] + offset[0], p2[1] + offset[1], z_offset + dh])

                    idx = len(cloud_vertices) + len(coast_vertices)
                    coast_vertices.extend(base)
                    coast_colors.extend([cream] * 8)

                    quads = [
                        (0, 1, 3, 2), (4, 5, 7, 6), (0, 2, 6, 4),
                        (1, 3, 7, 5), (2, 3, 7, 6), (0, 1, 5, 4)
                    ]
                    for quad in quads:
                        a, b, c, d = [idx + q for q in quad]
                        coast_faces.append([a, b, c])
                        coast_faces.append([a, c, d])

    coast_vertices = np.array(coast_vertices)[:, [0, 2, 1]]
    coast_vertices[:, 2] *= -1  # 方向を反転, 裏表は正しい。
    coast_faces = np.array(coast_faces)
                        
    # === 統合メッシュ出力 ===
    all_vertices = np.vstack([cloud_vertices, coast_vertices])
    all_faces    = np.vstack([cloud_faces, coast_faces])
    all_colors   = np.vstack([cloud_colors, coast_colors])

    mesh = trimesh.Trimesh(vertices=all_vertices, faces=all_faces,
                           vertex_colors=(all_colors * 255).astype(np.uint8),
                           process=False)
    mesh.export(out_path)

    # 範囲出力
    x_min, y_min, z_min = all_vertices.min(axis=0)
    x_max, y_max, z_max = all_vertices.max(axis=0)
    print(f"✅ 統合GLB出力完了: {out_path}")
    print(f"📐 X: {x_min:.2f} ～ {x_max:.2f}")
    print(f"📐 Y: {y_min:.2f} ～ {y_max:.2f}")
    print(f"📐 Z: {z_min:.2f} ～ {z_max:.2f}")

#################################################################################
from datetime import datetime, timedelta, timezone

def to_utc(year, month, day, hour, local_offset_hours=9):
    dt_local = datetime(year, month, day, hour)
    return dt_local - timedelta(hours=local_offset_hours)


REGIONS = {
    "japan": {
        "west": 120, "east": 148,
        "south": 24, "north": 46,
        "downsample": 2,
        "threshold_diff": 4
    },
    "southjapan": {
        "west": 116, "east": 146,
        "south": 20, "north": 40,
        "downsample": 2,
        "threshold_diff": 4
    },
    "hokkaido": {
        "west": 136, "east": 149,
        "south": 39, "north": 48,
        "downsample": 1,
        "threshold_diff": 3
    },
    "kantou": {
        "west": 132, "east": 144,
        "south": 32, "north": 40,
        "downsample": 1,
        "threshold_diff": 3
    },
}

def generate_model(region_name, year, month, day, hour, local_time=True):
    region = REGIONS[region_name]
    utc = (to_utc(year, month, day, hour) if local_time else
           datetime(year, month, day, hour))

    # ファイル読み込み
    tbb = load_tbb(
        loadFile(utc.year, utc.month, utc.day, utc.hour, 0),
        west=region["west"], east=region["east"],
        south=region["south"], north=region["north"],
        downsample=region["downsample"]
    )

    # 出力ファイル名
    timestamp = utc.strftime("%Y%m%d%H")
    out_name = f"CM_{region_name}_{timestamp}.glb"

    # モデル生成
    generate_cloud_and_coastline_glb(
        tbb, out_path=out_name, threshold_diff=region["threshold_diff"]
    )

#################################################################################
    
# 使用例（コマンドラインでの実行想定）
if __name__ == "__main__":

    import sys

    if len(sys.argv) != 3:
        print("Usage: python makeCloudGLB.py REGION YYYYMMDDHH")
        print("Example: python makeCloudGLB.py kantou 2024082800")
        sys.exit(1)

    region = sys.argv[1].lower()
    timestamp = sys.argv[2]

    if region not in REGIONS:
        print(f"Unknown region: {region}")
        print(f"Available regions: {', '.join(REGIONS.keys())}")
        sys.exit(1)

    if not (len(timestamp) == 10 and timestamp.isdigit()):
        print("Time must be in format YYYYMMDDHH (e.g. 2024082800)")
        sys.exit(1)

    year  = int(timestamp[0:4])
    month = int(timestamp[4:6])
    day   = int(timestamp[6:8])
    hour  = int(timestamp[8:10])

    generate_model(region, year, month, day, hour, local_time=False)
